export GADM, NUTS, makeregions, makeregions_nuts, makeoffshoreregions, saveregions, loadregions

struct GADM{T}
    parentregions::Vector{T}
    subregionnames::NTuple{N,T} where N
end
GADM(regionnames::T...) where T = GADM(T[], regionnames)
GADM(parentregions::Vector{T}, subregionnames::T...) where T = GADM(parentregions, subregionnames)

struct NUTS{T}
    subregionnames::NTuple{N,T} where N
end
NUTS(regionnames::T...) where T = NUTS(regionnames)

function saveregions(regionname, regiondefinitionarray; autocrop=true, uselandcovermask=true, bbox=[-90 -180; 90 180])
    datafolder = getconfig("datafolder")
    land = JLD.load(joinpath(datafolder, "landcover.jld"), "landcover")
    if !all(bbox .== [-90 -180; 90 180])
        autocrop = false         # ignore supplied autocrop option if user changed bbox
    end
    saveregions(regionname, regiondefinitionarray, land, autocrop, uselandcovermask, bbox)
end

function saveregions(regionname, regiondefinitionarray, landcover, autocrop, uselandcovermask, bbox)
    regions = makeregions(regiondefinitionarray)
    if uselandcovermask
        regions = regions .* (landcover.>0)
    end
    if autocrop
        # get indexes of the bounding box containing onshore region data with 3 degrees of padding
        lonrange, latrange = getbboxranges(regions, round(Int, 3/0.01))
    else
        latrange, lonrange = bbox2ranges(roundbbox(bbox,100), 100)      # TO DO: remove hardcoded raster density
    end
    regions = regions[lonrange, latrange]
    offshoreregions = makeoffshoreregions(regions, landcover[lonrange, latrange])
    regionlist = Symbol.(regiondefinitionarray[:,1])

    println("\nSaving regions and offshoreregions...")
    datafolder = getconfig("datafolder")
    JLD.save(joinpath(datafolder, "regions_$regionname.jld"), "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

function saveregions_global_gadm0()
    println("Creating a global region file to identify countries and land areas later (5-15 minutes)...\n")
    datafolder = getconfig("datafolder")
    g = readdlm(joinpath(datafolder, "gadmfields.csv"), ',', skipstart=1)
    gadm0 = unique(string.(g[:,2]))
    regiondefinitionarray = [gadm0 GADM.(gadm0)]
    # This map is used to identify country by pixel, so don't mask by landcover (i.e. set region=0 in lakes).
    saveregions("Global_GADM0", regiondefinitionarray, autocrop=false, uselandcovermask=false)
    println("\nGlobal region file saved.")
end

function loadregions(regionname)
    datafolder = getconfig("datafolder")
    jldopen(joinpath(datafolder, "regions_$regionname.jld"), "r") do file
        return read(file, "regions"), read(file, "offshoreregions"), read(file, "regionlist"),
                    read(file, "lonrange"), read(file, "latrange")
    end
end

function makeregions(regiondefinitionarray)
    regionnames, nutsdef, gadmdef = splitregiondefinitions(regiondefinitionarray)
    region = zeros(Int16, (36000,18000))    # hard code size for now
    if !all(isempty.(nutsdef))
        nuts, subregionnames = read_nuts()
        makeregions_nuts!(region, nuts, subregionnames, nutsdef)
    end
    if !all(isempty.(gadmdef))
        gadm, subregionnames = read_gadm()
        makeregions_gadm!(region, gadm, subregionnames, gadmdef)
    end
    return region
end

function regions2matlab(gisregion)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    datafolder = getconfig("datafolder")
    matopen(joinpath(datafolder, "regions_$gisregion.mat"), "w", compress=true) do file
        write(file, "regions", regions)
        write(file, "offshoreregions", offshoreregions)
        write(file, "regionlist", string.(regionlist))
        write(file, "lonrange", collect(lonrange))
        write(file, "latrange", collect(latrange))
    end
end

function splitregiondefinitions(regiondefinitionarray)
    regionnames = regiondefinitionarray[:,1]
    regiondefinitions = [regdef isa Tuple ? regdef : (regdef,) for regdef in regiondefinitionarray[:,2]]
    nutsdef = [Tuple(rd for rd in regdef if rd isa NUTS) for regdef in regiondefinitions]
    gadmdef = [Tuple(rd for rd in regdef if rd isa GADM) for regdef in regiondefinitions]
    return regionnames, nutsdef, gadmdef
end

function makeregions_gadm!(region, gadm, subregionnames, regiondefinitions)
    println("Making region index matrix...")
    regionlookup = build_inverseregionlookup(regiondefinitions)
    updateprogress = Progress(prod(size(region)), 1)
    for (i, g_uid) in enumerate(gadm)
        next!(updateprogress)
        g_uid == 0 && continue
        reg0, reg1, reg2 = subregionnames[g_uid,:]
        regid = lookup_regionnames(regionlookup, reg0, reg1, reg2)
        if regid > 0
            region[i] = regid
        end
    end
end

function makeregions_nuts!(region, nuts, subregionnames, regiondefinitions)
    println("Making region index matrix...")
    regionlookup = Dict(r => i for (i,tuptup) in enumerate(regiondefinitions)
                                    for ntup in tuptup for r in ntup.subregionnames)
    updateprogress = Progress(prod(size(region)), 1)
    for (i, g_uid) in enumerate(nuts)
        next!(updateprogress)
        g_uid == 0 && continue
        reg = subregionnames[g_uid]
        while length(reg) >= 2
            regid = get(regionlookup, reg, 0)
            if regid > 0
                region[i] = regid
                break
            end
            reg = reg[1:end-1]
        end
    end
end

function lookup_regionnames(regionlookup, reg0, reg1, reg2)
    v = get(regionlookup, (reg0, "*", "*"), 0)
    v > 0 && return v
    v = get(regionlookup, (reg0, reg1, "*"), 0)
    v > 0 && return v
    return get(regionlookup, (reg0, reg1, reg2), 0)
end

function build_inverseregionlookup(regiondefinitions)
    d = Dict{Tuple{String,String,String}, Int}()
    for reg = 1:length(regiondefinitions)
        for regdef in regiondefinitions[reg]
            parentregions, subregionnames = regdef.parentregions, regdef.subregionnames
            regions = ["*", "*", "*"]
            regions[1:length(parentregions)] = parentregions
            for s in subregionnames
                regions[length(parentregions)+1] = s
                d[regions...] = reg
            end
        end
    end
    return d
end

# Find the closest region pixel for each ocean pixel and major lake.
# Even VERY far offshore pixels will be allocated to whatever region is nearest, but
# those areas still won't be available for offshore wind power because of the
# requirement to be close enough to the electricity grid (or rather the grid proxy).
function makeoffshoreregions(regions, landcover)
    @time lakes = gridsplit(regions, majorlakes, Bool)  # = majorlakes(regions), but chunked calculation
    println("\nMaking offshore region index matrix...")
    closest_region = regions[feature_transform(regions.>0)]
    return closest_region .* ((regions .== 0) .| lakes) .* (landcover.==0)
end

# Use ImageSegmentation.jl to identify large lakes (large enough for offshore wind).
# This will incorrectly trigger for some very long and wide rivers, but that's not a
# big problem since near-shore areas are not allowed for offshore wind. 
function majorlakes(regions)
    println("Identifying major lakes (>1000 km2)...")
    println("...segmenting regions...")
    seg = fast_scanning(regions.>0, 0.1)
    println("...removing small segments and land areas...")
    large_water_segments = prune_segments(seg, i -> (segment_pixel_count(seg,i)<1000 || segment_mean(seg,i)>0.8),
                                            (i,j) -> -segment_pixel_count(seg,j))
    lakes = labels_map(large_water_segments)
    println("...removing the largest land segment...")
    pixelcount = countmap(lakes, alg=:dict)
    # the most common index will (almost certainly) be the sole non-lake segment
    mostcommonindex = sort(collect(pixelcount), by=x->x[2], rev=true)[1][1]
    println("...lake identification complete.")    
    return lakes .!= mostcommonindex
end
