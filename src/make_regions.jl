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

const NOREGION = typemax(Int16)

function saveregions(regionname, regiondefinitionarray; autocrop=true, bbox=[-90 -180; 90 180])
    datafolder = getconfig("datafolder")
    land = JLD.load(joinpath(datafolder, "landcover.jld"), "landcover")
    if !all(bbox .== [-90 -180; 90 180])
        autocrop = false         # ignore supplied autocrop option if user changed bbox
    end
    saveregions(regionname, regiondefinitionarray, land, autocrop, bbox)
end

function saveregions(regionname, regiondefinitionarray, landcover, autocrop, bbox)
    regions = makeregions(regiondefinitionarray)
    if autocrop
        # get indexes of the bounding box containing onshore region data with 3 degrees of padding
        lonrange, latrange = getbboxranges(regions, round(Int, 3/0.01))
    else
        latrange, lonrange = bbox2ranges(roundbbox(bbox,100), 100)          # TO DO: remove hardcoded raster density
    end
    landcover = landcover[lonrange, latrange]
    regions = regions[lonrange, latrange]

    if regionname != "Global_GADM0" && regionname != "Global_NUTS0"
        regiontype = determineregiontype(regiondefinitionarray)
        regiontype == :MIXED && error("Mixed GADM and NUTS region definitions not supported yet.")
        level0_regions, _, level0_regionlist = loadregions("Global_$(regiontype)0")
        level0_regions = level0_regions[lonrange, latrange]

        current_level0_regions = unique(level0_regions[regions .> 0])
        # return countmap(level0_regions[regions .> 0]), level0_regionlist

        println("\nMark non-region land areas in $(regiontype) dataset...")
        # set non-region land areas to NOREGION (a large positive integer), but skip pixels that have GADM region codes in active regions
        # (to avoid misalignment problems between GADM and NUTS)
        for (i, r) in enumerate(regions)
            if r == 0 && level0_regions[i] > 0 && !in(level0_regions[i], current_level0_regions)
                regions[i] = NOREGION
            end
        end
        # regions[(regions .== 0) .& (gadm_level0_regions .> 0)] .= NOREGION
    end

    # Find the closest region pixel for all non-region pixels (land and ocean)
    println("\nAllocate non-region pixels to the nearest region (for offshore wind and to fix misalignment of datasets)...")
    territory = regions[feature_transform(regions.>0)]

    # Use this to make a regions dataset that is pixel-compatible with the landcover dataset.
    regions = territory .* (landcover.>0)

    # Allocate ocean and (major) lake pixels to the region with the closest land region.
    # Even VERY far offshore pixels will be allocated to whatever region is nearest, but
    # those areas still won't be available for offshore wind power because of the
    # requirement to be close enough to the electricity grid (or rather the grid proxy).
    @time waterareas = gridsplit(landcover, identify_large_water, Bool)     # = identify_large_water(landcover), but chunked calculation
    offshoreregions = territory .* waterareas

    println("\nSaving regions and offshoreregions...")
    datafolder = getconfig("datafolder")
    regionlist = Symbol.(regiondefinitionarray[:,1])

    JLD.save(joinpath(datafolder, "regions_$regionname.jld"), "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

function determineregiontype(regiondefinitionarray)
    tuples = isa.(regiondefinitionarray[:,2], Tuple)
    list = regiondefinitionarray[.!tuples, 2]
    for tup in regiondefinitionarray[tuples, 2], t in tup
        push!(list, t)
    end
    gadm = sum(isa.(list, GADM{String}))
    nuts = sum(isa.(list, NUTS{String}))
    if gadm > 0 && nuts == 0
        return :GADM
    elseif nuts > 0 && gadm == 0
        return :NUTS
    elseif nuts > 0 && gadm > 0
        return :MIXED
    else
        error("Region definitions look fishy.")
    end
end

# Use ImageSegmentation.jl to identify large bodies of water (large enough for offshore wind).
# This will incorrectly trigger for some very long and wide rivers, but that's not a
# big problem since near-shore areas are not allowed for offshore wind. 
function identify_large_water(landcover)
    println("Identifying oceans and major lakes (>1000 km2)...")
    println("...segmenting...")
    seg = fast_scanning(landcover.>0, 0.1)
    println("...removing small segments and land areas...")
    large_segments = prune_segments(seg, i -> segment_pixel_count(seg,i) < 1000, (i,j) -> -segment_pixel_count(seg,j))  # merge small segments with their largest neighbors
    println("...classify as land or water...")
    sm = segment_mean(large_segments)
    waterareas = [(sm[i] < 0.5) for i in labels_map(large_segments)]
    println("...identification complete.")    
    return waterareas
end

function saveregions_global()
    datafolder = getconfig("datafolder")
    println("Creating a global GADM region file to identify countries and land areas later...\n")
    g = readdlm(joinpath(datafolder, "gadmfields.csv"), ',', skipstart=1)
    gadm0 = unique(string.(g[:,2]))
    regiondefinitionarray = [gadm0 GADM.(gadm0)]
    # This map is used to identify country by pixel, so don't mask by landcover (i.e. set region=0 in lakes).
    saveregions("Global_GADM0", regiondefinitionarray, autocrop=false)
    println("\nGlobal GADM region file saved.")

    println("\nCreating a global NUTS region file to identify countries and land areas later...\n")
    n = readdlm(joinpath(datafolder, "nutsfields.csv"), ',', skipstart=1)
    nuts0 = unique(string.(n[:,4]))
    regiondefinitionarray = [nuts0 NUTS.(nuts0)]
    # This map is used to identify country by pixel, so don't mask by landcover (i.e. set region=0 in lakes).
    saveregions("Global_NUTS0", regiondefinitionarray, autocrop=false)
    println("\nGlobal NUTS region file saved.")
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
