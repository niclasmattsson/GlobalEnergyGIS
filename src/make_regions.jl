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

function saveregions(regionname, regiondefinitionarray; autocrop=true, bbox=[-90 -180; 90 180], maxlakesize=1000, prunesize=100)
    datafolder = getconfig("datafolder")
    land = JLD.load(joinpath(datafolder, "landcover.jld"), "landcover")
    if !all(bbox .== [-90 -180; 90 180])
        autocrop = false         # ignore supplied autocrop option if user changed bbox
    end
    saveregions(regionname, regiondefinitionarray, land, autocrop, bbox, maxlakesize, prunesize)
end

function saveregions(regionname, regiondefinitionarray, landcover, autocrop, bbox, maxlakesize, prunesize)
    regions = makeregions(regiondefinitionarray)
    if autocrop
        # get indexes of the bounding box containing onshore region data with 3 degrees of padding
        lonrange, latrange = getbboxranges(regions, round(Int, 3/0.01))
    else
        latrange, lonrange = bbox2ranges(roundbbox(bbox,100), 100)          # TO DO: remove hardcoded raster density
    end
    landcover = landcover[lonrange, latrange]
    regions = regions[lonrange, latrange]

    # Find the closest region pixel for all non-region pixels (land and ocean)
    println("\nAllocate non-region pixels to the nearest region (for offshore wind)...")
    println("(This also fixes pixel-scale misalignments of GADM, NUTS and land cover datasets.)")
    territory = regions[feature_transform(regions.>0)]

    # Allocate land pixels with region==0 to the closest land region.
    # This ensures that the regions dataset is pixel-compatible with the landcover dataset.
    regions = territory .* (landcover .> 0)

    # Allocate ocean and lake pixels to the region with the closest land region.
    # Even VERY far offshore pixels will be allocated to whatever region is nearest, but
    # those areas still won't be available for offshore wind power because of the
    # requirement to be close enough to the electricity grid (or rather the grid proxy).
    offshoreregions = territory .* (landcover .== 0)

    println("\nSaving regions and offshoreregions...")
    datafolder = getconfig("datafolder")
    regionlist = Symbol.(regiondefinitionarray[:,1])

    JLD.save(joinpath(datafolder, "regions_$regionname.jld"), "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

# Use ImageSegmentation.jl to prune areas with small contiguous segments. 
function prune_areas(dataset, minpixelcount)
    println("...pruning segments smaller than $(minpixelcount) pixels...")
    seg = fast_scanning(dataset, 0.1)
    # merge small segments with their largest neighbors
    large_segments = prune_segments(seg, i -> segment_pixel_count(seg,i) < minpixelcount, (i,j) -> -segment_pixel_count(seg,j))
    indices = labels_map(large_segments)
    println("...rebuilding pruned dataset...")
    sm = segment_mean(large_segments)
    areas = [(sm[i] > 0.5) for i in indices]
    return areas
end

function saveregions_global(; args...)
    datafolder = getconfig("datafolder")
    println("Creating a global GADM region file to identify countries and land areas later...\n")
    g = readdlm(joinpath(datafolder, "gadmfields.csv"), ',', skipstart=1)
    gadm0 = unique(string.(g[:,2]))
    regiondefinitionarray = [gadm0 GADM.(gadm0)]
    # This map is used to identify country by pixel, so don't mask by landcover (i.e. set region=0 in lakes).
    saveregions("Global_GADM0", regiondefinitionarray; args..., autocrop=false)
    println("\nGlobal GADM region file saved.")
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
    rows, cols = size(region)
    updateprogress = Progress(cols, 1)
    for c in randperm(cols)
        for r = 1:rows
            gadm_uid = gadm[r,c]
            (gadm_uid == 0 || gadm_uid == 78413) && continue    # ignore Caspian Sea (weirdly classified as a region in GADM)
            reg0, reg1, reg2 = subregionnames[gadm_uid,:]
            regid = lookup_regionnames(regionlookup, reg0, reg1, reg2)
            region[r,c] = (regid > 0) ? regid : NOREGION
        end
        next!(updateprogress)
    end
end

function makeregions_nuts!(region, nuts, subregionnames, regiondefinitions)
    println("Making region index matrix...")
    regionlookup = Dict(r => i for (i,tuptup) in enumerate(regiondefinitions)
                                    for ntup in tuptup for r in ntup.subregionnames)
    rows, cols = size(region)
    updateprogress = Progress(cols, 1)
    for c in randperm(cols)
        for r = 1:rows
            nuts_id = nuts[r,c]
            nuts_id == 0 && continue
            reg = subregionnames[nuts_id]
            while length(reg) >= 2
                regid = get(regionlookup, reg, 0)
                if regid > 0
                    region[r,c] = regid
                    break
                end
                reg = reg[1:end-1]
            end
        end
        next!(updateprogress)
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
