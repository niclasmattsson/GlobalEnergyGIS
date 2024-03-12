export GADM, NUTS, makeregions, makeregions_nuts, makeoffshoreregions, saveregions, loadregions,
        saveregions_global, subregions, all_gadm_subregions

abstract type RegionType end

struct GADM{T} <: RegionType
    parentregions::Vector{T}
    subregionnames::NTuple{N,T} where N
end
GADM(regionnames::T...) where T = GADM(T[], regionnames)
GADM(parentregions::Vector{T}, subregionnames::T...) where T = GADM(parentregions, subregionnames)

struct NUTS{T} <: RegionType
    subregionnames::NTuple{N,T} where N
end
NUTS(regionnames::T...) where T = NUTS(regionnames)

const NOREGION = typemax(Int16)

function saveregions(regionname, subregionnames, regions::Matrix{Int32})
    land = JLD.load(in_datafolder("landcover.jld"), "landcover")
    saveregions(regionname, subregionnames, regions, :GADM, land, true, Int[;;])
end

function saveregions(regionname, regiondefinitionarray; autocrop=true, bbox=[-90 -180; 90 180])
    land = JLD.load(in_datafolder("landcover.jld"), "landcover")
    if !all(bbox .== [-90 -180; 90 180])
        autocrop = false         # ignore supplied autocrop option if user changed bbox
    end
    regions, regiontype = makeregions(regiondefinitionarray; allowmixed=(regionname=="Europe_background"))
    saveregions(regionname, regiondefinitionarray, regions, regiontype, land, autocrop, bbox)
end

function saveregions(regionname, regiondefinitionarray, regions, regiontype, landcover, autocrop, bbox)
    if autocrop
        # get indexes of the bounding box containing onshore region data with 6% of padding
        lonrange, latrange = getbboxranges(regions)
        padding = round(Int, maximum(size(regions[lonrange,latrange])) * 0.06)
        lonrange, latrange = getbboxranges(regions, padding)
    else
        latrange, lonrange = bbox2ranges(roundbbox(bbox,100), 100)          # TO DO: remove hardcoded raster density
    end
    landcover = landcover[lonrange, latrange]
    regions = regions[lonrange, latrange]

    if regionname != "Global_GADM0" && regionname != "Europe_background"
        if regiontype == :NUTS
            println("\nNUTS region definitions detected (using Europe_background region file)...")
            europeregions = loadregions("Europe_background")[1][lonrange, latrange]
            regions[(regions.==0) .& (europeregions.>0)] .= NOREGION
        elseif regiontype == :GADM
            println("\nGADM region definitions detected (using Global_GADM0 region file)...")
            globalregions = loadregions("Global_GADM0")[1][lonrange, latrange]
            regions[(regions.==0) .& (globalregions.>0)] .= NOREGION
        end
    end

    # Find the closest region pixel for all non-region pixels (land and ocean)
    println("\nAllocate non-region pixels to the nearest region (for offshore wind)...")
    territory = regions[feature_transform(regions.>0)]

    # Allocate ocean and lake pixels to the region with the closest land region.
    # Even VERY far offshore pixels will be allocated to whatever region is nearest, but
    # those areas still won't be available for offshore wind power because of the
    # requirement to be close enough to the electricity grid (or rather the grid proxy).
    offshoreregions = territory .* (landcover .== 0)

    if regionname != "Global_GADM0" && regionname != "Europe_background"
        # Allocate land pixels with region==0 to the closest land region.
        # This ensures that the regions dataset is pixel-compatible with the landcover dataset.
        regions = territory .* (landcover .> 0)
    end

    println("\nSaving regions and offshoreregions...")
    regionlist = Symbol.(regiondefinitionarray[:,1])

    JLD.save(in_datafolder("regions_$regionname.jld"), "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

function saveregions_global(; args...)
    println("\nCreating a global GADM region file to identify countries and land areas later...")
    g = readdlm(in_datafolder("gadmfields.csv"), ',', skipstart=1)
    gadm0 = unique(string.(g[:,2]))
    regiondefinitionarray = [gadm0 GADM.(gadm0)]
    saveregions("Global_GADM0", regiondefinitionarray; args..., autocrop=false)
    println("Global GADM region file saved.")

    println("\nCreating a 'background' NUTS region file to identify non-European land areas later...")
    regiondefinitionarray = [NUTS_Europe; non_NUTS_Europe]
    saveregions("Europe_background", regiondefinitionarray; args..., autocrop=false)
    println("\nEurope_background region file saved.")

    println("\nCreating a region file for the 44 countries with synthetic demand training data...")
    regiondefinitionarray = [syntheticdemandregions GADM.(syntheticdemandregions)]
    saveregions("SyntheticDemandRegions", regiondefinitionarray; args...)
    println("\nSyntheticDemandRegions file saved.")
end

function loadregions(regionname)
    jldopen(in_datafolder("regions_$regionname.jld"), "r") do file
        return read(file, "regions"), read(file, "offshoreregions"), read(file, "regionlist"),
                    read(file, "lonrange"), read(file, "latrange")
    end
end

function makeregions(regiondefinitionarray; allowmixed=false)
    regionnames, nutsdef, gadmdef = splitregiondefinitions(regiondefinitionarray)
    use_nuts, use_gadm = !all(isempty.(nutsdef)), !all(isempty.(gadmdef))
    regiontype =  (use_gadm && !use_nuts) ? :GADM  :
                  (use_nuts && !use_gadm) ? :NUTS  :
                  (use_nuts && use_gadm) ? :MIXED  :  :WEIRD
    !allowmixed && regiontype==:MIXED && error("Sorry, mixed NUTS & GADM definitions are not supported yet.")
    region = zeros(Int16, (36000,18000))    # hard code size for now
    if use_nuts
        nuts, subregionnames = read_nuts()
        makeregions_nuts!(region, nuts, subregionnames, nutsdef)
    end
    if use_gadm
        gadm, subregionnames = read_gadm()
        makeregions_gadm!(region, gadm, subregionnames, gadmdef)
    end
    return region, regiontype
end

function regions2matlab(gisregion)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    matopen(in_datafolder("regions_$gisregion.mat"), "w", compress=true) do file
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
            (gadm_uid == 0 || gadm_uid == 78413 || region[r,c] > 0) && continue    # ignore Caspian Sea (weirdly classified as a region in GADM)
            reg0, reg1, reg2 = subregionnames[gadm_uid,:]
            regid = lookup_regionnames(regionlookup, reg0, reg1, reg2)
            if regid > 0
                region[r,c] = regid
            end
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
            (nuts_id == 0 || region[r,c] > 0) && continue
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

function getsubregions(regtype::Type{GADM}, regionnames)
    gadm = readdlm(in_datafolder("gadmfields.csv"), ',', skipstart=1)[:, 2:4]
    for reg in regionnames
        gadm = gadm[gadm[:,1].==reg, 2:end]
    end
    return sort(unique(gadm[:,1]))
end

function getsubregions(regtype::Type{NUTS}, regionname)
    length(regionname) > 1 && error("Give only one string argument to match beginning of NUTS region names.")
    nuts = readdlm(in_datafolder("nutsfields.csv"), ',', skipstart=1)[:, 3]
    if isempty(regionname)
        return sort(unique(first.(nuts, 2)))
    else
        return sort(nuts[startswith.(nuts, regionname)])
    end
end

function subregions(regtype::Type{T} where T <: RegionType, regionnames::String...)
    reglist = join(regionnames, ", ")
    selected = string.(getsubregions(regtype, regionnames))
    selectedlist = join(selected, ", ")
    # isempty(regionnames) && println("Showing top level $regtype regions:")
    # println("$regtype($reglist): $selectedlist")
    return selected
end

function all_gadm_subregions(country::AbstractString, level::Int)
    level == 1 && return [country GADM(country)]
    return all_gadm_subregions(country, subregions(GADM, country), level)
end

function all_gadm_subregions(country::AbstractString, reg2::AbstractString, level::Int)
    reg2 == "" && return [country GADM(country)]
    level == 2 && return [reg2 GADM([country], reg2)]
    return vcat(all_gadm_subregions.(country, reg2, subregions(GADM, country, reg2))...)
end

function all_gadm_subregions(country::AbstractString, reg2::AbstractString, reg3::AbstractString)
    reg3 == "" && return [reg2 GADM([country], reg2)]
    return [reg3 GADM([country, reg2], reg3)]
end

all_gadm_subregions(countries::AbstractArray, level::Int) =
        vcat(all_gadm_subregions.(countries, level)...)

all_gadm_subregions(country::AbstractString, regions2::AbstractArray, level::Int) =
        vcat(all_gadm_subregions.(country, regions2, level)...)
