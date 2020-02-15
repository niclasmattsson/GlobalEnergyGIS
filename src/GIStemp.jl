tempoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :mindist => 3.3                     # minimum distance between population centers [in ERA5 pixels]
)

mutable struct TempOptions
    gisregion               ::String
    scenarioyear            ::String
    era_year                ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    mindist                 ::Float64
end

TempOptions() = TempOptions("","",0,0,0,0)

function TempOptions(d::Dict{Symbol,Any})
    options = TempOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function GIStemp(; optionlist...)

    options = TempOptions(merge(tempoptions(), optionlist))
    @unpack gisregion, era_year, mindist = options

    regions, offshoreregions, regionlist, lonrange, latrange, pop, meantemp, temp =
                read_temperature_datasets(options)

    erapop = rescale_population_to_ERA5_res(options, regions, offshoreregions, regionlist, lonrange, latrange, pop, temp)

    numreg = length(regionlist)
    popcenters = [findpopcenters(erapop[i,:,:], 3, mindist) for i = 1:numreg]

    hours = DateTime(era_year, 1, 1, 0) : Hour(1) : DateTime(era_year, 12, 31, 23)
    numhours = length(hours)
    numhours != size(temp,1) && error("Inconsistent number of hours.")

    datafolder = getconfig("datafolder")
    syntheticdemanddata = joinpath(datafolder, "syntheticdemand", "data", "julia")
    mkpath(syntheticdemanddata)

    for (r, reg) in enumerate(regionlist)
        filename = joinpath(syntheticdemanddata, "tempdata_$(era_year)_$(gisregion)_$(r)_$reg.csv")
        open(filename, "w") do file
            write(file, "time,temp1,temp2,temp3\n")
            for h = 1:numhours
                write(file, "$(Date(hours[h])) $(Time(hours[h]))")
                for i = 1:3
                    write(file, ",", string(temp[h, popcenters[r][i]...]))
                end
                write(file, "\n")
            end
        end
    end
end

function read_temperature_datasets(options)
    @unpack res, erares, gisregion, scenarioyear, era_year = options

    println("\nReading auxiliary datasets...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    datafolder = getconfig("datafolder")
    pop = JLD.load(joinpath(datafolder, "population_$scenarioyear.jld"), "population")[lonrange,latrange]

    println("Reading ERA5 temperature dataset...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    filename = joinpath(datafolder, "era5temp$era_year.h5")
    @time meantemp, temp = h5open(filename, "r") do file
        if length(eralonranges) == 1
            file["meantemp"][eralonranges[1], eralatrange],
                file["temp"][:,eralonranges[1], eralatrange]
        else
            [file["meantemp"][eralonranges[1], eralatrange]; file["meantemp"][eralonranges[2], eralatrange]],
                [file["temp"][:, eralonranges[1], eralatrange] file["temp"][:, eralonranges[2], eralatrange]]
        end
    end
    return regions, offshoreregions, regionlist, lonrange, latrange, pop, meantemp, temp
end


function rescale_population_to_ERA5_res(options, regions, offshoreregions, regionlist, lonrange, latrange, pop, temp)
    println("Calculating population centers at ERA5 resolution...")
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(options, lonrange, latrange)

    @unpack era_year, res, erares = options

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(temp)
    firsttime = DateTime(era_year, 1, 1)

    erapop = zeros(numreg, nlons, nlats)

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                if reg > 0 && reg != NOREGION
                    erapop[reg,i,j] += pop[r,c]
                end
            end
        end
        next!(updateprogress)
    end

    return erapop
end

# Find the (row,col) coordinates of the n largest population cells which are at least mindist cells apart from each other.
function findpopcenters(pop, n, mindist)
    n = min(n, sum(pop .> 0))   # make sure there are enough non-zero cells
    p = copy(pop)
    coords = CartesianIndex{2}[]
    for i = 1:n
        val, index = findmax(p)
        if val == 0   # if the largest cell has zero pop then we were too ambitious - redo with lower distance
            # println("lowering distance")
            return findpopcenters(pop, n, mindist/2)
        end
        push!(coords, index)
        fillcircle!(p, Tuple(index), mindist, 0)
    end
    # display([dist(coords[i], coords[j]) for i = 1:length(coords)-1, j = 2:length(coords) if j > i])
    return Tuple.(coords)
end

dist(a::CartesianIndex, b::CartesianIndex) = sqrt(sum((Tuple(a) .- Tuple(b)).^2))

function fillcircle!(a::AbstractMatrix, center::Tuple{Int,Int}, radius::Real, fillvalue)
    nrows, ncols = size(a)
    rowrange = max(1, ceil(Int, center[1] - radius)) : min(nrows, floor(Int, center[1] + radius))
    colrange = max(1, ceil(Int, center[2] - radius)) : min(nrows, floor(Int, center[2] + radius))
    for c in colrange, r in rowrange
        if (r - center[1])^2 + (c - center[2])^2 <= radius^2
            a[r,c] = fillvalue
        end
    end
end
