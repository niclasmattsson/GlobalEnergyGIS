function makedistances(gisregion; scenarioyear="ssp2_2050", res=0.01)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    numreg = length(regionlist)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    datafolder = getconfig("datafolder")
    outputfolder = joinpath(datafolder, "output")
    mkpath(outputfolder)

    pop = JLD.load(joinpath(datafolder, "population_$scenarioyear.jld"), "population")[lonrange,latrange]
    popdens = pop ./ cellarea'

    geocenters, popcenters = regioncenters(regions, numreg, popdens, lonrange, latrange, res)
    println("\nCalculating distances between centers...")
    distances = [greatcircledistance(popcenters[r1], popcenters[r2]) for r1 = 1:numreg, r2 = 1:numreg]
    connected, connectedoffshore = connectedregions(regions, offshoreregions, numreg)
    regionpop = [sum(pop.==r) for r = 1:numreg]

    matopen(joinpath(outputfolder, "distances_$gisregion.mat"), "w") do file
        write(file, "distances", distances)
        write(file, "connected", connected)
        write(file, "connectedoffshore", connectedoffshore)
        write(file, "regionlist", string.(regionlist))
        # write(file, "population", pop)
        # write(file, "demand", demand)
        write(file, "regioncenters_lat", [latlon[1] for latlon in popcenters])
        write(file, "regioncenters_lon", [latlon[2] for latlon in popcenters])
    end

    # Also save the region file in Matlab format
    # regions2matlab(gisregion)
    nothing
end

# returns great circle distance in km between points given as (lat,lon) tuples (in degrees).
greatcircledistance(point1::Tuple, point2::Tuple) = haversine(point1, point2, 6371.0)

# evaluate a vector at a float index usinglinear interpolation between elements
function lerp(vect::AbstractArray, floatindex::Float64)
    floatindex >= length(vect) && return vect[end]
    f, index = modf(floatindex)
    index = clamp(Int(index), 1, length(vect)-1)
    return vect[index]*(1-f) + vect[index+1]*f
end
# lerp(vect::AbstractArray, floatindex::Float64) = vect[floor(Int, floatindex)]

function regioncenters(regions, numreg, popdens, lonrange, latrange, res)
    println("Finding area- and population-weighted region centers...")
    lons = (-180+res/2:res:180-res/2)         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)          # latitude values (pixel center)
    geocenters = Vector{Tuple{Float64,Float64}}(undef, numreg)
    popcenters = Vector{Tuple{Float64,Float64}}(undef, numreg)
    updateprogress = Progress(numreg, 1)
    for r = 1:numreg
        regmask = (regions.==r)
        latindex = sum(latrange'.*regmask)/sum(regmask)
        lonindex = sum(lonrange.*regmask)/sum(regmask)
        geocenters[r] = lerp(lats, latindex), lerp(lons, lonindex)

        regpopdens = popdens.*regmask
        latindex = sum(latrange'.*regpopdens)/sum(regpopdens)
        lonindex = sum(lonrange.*regpopdens)/sum(regpopdens)
        popcenters[r] = lerp(lats, latindex), lerp(lons, lonindex)
        next!(updateprogress)
    end
    return geocenters, popcenters
end

function getneighbors(regions, r, numreg)
    regmask = (regions.==r)
    regmask_expanded = (imfilter(regmask, diskfilterkernel(2)) .> 0)    # expand mask with 1 pixel
    neighbormask = regmask_expanded .& .!regmask
    neighbors = unique(regions[neighbormask])
    return neighbors[(neighbors.>0) .& (neighbors.!=NOREGION) .& (neighbors.<=numreg) .& (neighbors.!=r)]
end

function connectedregions(regions, offshoreregions, numreg)
    println("\nDetermining onshore and offshore region connections...")
    connected, connectedoffshore = zeros(Bool, numreg, numreg), zeros(Bool, numreg, numreg)
    updateprogress = Progress(numreg, 1)
    for r = 1:numreg     
        neighbors = getneighbors(regions, r, numreg)
        connected[r, neighbors] .= true

        neighborsoffshore = getneighbors(offshoreregions, r, numreg)
        connectedoffshore[r, neighborsoffshore] .= .!connected[r, neighborsoffshore]
        next!(updateprogress)
    end
    return connected, connectedoffshore
end
