function makedistances(gisregion; scenarioyear="ssp2_2050", res=0.01)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    numreg = length(regionlist)

    println("Finding area- and population-weighted region centers...")
    geocenters, popcenters = getregioncenters(regions, numreg, lonrange, latrange, res, scenarioyear)
    println("\nCalculating distances between centers...")
    distances = [greatcircledistance(Tuple(popcenters[r1,:]), Tuple(popcenters[r2,:])) for r1 = 1:numreg, r2 = 1:numreg]

    println("\nFinding neighboring regions for transmission connections...")
    println("...onshore...")
    connected = connectedregions(regions, numreg)
    println("...offshore...")
    connectedoffshore = connectedregions(offshoreregions, numreg)
    connectedoffshore[connected] .= false
    # regionpop = [sum(pop.==r) for r = 1:numreg]

    println("\nSaving results...")
    matopen(in_datafolder("output", "distances_$gisregion.mat"), "w") do file
        write(file, "distances", distances)
        write(file, "connected", connected)
        write(file, "connectedoffshore", connectedoffshore)
        write(file, "regionlist", string.(regionlist))
        # write(file, "population", pop)
        # write(file, "demand", demand)
        write(file, "regioncenters_lon", popcenters[:,1])
        write(file, "regioncenters_lat", popcenters[:,2])
    end

    # Also save the region file in Matlab format
    # regions2matlab(gisregion)
    nothing
end

# returns great circle distance in km between points given as (lat,lon) tuples (in degrees).
greatcircledistance(point1::Tuple, point2::Tuple) = haversine(point1, point2, 6371.0)

function getregioncenters(regions, numreg, lonrange, latrange, res, scenarioyear)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    mkpath(in_datafolder("output"))

    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange]
    popdens = pop ./ cellarea'

    geocenters, popcenters = regioncenters(regions, numreg, popdens, lonrange, latrange, res)
    return geocenters, popcenters   # column order (lat,lon)
end

function regioncenters(regions, numreg, popdens, lonrange, latrange, res)
    lons = (-180+res/2:res:180-res/2)[lonrange]         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    rows, cols = size(regions)
    geocenters = zeros(numreg,2)
    popcenters = zeros(numreg,2)
    counts = zeros(Int, numreg)
    popdenssum = zeros(numreg)
    for c = 1:cols
        lat = lats[c]
        for r = 1:rows
            reg = regions[r,c]
            (reg == 0 || reg == NOREGION) && continue
            lon = lons[r]
            geocenters[reg,:] += [lat, lon]
            popcenters[reg,:] += popdens[r,c] .* [lat, lon]
            counts[reg] += 1
            popdenssum[reg] += popdens[r,c]
        end
    end
    geocenters ./= counts
    popcenters ./= popdenssum    
    return geocenters, popcenters   # column order (lat,lon)
end

function connectedregions(regions, numreg)
    connected = zeros(Bool, numreg, numreg)
    rows, cols = size(regions)
    for c = 1:cols
        for r = 1:rows
            reg = regions[r,c]
            (reg == 0 || reg == NOREGION) && continue
            for adjc = max(1, c-1):min(cols, c+1), adjr = max(1, r-1):min(rows, r+1)
                reg2 = regions[adjr,adjc]
                (reg2 == 0 || reg2 == NOREGION) && continue
                connected[reg,reg2] = true
            end
        end
    end
    for c = 1:numreg, r = 1:numreg
        connected[r,c] = connected[r,c] || connected[c,r]
    end 
    return connected
end
