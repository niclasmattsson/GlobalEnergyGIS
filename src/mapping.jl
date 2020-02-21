using GeoMakie, Makie, ColorSchemes

export createmaps

function createmap(gisregion, regions, regionlist, lons, lats, colors, source, dest, xs, ys,
                    geocenters, popcenters, connected, connectedoffshore; lines=false, labels=false, resolutionscale=1, textscale=1)
    nreg = length(regionlist)
    scale = maximum(size(regions))/6500

    regions[regions.==NOREGION] .= nreg + 1

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    aspect_ratio = (ymax - ymin) / (xmax - xmin)
    pngwidth = round(Int, resolutionscale*1.02*size(regions,1))                 # allow for margins (+1% on both sides)
    pngsize = pngwidth, round(Int, pngwidth * aspect_ratio)     # use aspect ratio after projection transformation

    println("...constructing map...")
    scene = surface(xs, ys; color=float.(regions), colormap=colors, shading=false, show_axis=false, scale_plot=false)
    ga = geoaxis!(scene, lons[1], lons[end], lats[end], lats[1]; crs=(src=source, dest=dest,))[end]
    ga.x.tick.color = RGBA(colorant"black", 0.4)
    ga.y.tick.color = RGBA(colorant"black", 0.4)
    ga.x.tick.width = scale
    ga.y.tick.width = scale

    if lines
        println("...drawing transmission lines...")
        for (i,conn) in enumerate([connected, connectedoffshore])
            for reg1 = 1:nreg, reg2 = 1:nreg
                reg1 <= reg2 && continue
                !conn[reg1, reg2] && continue
                # line = [lonlatpoint(popcenters[reg1,:]), lonlatpoint(popcenters[reg2,:])]     # straight line on projected map
                line = greatcircletrack(popcenters[reg1,:], popcenters[reg2,:], 50)             # great circle segments
                projectedline = Point2f0.(transform.(source, dest, line))
                color = (i == 1) ? :black : :white
                lines!(scene, projectedline, color=color, linewidth=resolutionscale*scale*3.5)
            end
        end
    end
    if labels
        for reg = 1:nreg
            pos = Point2f0(transform(source, dest, lonlatpoint(geocenters[reg,:])))
            text!(scene, string(regionlist[reg]); position=pos, align=(:center,:center), textsize=textscale*scale*50000)
        end
    end

    println("...saving...")
    datafolder = getconfig("datafolder")
    outputfolder = joinpath(datafolder, "output")
    mkpath(outputfolder)
    filename = joinpath(outputfolder, "$gisregion.png")
    isfile(filename) && rm(filename)
    Makie.save(filename, scene, resolution=pngsize)
    # return scene
end

function createmaps(gisregion; scenarioyear="ssp2_2050", lines=true, labels=true, resolutionscale=1, textscale=1, randseed=1)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    nreg = length(regionlist)

    println("Mapping colors to regions (avoid same color in adjacent regions)...")
    mergeregions = regions + offshoreregions
    mergeconnected = connectedregions(mergeregions, nreg)
    colorindices = greedycolor(mergeconnected, 1:7, 1:nreg, randseed=randseed)
    onshorecolors = [RGB(0.3,0.3,0.45); colorschemes[:Set2_7].colors[colorindices]; RGB(0.5,0.5,0.5)]
    offshorecolors = [RGB(0.5,0.5,0.5); colorschemes[:Set2_7].colors[colorindices]; RGB(0.3,0.3,0.45)]

    println("\nProjecting coordinates (Mollweide)...")
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    source = LonLat()
    dest = Projection("+proj=moll +lon_0=$(mean(lons))")
    xs, ys = xygrid(lons, lats)
    Proj4.transform!(source, dest, vec(xs), vec(ys))

    println("\nFinding interregional transmission lines...")
    geocenters, popcenters = getregioncenters(regions, nreg, lonrange, latrange, res, scenarioyear)   # column order (lat,lon)
    connected = connectedregions(regions, nreg)
    connectedoffshore = connectedregions(offshoreregions, nreg)
    connectedoffshore[connected] .= false

    println("\nOnshore map...")
    scene = createmap(gisregion, regions, regionlist, lons, lats, onshorecolors, source, dest, xs, ys,
        geocenters, popcenters, connected, connectedoffshore, lines=lines, labels=labels, resolutionscale=resolutionscale, textscale=textscale)
    println("\nOffshore map...")
    createmap("$(gisregion)_offshore", offshoreregions, regionlist, lons, lats, offshorecolors, source, dest, xs, ys,
        geocenters, popcenters, connected, connectedoffshore, lines=false, labels=false, resolutionscale=resolutionscale, textscale=textscale)
    return nothing
end

# first color in colorlist that is not in columncolors (which may have many repeated elements)
function firstcolor(columncolors, colorlist)
    colorcounts = zeros(Int, length(colorlist))
    for c in columncolors
        if c > 0
            colorcounts[colorlist[c]] += 1
        end
    end
    goodcolors = colorlist[colorcounts.==0]
    return isempty(goodcolors) ? 0 : rand(goodcolors)
end

function greedycolor(connected, colorlist, order; randseed=0)
    randseed > 0 && Random.seed!(randseed)
    n = size(connected, 1)
    colors = zeros(Int, n)
    for i in order
        columncolors = colors[connected[:,i]]
        col = firstcolor(columncolors, colorlist)
        if col == 0
            println("All colors used, retrying...")
            return greedycolor(connected, colorlist, order; randseed=randseed+1)
        end
        colors[i] = col
    end
    # display(countmap(colors))
    return colors
end

# reverse coordinate order 
lonlatpoint(latlon) = Point2f0(latlon[2], latlon[1])

# returns great circle track between points given as Point2f0 objects (lon, lat order, in degrees).
greatcircletrack(point1::AbstractArray, point2::AbstractArray, segments) = Point2f0.(greatcircletrack(Tuple(point1), Tuple(point2), segments))

# returns great circle track between points given as (lon, lat) tuples (in degrees).
greatcircletrack(point1::Tuple, point2::Tuple, segments) = [reverse(intermediatepoint(point1, point2, f)) for f in LinRange(0, 1, segments)]



# intermediatepoint() ported from Javascript to Julia
# Source: https://www.movable-type.co.uk/scripts/latlong.html, original name 'intermediatePointTo()'
#
# Returns the point on a great circle at given fraction between point1 and point2.
# 
# @param   {Tuple} point1 - (latitude,longitude) of starting point in degrees.
# @param   {Tuple} point2 - (latitude,longitude) of destination point in degrees.
# @param   {number} fraction - Fraction between the two points (0 = starting point, 1 = destination point).
# @returns {Tuple} Intermediate point between this point and destination point (latitude,longitude).
# 
# @example
#   p1 = (52.205, 0.119)
#   p2 = (48.857, 2.351)
#   pInt = intermediatePointTo(p1, p2, 0.25)    # 51.3721°N, 0.7073°E
function intermediatepoint(point1::Tuple, point2::Tuple, fraction)
    φ1, λ1 = point1
    φ2, λ2 = point2
    
    Δφ = φ2 .- φ1   # distance between points
    Δλ = λ2 .- λ1   # distance between points
    a = sind(Δφ/2)*sind(Δφ/2) + cosd(φ1)*cosd(φ2)*sind(Δλ/2)*sind(Δλ/2)
    δ = 2 * atan(sqrt(a), sqrt(1-a))   # radians

    A = sin((1-fraction)*δ) / sin(δ)
    B = sin(fraction*δ) / sin(δ)

    x = A*cosd(φ1)*cosd(λ1) + B*cosd(φ2)*cosd(λ2)
    y = A*cosd(φ1)*sind(λ1) + B*cosd(φ2)*sind(λ2)
    z = A*sind(φ1) + B*sind(φ2)

    lat = atand(z, sqrt(x^2 + y^2))
    lon = atand(y, x)
    return (lat, lon)
end
