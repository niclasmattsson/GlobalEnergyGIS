using FileIO, GeoMakie, ColorSchemes, Downloads, CairoMakie, Proj, GeoDataFrames, DelaunayTriangulation #, Rasters

GDF = GeoDataFrames
export createmaps, plotmap

using GeoMakie.GeoJSON, GeoMakie.GeoInterface

function createmap(gisregion, regions, regionlist, lons, lats, colors, source, dest,
                    landcenters, popcenters, connected, connectedoffshore;
                    lines=false, labels=false, resolutionscale=1, textscale=1, dotscale=1.5, legend=false, dots=nothing, project=true)
    nreg = length(regionlist)
    scale = maximum(size(regions))/6500

    regions[regions.==NOREGION] .= nreg + 1

    xmin, xmax = extrema(lons)
    ymin, ymax = extrema(lats)
    lims = ((floor(xmin)-0.1, ceil(xmax)+0.1), (floor(ymin)-0.1, ceil(ymax)+0.1))
    aspect_ratio = (ymax - ymin) / (xmax - xmin)
    pngwidth = max(1200, round(Int, resolutionscale*1.02*size(regions,1))) # allow for margins (+1% on both sides)
    pngsize = pngwidth, round(Int, pngwidth * aspect_ratio)     # use aspect ratio after projection transformation

    println("...constructing map...")
    fig = Figure(size=pngsize)
    ga = GeoAxis(fig[1, 1]; dest = "+proj=moll +lon_0=$(mean(lons))", limits=lims)  # or tmerc for Sweden

    if project  # no longer disables projection. heatmap! disables color interpolation but doesn't work in GLMakie, only CairoMakie
        surface!(ga, lons, lats, regions; colormap=cgrad(colors, categorical=true), shading=NoShading)
    else
        heatmap!(ga, lons, lats, regions; colormap=cgrad(colors, categorical=true), shading=NoShading)
    end

    if lines
        println("...drawing transmission lines...")
        for (i,conn) in enumerate([connected, connectedoffshore])
            for reg1 = 1:nreg, reg2 = 1:nreg
                reg1 <= reg2 && continue
                !conn[reg1, reg2] && continue
                # line = [lonlatpoint(popcenters[reg1,:]), lonlatpoint(popcenters[reg2,:])]     # straight line on projected map
                line = greatcircletrack(popcenters[reg1,:], popcenters[reg2,:], 20)             # great circle segments
                color = (i == 1) ? :black : :white
                lines!(Point2f.(line), color=color, linewidth=resolutionscale*scale*textscale*5, overdraw=true)
            end
        end
    end
    if labels
        for reg = 1:nreg
            text!(string(regionlist[reg]); position=lonlatpoint(popcenters[reg,:]), align=(:center,:center),
                    fontsize=textscale*scale*resolutionscale*100, overdraw=true)
        end
    end
    
    if dots != nothing
        tx, ty = dots
        scatter!(tx, ty, markersize=dotscale*resolutionscale, color=RGB(1,0,0))
    end

    println("...saving...")
    mkpath(in_datafolder("output"))
    filename = in_datafolder("output", "$gisregion.png")
    isfile(filename) && rm(filename)
    Makie.save(filename, fig)
    if legend
        makelegend(string.(regionlist), colors[2:end-1], scale=(4*scale)^0.7*resolutionscale)
        img = autocrop(load(filename))
        legend = autocrop(load(in_datafolder("output", "legend.png")))
        # FileIO.save(in_datafolder("output", "legend.png"), legend)
        height = max(size(img,1), size(legend,1)) + 100
        legendwidth = size(legend, 2)
        # display(legend')
        legend = ypad(legend', 362)'
        spacer = fill(RGB{N0f8}(1.0,1.0,1.0), (height, legendwidth÷5))
        combined = [spacer ypad(img, height) spacer ypad(legend, height) spacer]
        FileIO.save(filename, combined)
        rm(in_datafolder("output", "legend.png"))
    end
    nothing
end

function createmaps(gisregion; scenarioyear="ssp2_2050", lines=true, labels=true, resolutionscale=1, textscale=1, randseed=1, downsample=1)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    regions = regions[1:downsample:end, 1:downsample:end]
    offshoreregions = offshoreregions[1:downsample:end, 1:downsample:end]
    lonrange = lonrange[1:downsample:end]
    latrange = latrange[1:downsample:end]
    textscale *= downsample
    nreg = length(regionlist)

    println("Mapping colors to regions (avoid same color in adjacent regions)...")
    mergeregions = copy(regions)
    mergeregions[regions.==0] .= offshoreregions[regions.==0]
    mergeconnected = connectedregions(mergeregions, nreg)
    colorindices = greedycolor(mergeconnected, 1:7, 1:nreg, randseed=randseed)
    onshorecolors = [RGB(0.2,0.3,0.4); colorschemes[:Set2_7].colors[colorindices]; RGB(0.4,0.4,0.4)]
    offshorecolors = [RGB(0.4,0.4,0.4); colorschemes[:Set2_7].colors[colorindices]; RGB(0.2,0.3,0.4)]
    # onshorecolors = RGBA.([RGB(0.2,0.3,0.4); colorschemes[:Set2_7].colors[colorindices]; RGB(0.4,0.4,0.4)], 0.8)
    # offshorecolors = RGBA.([RGB(0.4,0.4,0.4); colorschemes[:Set2_7].colors[colorindices]; RGB(0.2,0.3,0.4)], 0.8)

    println("\nProjecting coordinates (Mollweide)...")
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    source = "+proj=longlat +datum=WGS84"
    dest = "+proj=moll +lon_0=$(mean(lons)) +ellps=WGS84"

    println("\nFinding interregional transmission lines...")
    geocenters, popcenters = getregioncenters(regions, nreg, lonrange, latrange, res, scenarioyear)   # column order (lat,lon)
    landcenters = find_landarea_near_popcenter(regions, nreg, popcenters, lonrange, latrange, res)
    connected = connectedregions(regions, nreg)
    connectedoffshore = connectedregions(offshoreregions, nreg)
    connectedoffshore[connected] .= false

    println("\nOnshore map...")
    createmap(gisregion, regions, regionlist, lons, lats, onshorecolors, source, dest,
        landcenters, popcenters, connected, connectedoffshore; lines, labels, resolutionscale, textscale)
    println("\nOffshore map...")
    createmap("$(gisregion)_offshore", offshoreregions, regionlist, lons, lats, offshorecolors, source, dest,
        landcenters, popcenters, connected, connectedoffshore, lines=false, labels=false, resolutionscale=resolutionscale, textscale=textscale)
    # exit()
    return nothing
end

xygrid(lons, lats) = [lon for lon in lons, lat in lats], [lat for lon in lons, lat in lats]

# ColorBrewer Set2_7:  https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
function maskmap(mapname, regions, regionlist, lonrange, latrange;
                    resolutionscale=1, textscale=1, randseed=1, legend=false, downsample=1)
    regions = regions[1:downsample:end, 1:downsample:end]
    lonrange = lonrange[1:downsample:end]
    latrange = latrange[1:downsample:end]
    nreg = length(regionlist)

    colors = [
        RGB([174,140,114]/255...),  # bad land type
        RGB([230,0,0]/255...),      # high population
        RGB([255,100,255]/255...),  # protected area
        RGB([120,170,80]/255...),   # no grid
        RGB([255,217,47]/255...),   # solar plant A
        RGB([255,150,0]/255...),    # solar plant B
        RGB([140,156,209]/255...),  # wind plant A
        RGB([110,136,169]/255...),  # wind plant B
    ]

    println("Mapping colors to regions (avoid same color in adjacent regions)...")
    connected = connectedregions(regions, nreg)
    # colorindices = (nreg > 7) ? greedycolor(connected, 1:7, 1:nreg, randseed=randseed) : collect(1:nreg)
    # colors = colorschemes[Symbol("Set2_$nreg")].colors[colorindices]
    onshorecolors = [RGB(0.3,0.3,0.45); colors; RGB(0.4,0.4,0.4)]

    println("\nProjecting coordinates (Mollweide)...")
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    source = "+proj=longlat +datum=WGS84"
    dest = "+proj=moll +lon_0=$(mean(lons)) +ellps=WGS84"

    println("\nOnshore map...")
    createmap(mapname, regions, regionlist, lons, lats, onshorecolors, source, dest,
        [], [], connected, connected, lines=false, labels=false, resolutionscale=resolutionscale, textscale=textscale, legend=legend)
    return nothing
end

function autocrop(img, padding::Int=0)
    rowrange = dataindexes_lat(vec(any(img .!= RGB(1.0,1.0,1.0), dims=2)), padding)
    colrange = dataindexes_lat(vec(any(img .!= RGB(1.0,1.0,1.0), dims=1)), padding)   # use _lat both times (_lon does it circularly) 
    return img[rowrange, colrange]
end

function ypad(img, newheight)
    height, width = size(img)
    height == newheight && return img
    height > newheight && error("Can't pad an image to a smaller size.")
    firstrow = (newheight - height) ÷ 2
    newimg = fill(RGB{N0f8}(1.0,1.0,1.0), (newheight, width))
    newimg[firstrow:(firstrow + height - 1), :] = img
    return newimg
end

function makelegend(labels, colors; scale=2)
    i = .!isempty.(labels)
    labels = labels[i]     # don't plot legend entries with empty labels
    colors = colors[i]
    len = length(labels)
    scale = 1
    markerpositions = Point2f.(0, len:-1:1) .* scale
    textpositions = Point2f.(0.06, len:-1:1) .* scale
    dummypositions = Point2f.(1.0, len:-1:1) .* scale
    fig = Figure()
    ax = Makie.Axis(fig[1, 1], aspect = AxisAspect(1), limits = (-0.15, 0.75, 0.5, 10.5)) # allow for 10 lines (to get constant text size for any number of labels)
    hidedecorations!.(ax)
    hidespines!(ax) 
    scatter!(ax,
        markerpositions,
        color = colors,
        marker = :rect,
        markersize = 0.13*scale,
        markerspace = :relative
    )
    scatter!(ax,    # anchor junk to the right to avoid auto axis scaling
        dummypositions,
        color = :black,
        marker = :rect,
        markersize = 0.1*scale,
        markerspace = :relative
    )
    annotations!(ax,
        labels,
        textpositions,
        align = (:left, :center),
        fontsize = 20*scale,
        # space = :data
    )  
    filename = in_datafolder("output", "legend.png")
    isfile(filename) && rm(filename)
    Makie.save(filename, fig)  #, resolution = scene.resolution.val .* scale)
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
            println("All colors used, retrying with randseed=$(randseed+1)...")
            return greedycolor(connected, colorlist, order; randseed=randseed+1)
        end
        colors[i] = col
    end
    # display(countmap(colors))
    return colors
end

# reverse coordinate order 
lonlatpoint(latlon) = Point2f(latlon[2], latlon[1])

# returns great circle track between points given as Point2f objects (lon, lat order, in degrees).
greatcircletrack(point1::AbstractArray, point2::AbstractArray, segments) = greatcircletrack(Tuple(point1), Tuple(point2), segments)

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

plotmap(x; args...) = heatmap(reverse(x, dims=2); scale_plot=false, args...)

function plottraining(countryname)
    df, _ = loadtrainingdata()
    cc = df[:, :country] .== countryname
    normdemand = loaddemanddata()[cc, :normdemand]
    # normalize temp1 to the range [0.5, 1.5] for the 5%-95% quantiles
    normtemp1 = (df[cc,:temp1] .- df[cc,:temp1_qlow]) ./ (df[cc,:temp1_qhigh] .- df[cc,:temp1_qlow]) .+ 0.5
    t = 1:8760
    plot(t, normdemand, color=:black, limits=FRect(0,0,9000,50))
    plot!(t, normtemp1, color=:blue)
end

function plottimezones()
    tzindices, tznames = loadtimezones(1:36000, 1:18000)
    nreg = length(tznames)
    cc = connectedregions(tzindices, nreg)
    colorindices = greedycolor(cc, 1:12, 1:nreg, randseed=494)
    cs = colorschemes[:Set3_12].colors[colorindices]
    plotmap(tzindices[1:5:end, 1:5:end], colormap=cs)
end

function plottimeoffsets()
    tzindices, tznames = loadtimezones(1:36000, 1:18000)
    tz = [n[1:3] == "Etc" ? TimeZone(n, TimeZones.Class(:LEGACY)) : TimeZone(n) for n in tznames]
    offsetlist = tzoffset.(tz)
    offsetvalues = sort(unique(offsetlist))
    offsetdata = [offsetlist[tzi] for tzi in tzindices]
    offsetindices = indexin(offsetdata, offsetvalues)
    nreg = length(offsetvalues)
    cc = connectedregions(offsetindices, nreg)
    colorindices = greedycolor(cc, 1:12, 1:nreg, randseed=1)
    cs = colorschemes[:Set3_12].colors[colorindices]
    plotmap(offsetindices[1:5:end, 1:5:end], colormap=cs)
end

function ehub500()
    buses = CSV.read(in_datafolder("Bus_data_EHUB500 - buses_original.csv"), DataFrame)
    buses = buses[.!(buses.type .== "TRAFO" .&& buses[!, "Voltage [kV]"] .> 220), :]
    unique!(buses, ["x-coordinate", "y-coordinate"])
    xy = Matrix(buses[:, ["x-coordinate", "y-coordinate"]])
    # @show extrema(xy[:,1])
    # @show extrema(xy[:,2])
    bbox = (4.5, 28.9, 55.1, 70.7)

    tri = triangulate(xy')
    vorn = voronoi(tri)

    # fig = Figure()
    # ax = GeoMakie.Axis(fig[1, 1], title="(f): Voronoi tessellation", width=400, height=400)
    # voronoiplot!(ax, vorn, show_generators=false)
    # display(fig)

    n = num_polygons(vorn)
    poly = [Point.(DelaunayTriangulation.get_polygon_coordinates(vorn, i, bbox)) |> GeoMakie.Polygon for i = 1:n] |> GeoMakie.MultiPolygon

    df = DataFrame(geometry=poly, FID=1:n, bus_id=buses.bus_id)
    GDF.write(in_datafolder("ehub500.shp"), df)
    GDF.write(in_datafolder("ehub500.geojson"), df)

    rasterize_ehub500()

    ehub500 = readraster(in_datafolder("ehub500.tif"))
    saveregions("ehub500", buses.bus_id, ehub500)

    return df
end

function ginfo(infile)
    gdalinfo_path() do gdalinfo
        run(`$gdalinfo $infile`)
    end
end

function rasterize_ehub500()
    println("\nRasterizing ehub500 shapefile...")
    shapefile = in_datafolder("ehub500.shp")
    outfile = in_datafolder("ehub500.tif")
    options = "-a FID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize(shapefile, outfile, split(options, ' '))
 
    # println("Creating .csv file for regional index and name lookup...")
    # sql = "select name from ehub500"
    # # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    # outfile = in_datafolder("ehub500.csv")
    # ogr2ogr_path() do ogr2ogr
    #     @time run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
    # end

    nothing
end

function inpolys(point, polygons)
    for (i, poly) in enumerate(polygons)
        if ArchGDAL.contains(poly, point)
            return i
        end
    end
    return 0
end

function readhydro()
    infile = in_datafolder("hydro-power-database", "data", "jrc-hydro-power-plant-database.csv")
    hydro = CSV.read(infile, DataFrame)
    df = GDF.read(in_datafolder("ehub500.geojson"))

    hydro.bus_id = copy(hydro.GEO)
    select!(hydro, 1:2, [:bus_id, :installed_capacity_MW, :storage_capacity_MWh, :avg_annual_generation_GWh], :)

    for hplant in eachrow(hydro)
        if hplant.lat < 54
            hplant.bus_id = missing
            continue
        end
        coords = ArchGDAL.createpoint(hplant.lon, hplant.lat)
        ndx = inpolys(coords, df.geometry)
        hplant.bus_id = (ndx > 0) ? df.bus_id[ndx] : missing
    end

    hydro = hydro[.!ismissing.(hydro.bus_id), :]
    hydro.installed_capacity_MW .= replace(hydro.installed_capacity_MW, missing => 0.0)
    hydro.storage_capacity_MWh .= replace(hydro.storage_capacity_MWh, missing => 0.0)
    hydro.avg_annual_generation_GWh .= replace(hydro.avg_annual_generation_GWh, missing => 0.0)

    allbuses = CSV.read(in_datafolder("Bus_data_EHUB500 - buses_original.csv"), DataFrame)[!, :bus_id] |> sort
    extrabuses = setdiff(allbuses, hydro.bus_id)

    return extrabuses

    hydro_gdf = groupby(hydro, :bus_id)
    hydro_sums = DataFrames.combine(hydro_gdf, :installed_capacity_MW => sum, :storage_capacity_MWh => sum, :avg_annual_generation_GWh => sum)
    sort!(hydro_sums, :bus_id)

    CSV.write(in_datafolder("hydro_nordic.csv"), hydro)
    CSV.write(in_datafolder("hydro_nordic_sums.csv"), hydro_sums)

    return hydro, hydro_sums
end

function matlab2ehub()
    df_ehub = GDF.read(in_datafolder("ehub500.geojson"))
    soldata = matread(in_datafolder("output", "GISdata_solar2019_ehub500.mat"))
    winddata = matread(in_datafolder("output", "GISdata_wind2019_ehub500.mat"))

    allbuses = CSV.read(in_datafolder("Bus_data_EHUB500 - buses_original.csv"), DataFrame)[!, :bus_id] |> sort
    extrabuses = setdiff(allbuses, df_ehub.bus_id)
    df_extra = DataFrame(zeros(8760, length(extrabuses)), string.(extrabuses))

    hours = 'h' .* lpad.(string.(1:8760), 4, '0')
    df_hours = DataFrame("" => hours)
    busid = string.(df_ehub.bus_id)

    df_pvp = DataFrame(max.(0, soldata["CFtime_pvplantA"])[:,:,1], busid)
    df_pvr = DataFrame(max.(0, soldata["CFtime_pvrooftop"])[:,:,1], busid)
    df_won = DataFrame(max.(0, winddata["CFtime_windonshoreA"])[:,:,1], busid)
    df_wof = DataFrame(max.(0, winddata["CFtime_windoffshore"])[:,:,1], busid)

    caps = [soldata["capacity_pvplantA"] soldata["capacity_pvrooftop"] winddata["capacity_onshoreA"] winddata["capacity_offshore"]]
    dfcap = DataFrame([busid caps], ["bus_id", "PV", "PVR", "WON", "WOFF"])
    extra2 = DataFrame([string.(extrabuses) zeros(length(extrabuses), 4)], ["bus_id", "PV", "PVR", "WON", "WOFF"])
    df_caps = sort(vcat(dfcap, extra2), :bus_id, by = x->parse(Int, x))

    df_pvp = select(hcat(df_hours, df_pvp, df_extra), [""; string.(allbuses)])
    df_pvr = select(hcat(df_hours, df_pvr, df_extra), [""; string.(allbuses)])
    df_won = select(hcat(df_hours, df_won, df_extra), [""; string.(allbuses)])
    df_wof = select(hcat(df_hours, df_wof, df_extra), [""; string.(allbuses)])

    df_pvp[!, 2:end] .= replace(Matrix(df_pvp[!, 2:end]), NaN => 0.0)
    df_pvr[!, 2:end] .= replace(Matrix(df_pvr[!, 2:end]), NaN => 0.0)
    df_won[!, 2:end] .= replace(Matrix(df_won[!, 2:end]), NaN => 0.0)
    df_wof[!, 2:end] .= replace(Matrix(df_wof[!, 2:end]), NaN => 0.0)

    CSV.write(in_datafolder("ehub_profile_pv.csv"), df_pvp, delim=';', decimal=',')
    CSV.write(in_datafolder("ehub_profile_pvr.csv"), df_pvr, delim=';', decimal=',')
    CSV.write(in_datafolder("ehub_profile_won.csv"), df_won, delim=';', decimal=',')
    CSV.write(in_datafolder("ehub_profile_wof.csv"), df_wof, delim=';', decimal=',')
    CSV.write(in_datafolder("ehub_capacities.csv"), df_caps, delim=';', decimal=',')

    nothing
end
