function vtest()
    fig = Figure()

    ## Voronoi tessellation: Make tessellations from their dual triangulation
    pts = 25randn(2, 500)
    tri = triangulate(pts)
    vorn = voronoi(tri)

    ax = GeoMakie.Axis(fig[2, 2], title="(f): Voronoi tessellation", titlealign=:left, width=400, height=400)
    voronoiplot!(ax, vorn, show_generators=false)

    ## Clipped Voronoi tessellation 
    vorn = voronoi(tri, true)
    ax = GeoMakie.Axis(fig[2, 3], title="(g): Clipped Voronoi tessellation", titlealign=:left, width=400, height=400)
    voronoiplot!(ax, vorn, show_generators=false, color=:white)

    ## Centroidal Voronoi tessellation (CVT)
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1])
    refine!(tri; max_area=1e-3, min_angle=29.871)
    vorn = voronoi(tri)
    smooth_vorn = centroidal_smooth(vorn; maxiters=2500)
    ax = GeoMakie.Axis(fig[2, 4], title="(h): Centroidal Voronoi tessellation", titlealign=:left, width=400, height=400)
    voronoiplot!(ax, smooth_vorn, show_generators=true, markersize=4, colormap=:jet)

    resize_to_layout!(fig)
    fig
end

function vtest2()
    pts = 25 * randn(2, 500)
    tri = triangulate(pts)
    vorn = voronoi(tri)
    clip_vorn = voronoi(tri, true)
    smooth_vorn = centroidal_smooth(clip_vorn)

    cmap = Makie.cgrad(:jet)
    colors = get_polygon_colors(vorn, cmap)
    fig = Figure(fontsize=38, resolution=(1540, 485))
    for (j, vor) in enumerate((vorn, clip_vorn, smooth_vorn))
        ax = Makie.Axis(fig[1, j], width=400, height=400)
        voronoiplot!(ax, vor, strokecolor=:red, strokewidth=0.2, polygon_color=colors, markersize=4)
        xlims!(ax, -100, 100)
        ylims!(ax, -100, 100)
    end
    return fig
end

# Helper method deleted from DelaunayTriangulation after v0.6
function get_polygon_colors(vorn::VoronoiTessellation, cmap)
    F = DelaunayTriangulation.number_type(vorn)
    gtr = [get_generator(vorn, i) for i in each_generator(vorn)]
    gtr_mat = reinterpret(reshape, F, gtr)
    colors = get(cmap, gtr_mat, :extrema)
    return [(a + b) / 2 for (a, b) in eachcol(colors)]
end

function mtest()
    f = Figure()
    ax = GeoMakie.Axis(f[1,1])
    hm = heatmap!(ax, 1:4, 1:4, repeat([2 2 3 1], outer=4), colormap=cgrad(:viridis, 3, categorical=true), shading=NoShading)
    Colorbar(f[:, 2], hm)
    f
end

function mtest2()
    f = Figure()
    ax = GeoAxis(f[1,1])
    hm = surface!(ax, 11:14, 51:54, repeat([2 2 3 1], outer=4), colormap=cgrad(:viridis, 3, categorical=true), shading=NoShading)
    Colorbar(f[:, 2], hm)
    f
end

function geotest_old()
    source = "+proj=longlat +datum=WGS84"
    dest = "+proj=natearth2"
    ptrans = Makie.PointTrans{2}(Proj.Transformation(source, dest, always_xy=true))

    fig = Figure()
    ax = GeoMakie.Axis(fig[1,1], aspect = DataAspect())

    # all input data coordinates are projected using this function
    ax.scene.transformation.transform_func[] = ptrans

    # draw projected grid lines and set limits accordingly
    lats = -90:10.0:90
    lons = -180:10.0:180
    lons = collect(lons)
    lons[end] = prevfloat(lons[end])  # avoid PROJ wrapping 180 to -180
    sz = length(lons), length(lats)
    points = map(CartesianIndices(sz)) do xy
        x, y = Tuple(xy)
        Point2f(lons[x], lats[y])
    end
    limits = Rect2f(Makie.apply_transform(ptrans, points))
    limits!(ax, limits)
    wireframe!(ax, lons, lats, zeros(sz), color=(:gray, 0.2), transparency=true)

    # add black polygons for land area
    url = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
    land = Downloads.download(url * "ne_110m_land.geojson", IOBuffer())
    land_geo = GeoJSON.read(seekstart(land))
    n = length(GeoInterface.features(land_geo))
    poly!(ax, land_geo, color=1:n)

    # add grey dots for populated places
    pop = Downloads.download(url * "ne_10m_populated_places_simple.geojson", IOBuffer())
    pop_geo = GeoJSON.read(seekstart(pop))
    # scatter!(ax, GeoMakie.geo2basic(pop_geo), color="lightgrey", markersize=1.2)

    hidedecorations!(ax)    # x & y axis numbers
    hidespines!(ax)         # outer rectangle

    # return fig
    Makie.save("geomakie_figtest.png", ax.scene, resolution=(1000,600))
    screen = GLMakie.singleton_screen(false)
    close(screen)
    GLMakie.destroy!(screen)
end

function geotest_new()
    # GLMakie.activate!(px_per_unit = 4) # hide

    source = "+proj=longlat +datum=WGS84"
    dest = "+proj=natearth2"

    land = GeoMakie.assetpath("ne_110m_land.geojson")
    land_geo = GeoJSON.read(read(land, String))
    pop = GeoMakie.assetpath("ne_10m_populated_places_simple.geojson")
    pop_geo = GeoJSON.read(read(pop, String))

    begin
        fig = Figure(size = (1000,500))
        ga = GeoAxis(
            fig[1, 1];
            source = source,
            dest = dest
        )

        ga.xticklabelsvisible[] = false
        ga.yticklabelsvisible[] = false
        poly!(ga, land_geo, color=:black)
        popisqrt = map(pop_geo) do geo
            popi = geo.pop_max
            popi > 0 ? sqrt(popi) : 0.0
        end
        mini, maxi = extrema(popisqrt)
        msize = map(popisqrt) do popi
            normed = (popi .- mini) ./ (maxi - mini)
            return (normed * 20) .+ 1
        end
        scatter!(ga, pop_geo; color=popisqrt, markersize=msize)
        fig
    end
end

# Hack for missing method in GeoMakie/MakieCore, https://github.com/MakieOrg/GeoMakie.jl/issues/201
Makie.point_iterator(::Makie.MakieCore.Text{Tuple{Vector{Makie.GlyphCollection}}}) = Point2f[]

function geotest2()
    # GLMakie.activate!(px_per_unit = 4) # hide
    img = rotr90(GeoMakie.earth())
    fig = Figure(size = (1200, 600))
    ga2 = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "Image of Earth")
    surface!(ga2, (-180f0)..(180f0), -90f0..90f0, zeros(size(img)); shading = NoShading, color = img)
    text!(ga2, 30, 45, text = "center", align = (:center, :center))
    save("geotest2.png", fig; px_per_unit=2)    # needs point_iterator() above 
    fig
end

function geotest3() # README example, not quite working
    # GLMakie.activate!(px_per_unit = 4) # hide

    fieldlons = -180:180; fieldlats = -90:90
    field = [exp(cosd(lon)) + 3(lat/90) for lon in fieldlons, lat in fieldlats]

    img = rotr90(GeoMakie.earth())
    land = GeoMakie.land()

    fig = Figure(size = (1000, 1000))

    ga1 = GeoAxis(fig[1, 1]; dest = "+proj=ortho", title = "Orthographic\n ")
    ga2 = GeoAxis(fig[1, 2]; dest = "+proj=moll", title = "Image of Earth\n ")
    ga3 = GeoAxis(fig[2, 1]; title = "Plotting polygons")
    ga4 = GeoAxis(fig[2, 2]; dest = "+proj=natearth", title = "Auto limits") # you can plot geodata on regular axes too

    surface!(ga1, fieldlons, fieldlats, field; colormap = :rainbow_bgyrm_35_85_c69_n256, shading = NoShading)
    lines!(ga1, GeoMakie.coastlines(); transformation = (; translation = (0, 0, 1)))    # maybe transformation only need in GLMakie, not CairoMakie
    image!(ga2, -180..180, -90..90, img; interpolate = false)
    poly!(ga3, land[50:100]; color = 1:51, colormap = (:plasma, 0.5))
    poly!(ga4, land[22]); #datalims!(ga4)

    fig
end

function geotest4()
    # GLMakie.activate!(px_per_unit = 4) # hide

    lons = -180:180
    lats = -90:90
    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]
    
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1],
        dest="+proj=ortho",
        title = "Orthographic projection",
        xticklabelcolor=:red, xgridcolor=:red,
    )
    lp = lines!(ga, GeoMakie.coastlines(); transformation = (; translation = (0, 0, 1)))
    sp = surface!(ga, lons, lats, zeros(size(field)); color=field, shading = NoShading, colormap=:rainbow_bgyrm_35_85_c69_n256)
    cb = Colorbar(fig[1, 2], sp)
    fig
end

function geotest5()
    # GLMakie.activate!(px_per_unit = 4) # hide

    lons = -180:180
    lats = -90:90
    field = [exp(cosd(l)) + 3(y / 90) for l in lons, y in lats]
    
    fig = Figure()
    ax1 = GeoAxis(fig[1, 1], dest = "+proj=vitk1 +lat_1=45 +lat_2=55",title = "vitk1", xgridcolor=:red)
    ax2 = GeoAxis(fig[1, 2]; dest="+proj=wintri", title = "wintri")
    
    surface!(ax1, lons, lats, field; shading = NoShading, colormap = (:plasma, 0.45))
    surface!(ax2, lons, lats, field; shading = NoShading, colormap = (:plasma, 0.45))
    
    fig
end

function geotest6()
    # GLMakie.activate!(px_per_unit = 4) # hide

    lons = -180:180
    lats = -90:90
    # Create some field of values across `lons` and `lats`.
    #
    # This grid can be of any density, but note that the
    # time it takes to plot scales with the grid size!
    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]
    
    fig = Figure()
    ax = GeoAxis(fig[1,1])
    contourf!(ax, lons, lats, field)
    fig
end

function geotest7()
    # GLMakie.activate!(px_per_unit = 4) # hide

    projections = ["+proj=adams_hemi", "+proj=adams_ws1", "+proj=adams_ws2",
    "+proj=aea +lat_1=29.5 +lat_2=42.5", "+proj=aeqd", "+proj=airy", "+proj=aitoff",
    "+proj=apian", "+proj=august", "+proj=bacon", "+proj=bertin1953", "+proj=bipc +ns",
    "+proj=boggs", "+proj=bonne +lat_1=10", "+proj=cass", "+proj=cea",
    "+proj=chamb +lat_1=10 +lon_1=30 +lon_2=40", "+proj=collg", "+proj=comill",
    "+proj=crast", "+proj=denoy", "+proj=eck1", "+proj=eck2", "+proj=eck3",
    "+proj=eck4", "+proj=eck5", "+proj=eck6", "+proj=eqc", "+proj=eqdc +lat_1=55 +lat_2=60",
    "+proj=eqearth", "+proj=euler +lat_1=67 +lat_2=75", "+proj=fahey", "+proj=fouc", "+proj=fouc_s",
    "+proj=gall", "+proj=geos +h=35785831.0 +lon_0=-60 +sweep=y", "+proj=gins8", "+proj=gn_sinu +m=2 +n=3",
    "+proj=goode", "+proj=guyou", "+proj=hammer", "+proj=hatano",
    "+proj=igh", "+proj=igh_o +lon_0=-160", "+proj=imw_p +lat_1=30 +lat_2=-40", "+proj=isea",
    "+proj=kav5", "+proj=kav7", "+proj=laea", "+proj=lagrng", "+proj=larr", "+proj=lask",
    "+proj=lcca +lat_0=35", "+proj=leac", "+proj=loxim",
    "+proj=lsat +ellps=GRS80 +lat_1=-60 +lat_2=60 +lsat=2 +path=2", "+proj=mbt_s", "+proj=mbt_fps",
    "+proj=mbtfpp", "+proj=mbtfpq", "+proj=mbtfps", "+proj=merc", "+proj=mill", "+proj=misrsom +path=1",
    "+proj=moll", "+proj=murd1 +lat_1=30 +lat_2=50",
    "+proj=murd3 +lat_1=30 +lat_2=50", "+proj=natearth", "+proj=natearth2",
    "+proj=nell", "+proj=nell_h", "+proj=nicol",
    "+proj=ob_tran +o_proj=mill +o_lon_p=40 +o_lat_p=50 +lon_0=60", "+proj=ocea", "+proj=oea +m=1 +n=2",
    "+proj=omerc +lat_1=45 +lat_2=55", "+proj=ortel", "+proj=ortho", "+proj=patterson", "+proj=poly",
    "+proj=putp1", "+proj=putp2", "+proj=putp3", "+proj=putp3p", "+proj=putp4p", "+proj=putp5",
    "+proj=putp5p", "+proj=putp6", "+proj=putp6p", "+proj=qua_aut", "+proj=robin", "+proj=rouss",
    "+proj=rpoly", "+proj=sinu", "+proj=times", "+proj=tissot +lat_1=60 +lat_2=65", "+proj=tmerc",
    "+proj=tobmerc", "+proj=tpeqd +lat_1=60 +lat_2=65", "+proj=urm5 +n=0.9 +alpha=2 +q=4",
    "+proj=urmfps +n=0.5", "+proj=vandg", "+proj=vandg2", "+proj=vandg3", "+proj=vandg4",
    "+proj=vitk1 +lat_1=45 +lat_2=55", "+proj=wag1", "+proj=wag2", "+proj=wag3", "+proj=wag4",
    "+proj=wag5", "+proj=wag6", "+proj=wag7", "+proj=webmerc +datum=WGS84", "+proj=weren",
    "+proj=wink1", "+proj=wink2", "+proj=wintri"]
    let k = 1
        fig = Figure(size=(1500, 1500))
        @time for i in 1:10, j in 1:3
            try
                ga = GeoAxis(
                    fig[i, j];
                    dest=projections[k],
                    title="$(projections[k])",
                )
                lines!(ga, GeoMakie.coastlines())
            catch error
                println("Error at iteration $k")
                break
            end
    
            k += 1
        end
        fig
    end
end

function geotest8()
    # https://datahub.io/core/geo-countries#curl # download data from here
    path = GeoMakie.assetpath("vector", "countries.geo.json")
    json_str = read(path, String)
    worldCountries = GeoJSON.read(json_str)
    n = length(worldCountries)
    lons = -180:180
    lats = -90:90
    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]

    fig = Figure(size = (1200,800), fontsize = 22)

    ax = GeoAxis(
        fig[1,1];
        dest="+proj=wintri",
        title = "World Countries",
        tellheight = true,
    )

    hm1 = surface!(ax, lons, lats, field; shading = NoShading)
    translate!(hm1, 0, 0, -10)

    hm2 = poly!(
        ax, worldCountries;     # vector data (probably)!
        color= 1:n,
        colormap = Reverse(:plasma),
        strokecolor = :black,
        strokewidth = 0.25
    )

    cb = Colorbar(fig[1,2]; colorrange = (1, n), colormap = Reverse(:plasma), label = "variable, color code", height = Relative(0.65))

    fig
end
