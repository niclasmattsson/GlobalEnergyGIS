using GeoMakie, Makie

export savemaps, savemaptest

function savemaptest(gisregion)
    lons = LinRange(-179.5, 179.5, 360)
    lats = LinRange(-89.5, 89.5, 180)

    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]

    source = LonLat()
    dest = Projection("+proj=moll +lon_0=0")

    xs, ys = xygrid(lons, lats)
    Proj4.transform!(source, dest, vec(xs), vec(ys))

    scene = surface(xs, ys; color = field, shading = false, show_axis = false, scale_plot = false)
    geoaxis!(scene, -180, 180, -90, 90; crs = (src = source, dest = dest,))

    coastlines!(scene, 1; crs = (src = source, dest = dest,))

    Makie.save("$gisregion.png", scene, resolution=(3840,2160))
    nothing
end

function savemap(gisregion, regions, lons, lats)
    source = LonLat()
    dest = Projection("+proj=moll +lon_0=0")

    xs, ys = xygrid(lons, lats)
    Proj4.transform!(source, dest, vec(xs), vec(ys))

    scene = surface(xs, ys; color=regions, shading=false, show_axis=false, scale_plot=false)
    geoaxis!(scene, lons[1], lons[end], lats[end], lats[1]; crs=(src=source, dest=dest,))

    filename = "$gisregion.png"
    isfile(filename) && rm(filename)
    # display(size(regions))
    Makie.save(filename, scene, resolution=size(regions))
    scene
end

function savemaps(gisregion)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    regions[regions.==NOREGION] .= length(regionlist) + 1
    # regions[offshoreregions.==NOREGION] .= length(regionlist) + 1

    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    savemap(gisregion, float.(regions), lons, lats)
    # savemap("$(gisregion)_offshore", float.(offshoreregions), lons, lats)
end
