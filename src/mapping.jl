using GeoMakie, Makie

export savemaps

function savemaps(gisregion)
	# regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)

	lons = LinRange(-179.5, 179.5, 360)
	lats = LinRange(-89.5, 89.5, 180)

	field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]

	source = LonLat()
	dest = WinkelTripel()

	xs, ys = xygrid(lons, lats)
	Proj4.transform!(source, dest, vec(xs), vec(ys))

	scene = surface(xs, ys; color = field, shading = false, show_axis = false, scale_plot = false)

	geoaxis!(scene, -180, 180, -90, 90; crs = (src = source, dest = dest,))

	coastlines!(scene, 1; crs = (src = source, dest = dest,))
end
