using StatsMakie

export shadedmap, hist, GISturbines_density, aggregate

function turbinecharts()
    td, pd = GISturbines_density(gisregion="Denmark5");
    nz = (pd.>1) .| (td.>0);
    GE.plot(GE.histogram(nbins = 50), td[td.>0])
    GE.plot(GE.histogram(nbins = 50), log10.(1 .+ pd[pd.>0]))

    GE.StatsMakie.scatter(td[nz], log10.(pd[nz]), markersize=0.2)

    # shadedmap("Denmark5", log.(1 .+tdm), downsample=2, colorscheme=:magma)
    # shadedmap("Denmark5", log10.(1 .+pd), downsample=2, colorscheme=:magma)
end

function hist(data; normalize=false, args...)
    h = fit(Histogram, data; args...)
    weights = normalize ? h.weights/sum(h.weights) : h.weights
    barplot(h.edges[1][1:end-1], weights)
end

function GISturbines_density(; agg=1, optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    @unpack res, gisregion, scenarioyear = options

    println("\nReading regions and population dataset...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange]
    popdens = aggregate(pop, agg) ./ aggregate(cellarea, agg)'      # persons/km2

    df_DK = DataFrame!(CSV.File(in_datafolder("turbines_DK.csv")))
    df_USA = DataFrame!(CSV.File(in_datafolder("turbines_USA.csv")))
    turbines = vcat(df_DK, df_USA)

    turbinecapacity = aggregate_turbine_capacity(gisregion, turbines, regions, regionlist, lonrange, latrange)
    turbinedensity = aggregate(turbinecapacity, agg) ./ aggregate(cellarea, agg)'    # MW/km2

    # println("\nSaving...")

    # mkpath(in_datafolder("output"))
    # matopen(in_datafolder("output", "GISdata_turbines_$gisregion.mat"), "w") do file
    #     write(file, "turbinecapacity", turbinecapacity)
    # end
    return turbinedensity, popdens
end

# function aggregate(a::AbstractArray, n::Int)
#     sz = size(a) .÷ n
#     out = similar(a, sz)
#     R = CartesianIndices(A)
#     Ifirst, Ilast = first(R), last(R)
#     I1 = oneunit(Ifirst)
#     for I in R
#         n, s = 0, zero(eltype(out))
#         for J in max(Ifirst, I-I1):min(Ilast, I+I1)
#             s += A[J]
#             n += 1
#         end
#         out[I] = s/n
#     end
#     out
# end

function aggregate(a::AbstractArray, n::Int)
    n == 1 && return a
    sz = size(a)
    dims = length(sz)
    out = similar(a, size(a) .÷ n)
    R = CartesianIndices(a)[ntuple(i -> 1:n:sz[i]-n+1, dims)...]   # index to every nth element in all dims
    Rsub = CartesianIndices(ntuple(i -> 0:n-1, dims)) # 
    for (i, I) in enumerate(R)
        s = zero(eltype(out))
        for J in I .+ Rsub
            s += a[J]
        end
        out[i] = s
    end
    out
end


function aggregate_turbine_capacity(gisregion, turbines, regions, regionlist, lonrange, latrange)
    println("Calculate turbine capacity per pixel...")

    numreg = length(regionlist)
    turbinecapacity = zeros(size(regions))       # MW

    for i = 1:size(turbines,1)
        reg, flon, flat = getregion_and_index(turbines[i,:lon], turbines[i,:lat], regions, lonrange, latrange)
        (reg == 0 || reg == NOREGION) && continue
        
        turbinecapacity[flon, flat] += turbines[i,:capac]/1000
    end

    return turbinecapacity
end

function getregion_and_index(lon, lat, regions, lonrange=1:36000, latrange=1:18000)
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    flon = findfirst(round(lon-res2, digits=5) .< lons .<= round(lon+res2, digits=5))
    flat = findfirst(round(lat-res2, digits=5) .< lats .<= round(lat+res2, digits=5))
    (flon == nothing || flat == nothing) && return 0, flon, flat
    # println(lons[flon],", ",lats[flat])
    # println(flon, " ", flat)
    if regions[flon, flat] > 0
        return regions[flon, flat], flon, flat
    else
        return NOREGION, flon, flat
    end
end

function shadedmap(gisregion, plotdata, regionlist, lons, lats, colors, source, dest, xs, ys;
                    resolutionscale=1)
    nreg = length(regionlist)
    scale = maximum(size(plotdata))/6500

    # data[regions.==NOREGION] .= nreg + 1

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    aspect_ratio = (ymax - ymin) / (xmax - xmin)
    pngwidth = round(Int, resolutionscale*1.02*size(plotdata,1))     # allow for margins (+1% on both sides)
    pngsize = pngwidth, round(Int, pngwidth * aspect_ratio)     # use aspect ratio after projection transformation

    println("...constructing map...")
    scene = surface(xs, ys; color=float.(plotdata), colormap=colors,
                            shading=false, show_axis=false, scale_plot=false, interpolate=false)
    ga = geoaxis!(scene, lons[1], lons[end], lats[end], lats[1]; crs=(src=source, dest=dest,))[end]
    ga.x.tick.color = RGBA(colorant"black", 0.4)
    ga.y.tick.color = RGBA(colorant"black", 0.4)
    ga.x.tick.width = scale
    ga.y.tick.width = scale

    println("...saving...")
    mkpath(in_datafolder("output"))
    filename = in_datafolder("output", "$(gisregion)_shadedmap.png")
    isfile(filename) && rm(filename)
    Makie.save(filename, scene, resolution=pngsize)
end

function shadedmap(gisregion, plotdata; resolutionscale=1, downsample=1, colorscheme=:viridis)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    plotdata = plotdata[1:downsample:end, 1:downsample:end]
    lonrange = lonrange[1:downsample:end]
    latrange = latrange[1:downsample:end]
    nreg = length(regionlist)

    colors = colorschemes[colorscheme].colors  #[RGB(0.2,0.3,0.4); colorschemes[:viridis].colors; RGB(0.4,0.4,0.4)]

    println("\nProjecting coordinates (Mollweide)...")
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    source = LonLat()
    dest = Projection("+proj=moll +lon_0=$(mean(lons))")
    xs, ys = xygrid(lons, lats)
    Proj4.transform!(source, dest, vec(xs), vec(ys))

    println("\nOnshore map...")
    shadedmap(gisregion, plotdata, regionlist, lons, lats, colors, source, dest, xs, ys,
        resolutionscale=resolutionscale)
    # exit()
    return nothing
end
