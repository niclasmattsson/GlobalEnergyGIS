# using StatsMakie

export shadedmap, hist, GISturbines_density, aggregate, exportGISturbinedata, landtypes

const landtypes = ["Evergreen Needleleaf", "Evergreen Broadleaf",
                "Deciduous Needleleaf", "Deciduous Broadleaf", "Mixed Forests",
                "Closed Shrublands", "Open Shrublands", "Woody Savannas", "Savannas", "Grasslands",
                "Permanent Wetlands", "Croplands", "Urban", "Cropland/Natural", "Snow/Ice", "Barren"]

# const landtypes = ["Forest", "Forest", "Forest", "Forest", "Forest", "Shrubland", "Shrubland", "Forest",
#                 "Savanna", "Grassland", "Wetland", "Cropland", "Urban", "Cropland", "Snow/Ice", "Barren"]

function exportGISturbinedata(; mincapac=0, minclass=0, firstyear=1978, plotmasks=false, optionlist...)
    cols = [:lon, :lat, :capac, :onshore, :elec2018]
    df_DK = DataFrame!(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
    df_SE = DataFrame!(CSV.File(in_datafolder("turbines_SE.csv")))[:, cols[1:4]]
    df_DE = DataFrame!(CSV.File(in_datafolder("turbines_DE.csv")))[:, cols[1:3]]
    df_USA = DataFrame!(CSV.File(in_datafolder("turbines_USA.csv")))[:, cols[1:3]]
    df_SE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :onshore] .= missing    # MWh/year
    df_USA[:, :elec2018] .= missing    # MWh/year
    df_USA[:, :onshore] .= missing    # MWh/year
    turbines = vcat(df_DK, df_SE, df_DE, df_USA)
    turbines = turbines[turbines[!,:capac].>=mincapac, :]

    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, downsample_masks, onshore_density, area_onshore, res = options
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected_mat, lonrange, latrange =
                read_datasets(options)
    mask_onshoreA, mask_onshoreB, mask_offshore,
        gridA, lowpopdens, goodland, protected_area =
            return_wind_masks(options, regions, offshoreregions, gridaccess, popdens,
                            topo, land, protected_mat, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)
    mask_onshore = mask_onshoreA .| mask_onshoreB
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    windatlas = getwindatlas()[lonrange,latrange]
    onshoreclass, offshoreclass = makewindclasses(options, windatlas)

    nreg = length(regionlist)
    capac, offcapac, elec, pop = zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg)
    onshore, offshore = zeros(nreg), zeros(nreg)
    n_onshore, n_offshore = zeros(Int, nreg), zeros(Int, nreg)
    landcount, commonland, freq = zeros(Int, nreg, 16), fill("", nreg, 3), zeros(nreg, 3)
    maskedturbines, nogrid, highpop, badland, protected =
        zeros(Int, nreg), zeros(Int, nreg), zeros(Int, nreg), zeros(Int, nreg), zeros(Int, nreg)
    capac_ok, area, area_ok, maxcapac, maxcapac_ok, class =
        zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg)

    for row in eachrow(turbines)
        reg, flon, flat = getregion_and_index(row[:lon], row[:lat], regions, lonrange, latrange)
        if (reg == 0 || reg == NOREGION) 
            reg, flon, flat = getregion_and_index(row[:lon], row[:lat], offshoreregions, lonrange, latrange)
            (reg == 0 || reg == NOREGION) && continue
        end
        if ismissing(row[:onshore]) || row[:onshore]
            onshoreclass[flon, flat] < minclass && continue
            capac[reg] += row[:capac] / 1000
            capac_ok[reg] += row[:capac] / 1000 * mask_onshore[flon, flat]
            class[reg] += onshoreclass[flon, flat] 
            onshore[reg] += mask_onshore[flon, flat]
            n_onshore[reg] += 1
            if land[flon, flat] != 0
                landcount[reg, land[flon, flat]] += 1
            end
            pop[reg] += capac[reg]*popdens[flon, flat]  # capacity-weighted
            if !mask_onshore[flon, flat]
                maskedturbines[reg] += 1
                nogrid[reg] += !gridA[flon, flat]
                highpop[reg] += !lowpopdens[flon, flat]
                badland[reg] += !goodland[flon, flat]
                protected[reg] += protected_area[flon, flat]
            end
        else
            offshoreclass[flon, flat] < minclass && continue
            offcapac[reg] += row[:capac] / 1000
            offshore[reg] += mask_offshore[flon, flat]
            n_offshore[reg] += 1
        end
        elec[reg] += coalesce(row[:elec2018] / 1000, 0.0)
    end
    freg = (regions.>0) .& (regions.!=NOREGION)
    all_onshore = sum(mask_onshore[freg]) / sum(freg)
    fregoff = (offshoreregions.>0) .& (offshoreregions.!=NOREGION)
    all_offshore = sum(mask_offshore[fregoff]) / sum(fregoff)
    println("\nTotal onshore OK: ", percentstring(sum(onshore)/sum(n_onshore)),
            " ($(percentstring(all_onshore)))")
    println("Total offshore OK: ", percentstring(sum(offshore)/sum(n_offshore)),
            " ($(percentstring(all_offshore)))")
    regionalpopdens = sum(popdens[freg]) / sum(freg)
    println("Masked turbines: ", percentstring(sum(maskedturbines)/sum(n_onshore)))
    println("\tno grid: ", percentstring(sum(nogrid)/sum(maskedturbines)))
    println("\thigh pop: ", percentstring(sum(highpop)/sum(maskedturbines)))
    println("\tbad land: ", percentstring(sum(badland)/sum(maskedturbines)))
    println("\tprotected: ", percentstring(sum(protected)/sum(maskedturbines)))
    println("Total population density (persons/km2): ", round(sum(pop)/sum(n_onshore)/sum(capac), digits=1),
            " ($(round(regionalpopdens, digits=1)))")
    onshore = round.(onshore ./ n_onshore * 100, digits=1)
    class = round.(class ./ n_onshore, digits=1)
    offshore = round.(offshore ./ n_offshore * 100, digits=1)
    pop = round.(pop ./ n_onshore ./ capac, digits=1)
    masked = round.(maskedturbines ./ n_onshore * 100, digits=1)
    nogrid = round.(nogrid ./ maskedturbines * 100, digits=1)
    highpop = round.(highpop ./ maskedturbines * 100, digits=1)
    badland = round.(badland ./ maskedturbines * 100, digits=1)
    protected = round.(protected ./ maskedturbines * 100, digits=1)    
    for reg = 1:nreg
        area[reg] += sum((regions .== reg) .* cellarea')
        maxcapac[reg] = area[reg] * onshore_density * area_onshore
        area_ok[reg] += sum(((regions .== reg) .& mask_onshore) .* cellarea')
        maxcapac_ok[reg] = area_ok[reg] * onshore_density * area_onshore
        ss = sortslices([landcount[reg, :] collect(1:16)], dims=1, rev=true)
        commonland[reg,:] = getindex.(Ref(landtypes), ss[1:3, 2])
        freq[reg,:] = round.(ss[1:3, 1]/n_onshore[reg] * 100, digits=1)
    end
    println("Turbine density of unmasked land (MW/km2): ",
            round(sum(capac_ok)/sum(area_ok), digits=3),
            " ($(round(sum(capac)/sum(area), digits=3)))")
    println("Exploited share of potential capacity (%): ",
            round(sum(capac_ok)/sum(maxcapac_ok)*100, digits=1),
            " ($(round(sum(capac)/sum(maxcapac)*100, digits=1)))")
    ss = sortslices([sum(landcount, dims=1)' collect(1:16)], dims=1, rev=true)
    for i = 1:3
        println("Landtype $i: $(landtypes[ss[i,2]]) ($(percentstring(ss[i,1]/sum(n_onshore))))")
    end
    df = DataFrame(region=regionlist, capac=round.(capac, digits=2),
            offcapac=round.(offcapac, digits=2), elec2018=round.(elec, digits=2),
            onshore=replace(onshore, NaN=>missing), offshore=replace(offshore, NaN=>missing),
            masked=masked, nogrid=nogrid, highpop=highpop, badland=badland, protected=protected,
            dens_tot=round.(capac./area, digits=3), dens_ok=round.(capac_ok./area_ok, digits=3), 
            exploit_tot=round.(capac./maxcapac*100, digits=1), exploit_ok=round.(capac_ok./maxcapac_ok*100, digits=1), 
            popdens=pop, land1=commonland[:,1], freq1=freq[:,1], land2=commonland[:,2], freq2=freq[:,2],
            land3=commonland[:,3], freq3=freq[:,3], area=round.(Int, area), class=class)
    CSV.write(in_datafolder("output", "regionalwindGIS_$gisregion.csv"), df)
    df
end

percentstring(x) = "$(round(100*x, digits=1))%"

function turbinecharts(gisregion)
    td, pd = GISturbines_density(gisregion=gisregion);
    nz = (pd.>1) .| (td.>0);
    scene = plot(histogram(nbins = 50), td[td.>0])
    xlabel!("turbine density (MW/km2)")
    # scene = plot(histogram(nbins = 50), log10.(pd[nz]))
    # xlabel!("pop. density (log10 persons/km2)")

    # scene = StatsMakie.scatter(td[nz], log10.(pd[nz]), markersize=0.2)
    # xlabel!("turbine density (MW/km2)")
    # ylabel!("pop. density (log10 persons/km2)")

    scene

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
    popdens = aggregate_array(pop, agg) ./ aggregate_array(cellarea, agg)'      # persons/km2

    df_DK = DataFrame!(CSV.File(in_datafolder("turbines_DK.csv")))[:, [:lon, :lat, :capac]]
    df_USA = DataFrame!(CSV.File(in_datafolder("turbines_USA.csv")))[:, [:lon, :lat, :capac]]
    df_SE = DataFrame!(CSV.File(in_datafolder("turbines_SE.csv")))[:, [:lon, :lat, :capac]]
    df_UK = DataFrame!(CSV.File(in_datafolder("turbines_UK.csv")))[:, [:lon, :lat, :capac]]
    turbines = vcat(df_DK, df_USA, df_SE, df_UK)

    turbinecapacity = aggregate_turbine_capacity(gisregion, turbines, regions, regionlist, lonrange, latrange)
    turbinedensity = aggregate_array(turbinecapacity, agg) ./ aggregate_array(cellarea, agg)'    # MW/km2

    # println("\nSaving...")

    # mkpath(in_datafolder("output"))
    # matopen(in_datafolder("output", "GISdata_turbines_$gisregion.mat"), "w") do file
    #     write(file, "turbinecapacity", turbinecapacity)
    # end
    return turbinedensity, popdens
end

# function aggregate_array(a::AbstractArray, n::Int)
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

function aggregate_array(a::AbstractArray, n::Int)
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

function return_wind_masks(options, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, persons_per_km2,
                min_shore_distance, max_depth, classB_threshold = options

    println("Creating masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A 
    gridA = (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    mask_onshoreA = gridA .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    mask_onshoreB = (gridB .& .!gridA) .& (popdens .< persons_per_km2) .& goodland .& .!protected_area

    # shoreline mask for offshore wind
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    # all mask conditions
    mask_offshore = gridB .& .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!protected_area

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        masks[(masks .== 0) .& .!gridA .& gridB] .= 8
        masks[(masks .== 0) .& isregion] .= 7
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "no grid", "", "", "wind plant A", "wind plant B"]
        maskmap("$(gisregion)_masks_wind", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)

        # drawmap(.!shore .& (offshoreregions .> 0))
        # drawmap((topo .> -max_depth) .& (offshoreregions .> 0))
        # drawmap(.!protected_area .& (offshoreregions .> 0))
        # drawmap(mask_offshore)
    end
    
    return mask_onshoreA, mask_onshoreB, mask_offshore,
            gridA, (popdens .< persons_per_km2), goodland, protected_area
end


function make_turbineCSV_for_postdocs(; mincapac=0, minclass=0, firstyear=1978, plotmasks=false, optionlist...)
    cols = [:lon, :lat, :capac, :year, :onshore, :elec2018]
    turbines = DataFrame!(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
    turbines = turbines[turbines[:capac].>=mincapac, :]
    nrows = size(turbines, 1)

    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, downsample_masks, onshore_density, area_onshore, res = options
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected_mat, lonrange, latrange =
                read_datasets(options)
    mask_onshoreA, mask_onshoreB, mask_offshore,
        gridA, lowpopdens, goodland, protected_area =
            return_wind_masks(options, regions, offshoreregions, gridaccess, popdens,
                            topo, land, protected_mat, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)
    mask_onshore = mask_onshoreA .| mask_onshoreB
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    windatlas = getwindatlas()[lonrange,latrange]

    nreg = length(regionlist)

    pop, area, area_ok, maxcapac, maxcapac_ok =
        zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg)

    turbines[!, :municipality] = fill("", nrows)
    turbines[!, :avg_windspeed] = fill(0.0, nrows)
    turbines[!, :landtype] = fill("", nrows)
    turbines[!, :popdensity] = fill(0.0, nrows)
    turbines[!, :masks_ok] = fill(false, nrows)

    for row in eachrow(turbines)
        reg, flon, flat = getregion_and_index(row[:lon], row[:lat], regions, lonrange, latrange)
        if (reg == 0 || reg == NOREGION) 
            reg, flon, flat = getregion_and_index(row[:lon], row[:lat], offshoreregions, lonrange, latrange)
            (reg == 0 || reg == NOREGION) && continue
        end
        row[:municipality] = string(regionlist[reg])
        row[:avg_windspeed] = windatlas[flon,flat]
        row[:landtype] = land[flon,flat] > 0 ? landtypes[land[flon,flat]] : "Water"
        row[:popdensity] = popdens[flon, flat]
        row[:masks_ok] = row[:onshore] ? mask_onshore[flon, flat] : mask_offshore[flon, flat]
    end

    for reg = 1:nreg
        area[reg] += sum((regions .== reg) .* cellarea')
        area_ok[reg] += sum(((regions .== reg) .& mask_onshore) .* cellarea')
        pop[reg] = sum((regions .== reg) .* popdens) / area[reg]
    end

    munic = DataFrame(region=regionlist, avg_popdens=pop, area=area, area_masks_ok=area_ok)
    CSV.write("turbines_DK.csv", turbines)
    CSV.write("municipalities_DK.csv", munic)
    turbines, munic
end
