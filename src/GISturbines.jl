using CategoricalArrays, XLSX

export shadedmap, myhist, GISturbines_density, aggregate, exportGISturbinedata, landtypes,
        makePixelDataframe, savePixelData, grouppopulationdensity, groupwindspeeds,
        scatter, plot, plot!, Point2f0, RGBA, FRect, plotfix!, heatmap, colorlegend, vbox

const landtypes = ["Evergreen Needleleaf", "Evergreen Broadleaf",
                "Deciduous Needleleaf", "Deciduous Broadleaf", "Mixed Forests",
                "Closed Shrublands", "Open Shrublands", "Woody Savannas", "Savannas", "Grasslands",
                "Permanent Wetlands", "Croplands", "Urban", "Cropland/Natural", "Snow/Ice", "Barren"]

# const landtypes = ["Forest", "Forest", "Forest", "Forest", "Forest", "Shrubland", "Shrubland", "Forest",
#                 "Savanna", "Grassland", "Wetland", "Cropland", "Urban", "Cropland", "Snow/Ice", "Barren"]

function map_protected(; downsample=1, resolutionscale=1, textscale=1, optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, res, protected_codes = options
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected, lonrange, latrange = read_datasets(options)

    n2000_full, ex_n2000 = readraster(in_datafolder("natura2000.tif"), :getextent)
    latrange_n2000, lonrange_n2000 = bbox2ranges(extent2bbox(ex_n2000), 100)
    lons_n2000 = (lonrange[1]:lonrange[end]) .- lonrange_n2000[1] .+ 1
    lats_n2000 = latrange .- latrange_n2000[1] .+ 1
    n2000 = n2000_full[lons_n2000, lats_n2000]

    wdpa = zeros(Bool, size(protected))
    for i in protected_codes
        wdpa[protected .== i] .= true
    end

    reg = (regions.>0) .& (regions .< NOREGION)

    mask = zeros(Int16, size(regions))
    mask[reg] .= 1
    mask[wdpa .& reg] .= 2
    mask[reg .& (n2000 .> 0)] .= 3
    mask[reg .& wdpa .& (n2000 .> 0)] .= 4
    mask[regions.==NOREGION] .= NOREGION

    legendtext = ["unprotected", "WDPA", "Natura2000", "both", "wind turbines"]

    mask = mask[1:downsample:end, 1:downsample:end]
    lonrange = lonrange[1:downsample:end]
    latrange = latrange[1:downsample:end]
    nreg = length(legendtext)

    colors = RGBA.([RGB(0.2,0.3,0.4), RGB(0.7,0.7,0.7), RGB(0,.7,0), RGB(0,0,.7), RGB(.9,.9,0), RGB(1,0,0), RGB(0.4,0.4,0.4)], 0.8)

    println("Mapping colors to regions (avoid same color in adjacent regions)...")
    connected = connectedregions(mask, nreg)
    # colorindices = (nreg > 7) ? greedycolor(connected, 1:7, 1:nreg, randseed=randseed) : collect(1:nreg)
    # colors = colorschemes[Symbol("Set2_$nreg")].colors[colorindices]

    println("\nProjecting coordinates (Mollweide)...")
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    source = Projection("+proj=longlat +datum=WGS84")
    dest = Projection("+proj=moll +lon_0=$(mean(lons)) +ellps=WGS84")
    xs, ys = xygrid(lons, lats)
    Proj4.transform!(source, dest, vec(xs), vec(ys))

    cols = [:lon, :lat, :capac, :year, :onshore, :elec2018]
    df = DataFrame(CSV.File(in_datafolder("turbines_DE.csv")))[:, cols[1:4]]
    tx, ty = df.lon, df.lat
    Proj4.transform!(source, dest, vec(tx), vec(ty))

    println("\nOnshore map...")
    createmap("$(gisregion)_protected", mask, legendtext, lons, lats, colors, source, dest, xs, ys,
        [], [], connected, connected; lines=false, labels=false, resolutionscale=resolutionscale,
        textscale=textscale, legend=true, dots=(tx,ty))
end

function makePixelDataframe(; optionlist...)
    cols = [:lon, :lat, :capac, :year, :onshore, :elec2018]
    df_DK = DataFrame(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
    df_SE = DataFrame(CSV.File(in_datafolder("turbines_SE.csv")))[:, cols[1:5]]
    df_DE = DataFrame(CSV.File(in_datafolder("turbines_DE.csv")))[:, cols[1:4]]
    df_USA = DataFrame(CSV.File(in_datafolder("turbines_USA.csv")))[:, cols[1:4]]
    df_SE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :elec2018] .= missing    # MWh/year
    df_DE[:, :onshore] .= missing    # MWh/year
    df_USA[:, :elec2018] .= missing    # MWh/year
    df_USA[:, :onshore] .= missing    # MWh/year
    df_SE[!,:year] = [ismissing(d) ? 2000 : year(d) for d in df_SE[!,:year]]
    turbines = vcat(df_DK, df_SE, df_DE, df_USA)

    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, res = options
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected, lonrange, latrange = read_datasets(options)

    windatlas = getwindatlas()[lonrange,latrange]
    miuu = zeros(Float32, size(windatlas))

    if gisregion == "SwedenGADM3"
        miuu_raw, ex_miuu = readraster(in_datafolder("miuu_windatlas.tif"), :getextent)
        latrange_miuu, lonrange_miuu = bbox2ranges(extent2bbox(ex_miuu), 100)
        # MIUU has a smaller extent than SwedenGADM3
        lons_miuu = lonrange_miuu .- lonrange[1] .+ 1
        lats_miuu = latrange_miuu .- latrange[1] .+ 1
        miuu[lons_miuu, lats_miuu] = miuu_raw
    end

    n2000, ex_n2000 = readraster(in_datafolder("natura2000.tif"), :getextent)
    latrange_n2000, lonrange_n2000 = bbox2ranges(extent2bbox(ex_n2000), 100)
    lons_n2000 = (lonrange[1]:lonrange[end]) .- lonrange_n2000[1] .+ 1
    lats_n2000 = latrange .- latrange_n2000[1] .+ 1
    n2000 = occursin("USA", gisregion) ? zeros(UInt8, size(windatlas)) : n2000[lons_n2000, lats_n2000]

    nreg = length(regionlist)
    nlon, nlat = size(windatlas)
    turbinecapac = zeros(Int, nlon, nlat)
    nturbines = zeros(Int, nlon, nlat)
    turbineyear = zeros(Int, nlon, nlat)

    for row in eachrow(turbines)
        reg, flon, flat = getregion_and_index(row[:lon], row[:lat], regions, lonrange, latrange)
        if (reg == 0 || reg == NOREGION) 
            reg, flon, flat = getregion_and_index(row[:lon], row[:lat], offshoreregions, lonrange, latrange)
            (reg == 0 || reg == NOREGION) && continue
        end
        nturbines[flon,flat] += 1
        turbineyear[flon,flat] = max(turbineyear[flon,flat], row[:year])
        turbinecapac[flon,flat] += row[:capac]
    end

    region = max.(regions, offshoreregions)
    lons = (-180+res/2:res:180-res/2)[lonrange]   # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]    # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)
    lon = repeat(lons, 1, nlat)
    lat = repeat(lats', nlon)
    area = repeat(cellarea', nlon)

    r = (region .> 0) .& (region .!= NOREGION)

    return DataFrame(lat=lat[r], lon=lon[r], munic=region[r],
        area=round.(area[r], sigdigits=4), landtype=land[r], protected=protected[r], natura2000=Int.(n2000[r]),
        popdens=round.(popdens[r], sigdigits=4), windspeed=round.(windatlas[r], sigdigits=4), miuu=round.(miuu[r], sigdigits=4),
        nturbines=nturbines[r], turbinecapac=turbinecapac[r], turbineyear=turbineyear[r])
end

function savePixelData()
    gisregions = ["SwedenGADM3", "Denmark83", "GermanyGADM3", "USA_Texas", "USA_Iowa", "USA_Oklahoma",
        "USA_Kansas", "USA_Illinois", "USA_South Dakota", "USA_North Dakota", "USA_California",
        "USA_Minnesota", "USA_Colorado", "USA_Oregon"]
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    println("Building municipality list...")
    munic, nmunic = Dict(), Dict()
    for (gisregion, country) in zip(gisregions, countries)
        _, _, regionlist, _, _ = loadregions(gisregion)
        munic[country] = regionlist
        nmunic[gisregion] = length(regionlist)
    end
    dfmunic = vcat((DataFrame(country=country, municipality=munic[country])
                for country in countries)...)
    insertcols!(dfmunic, 1, :index => 1:nrow(dfmunic))
    CSV.write("D:/GISdata/municipalities.csv", dfmunic)
    println("Building main DataFrame...")
    println(countries[1])
    dfall = makePixelDataframe(gisregion=gisregions[1])
    countrynum = 1
    insertcols!(dfall, 3, :country => countrynum)
    len = nmunic[gisregions[1]]
    for gisregion in gisregions[2:end]
        println(gisregion)
        df = makePixelDataframe(gisregion=gisregion)
        df[!,:munic] .+= len
        countrynum += 1
        insertcols!(df, 3, :country => countrynum)
        dfall = vcat(dfall, df)
        len += nmunic[gisregion]
    end
    println("Saving DataFrame...")
    CSV.write("D:/GISdata/windpixeldata.csv", dfall)
end

extent2bbox(ex) = [ex[2] ex[1]; ex[4] ex[3]]

function protected_vs_natura2000(; optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected, lonrange, latrange = read_datasets(options)

    n2000, ex_n2000 = readraster(in_datafolder("natura2000.tif"), :getextent)
    latrange_n2000, lonrange_n2000 = bbox2ranges(extent2bbox(ex_n2000), 100)
    lons_n2000 = (lonrange[1]:lonrange[end]) .- lonrange_n2000[1] .+ 1
    lats_n2000 = latrange .- latrange_n2000[1] .+ 1
    n2000 = n2000[lons_n2000, lats_n2000]

    onshore = (regions .> 0) .& (regions .!= NOREGION) .& (land .> 0)
    offshore = (offshoreregions .> 0) .& (offshoreregions .!= NOREGION) .& (land .== 0)

    count_onshore = zeros(Int, 11, 4)
    count_offshore = zeros(Int, 11, 4)
    for i = 1:length(protected)
        if onshore[i]
            count_onshore[protected[i]+1, n2000[i]+1] += 1
        elseif offshore[i]
            count_offshore[protected[i]+1, n2000[i]+1] += 1
        end
    end
    display(count_onshore)
    display(count_offshore)

    return protected, n2000
end

function miuu_vs_windatlas()
    options = windoptions()
    options[:gisregion] = "SwedenGADM3"
    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land,
        protected, lonrange, latrange = read_datasets(options)

    windatlas = getwindatlas()[lonrange,latrange]
    onshore = (regions .> 0) .& (regions .!= NOREGION) .& (land .> 0)
    offshore = (offshoreregions .> 0) .& (offshoreregions .!= NOREGION) .& (land .== 0)

    miuu, ex_miuu = readraster(in_datafolder("miuu_windatlas.tif"), :getextent)
    latrange_miuu, lonrange_miuu = bbox2ranges(extent2bbox(ex_miuu), 100)
    # MIUU has a smaller extent than SwedenGADM3
    lons_miuu = lonrange_miuu .- lonrange[1] .+ 1
    lats_miuu = latrange_miuu .- latrange[1] .+ 1

    windatlas = windatlas[lons_miuu, lats_miuu]
    onshore = onshore[lons_miuu, lats_miuu] .& (miuu .> 0)
    offshore = offshore[lons_miuu, lats_miuu] .& (miuu .> 0)    

    return miuu, windatlas, onshore, offshore
end

function plotfix!(scene)
    xlabel!(scene, "Global Wind Atlas [m/s]")
    ylabel!(scene, "MIUU [m/s]")
    plot!(scene, [0,14],[0,14], color=:red)
end

#=
mm,ww,on,off = GE.miuu_vs_windatlas();
p = heatmap(reverse(ww.*(mm.>0),dims=2), resolution=(950,950), colorrange=(3,10), highclip=:red, lowclip=:black)
p = heatmap(reverse(mm,dims=2), resolution=(950,950), colorrange=(3,10), highclip=:red, lowclip=:black)
        # cl = colorlegend(p[end], resolution=(950,950)); vbox(p,vl)
opt = (color=RGBA(0,0,0,0.2), markersize=0.01, limits=FRect(0,0,14,14), textsize=10)
p = scatter(Point2f0.(ww[on], mm[on]+randn(sum(on))*.05); opt..., color=RGBA(0,0,0,1)); plotfix!(p)
p = scatter(Point2f0.(ww[off], mm[off]+randn(sum(off))*.05); opt..., color=RGBA(0,0,0,1)); plotfix!(p)

pp,nn = GE.protected_vs_natura2000(gisregion="SwedenGADM3");
p = heatmap(reverse(nn,dims=2), resolution=(950,950), highclip=:red, lowclip=:black)
p = heatmap(reverse(pp,dims=2), resolution=(950,950), highclip=:red, lowclip=:black)

=#

function analyze_protected(; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    df.onshore = df.landtype .> 0
    protected_names = [
       "(Ia) Strict Nature Reserve"
       "(Ib) Wilderness Area"
       "(II) National Park"
       "(III) Natural Monument"
       "(IV) Habitat/Species Management"
       "(V) Protected Landscape/Seascape"
       "(VI) Managed Resource Protected Area"
       "Not Reported"
       "Not Applicable"
       "Not Assigned"
    ]
    protected_type = Dict(i => n for (i,n) in enumerate(protected_names))
    protected_type[0] = "UNPROTECTED"

    gdf = groupby(df, [:countryname, :onshore, :protected])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, [:countryname, :onshore])
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.protected_type = getindex.(Ref(protected_type), cdf.protected)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :protected, order(:onshore, rev=true)])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x)), order(:onshore, rev=true)])
    # CSV.write("protectedinfo.csv", cdf)
    # CSV.write("protectedinfo_tot.csv", cdf_tot)

    filename = "protected WDPA.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end

    return cdf, cdf_tot
end

function analyze_natura2000(; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    df.onshore = df.landtype .> 0
    protected_names = [
       "A: SPAs"
       "B: SCIs and SACs"
       "C: both categories A + B"
    ]
    protected_type = Dict(i => n for (i,n) in enumerate(protected_names))
    protected_type[0] = "UNPROTECTED"

    gdf = groupby(df, [:countryname, :onshore, :natura2000])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, [:countryname, :onshore])
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.protected_type = getindex.(Ref(protected_type), cdf.natura2000)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :natura2000, order(:onshore, rev=true)])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x)), order(:onshore, rev=true)])

    filename = "protected Natura2000.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end
    # CSV.write("protectedinfo.csv", cdf)
    # CSV.write("protectedinfo_tot.csv", cdf_tot)

    return cdf, cdf_tot
end

function analyze_landtype(; firstyear=1978, lastyear=2021, minspeed=0)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame

    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear) .& (df.windspeed .>= minspeed)
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    df.countryname = countries[df.country]
    landtype_names = [
        "Water" 
        "Evergreen Needleleaf Forests"
        "Evergreen Broadleaf Forests"
        "Deciduous Needleleaf Forests"
        "Deciduous Broadleaf Forests"
        "Mixed Forests"
        "Closed Shrublands"
        "Open Shrublands"
        "Woody Savannas"
        "Savannas"
        "Grasslands"
        "Permanent Wetlands"
        "Croplands"
        "Urban"
        "Cropland/Natural"
        "Snow/Ice"
        "Barren"
    ]
    landtype = Dict(i-1 => n for (i,n) in enumerate(landtype_names))

    gdf = groupby(df, [:countryname, :landtype])
    cdf = combine(gdf,
            :area => (a -> sum(a)/1000) => :area,
            [:turbinecapac, :inyears] => ((c,y) -> sum(c.*y)/1000) => :capac
        )
    gdf_tot = groupby(cdf, :countryname)
    cdf_tot = combine(gdf_tot, :area => sum, :capac => sum)
    cdf.turbinedensity = cdf.capac ./ cdf.area
    cdf.landtype_name = getindex.(Ref(landtype), cdf.landtype)
    cdf_tot.turbinedensity = cdf_tot.capac_sum ./ cdf_tot.area_sum

    sort!(cdf, [order(:countryname, by = x -> findfirst(countries.==x)), :landtype])
    sort!(cdf_tot, [order(:countryname, by = x -> findfirst(countries.==x))])

    filename = "landtypes minspeed $minspeed.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [cdf[:,i] for i=1:size(cdf,2)], names(cdf), anchor_cell=XLSX.CellRef("A4"))
        XLSX.writetable!(sheet, [cdf_tot[:,i] for i=1:size(cdf_tot,2)], names(cdf_tot), anchor_cell=XLSX.CellRef("L4"))
    end

    return cdf, cdf_tot
end

function grouppopulationdensity(country; firstyear=1978, lastyear=2021)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame
    df = df[df.country .== country, :]
    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    df.nturbines = df.nturbines .* df.inyears
    df.capac = df.turbinecapac/1000 .* df.inyears
    pops = [[i*10.0^j for j = -1:3 for i in [1,2,5]]; 100_000]
    popmin = [0; pops[1:end-1]]
    df.poprange = CategoricalArrays.cut(df.popdens, pops, extend=true)
    gdf = groupby(df, :poprange)
    cdf = combine(gdf, :nturbines .=> sum, :capac .=> sum, nrow => :pixels, :area .=> sum, 
        renamecols = false)    
    allranges = CategoricalArrays.cut(popmin, pops, extend=true)
    df_out = similar(cdf, length(allranges))
    df_out.poprange = allranges
    df_out[:,2:end] .= 0
    insertcols!(df_out, 2, :popmin => popmin, :popmax => pops)
    indexes = [findfirst(allranges .== r) for r in cdf.poprange]
    df_out[indexes, 4:end] = cdf[:, 2:end]
    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]
    
    filename = "wind_popdens 1978-2021.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$(countries[country]) $firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [df_out[:,i] for i=2:size(df_out,2)], names(df_out)[2:end], anchor_cell=XLSX.CellRef("A1"))
    end
    # CSV.write(filename, df_out[:, 2:end])
end

# just for the histogram data Fredrik wanted
function groupwindspeeds(country; firstyear=1978, lastyear=2021, usemiuu=false)
    df = CSV.File("D:/GISdata/windpixeldata.csv") |> DataFrame
    df = df[df.country .== country, :]
    df.inyears = (df.turbineyear.>=firstyear) .& (df.turbineyear.<=lastyear)
    df.is_onshore = (df.landtype .> 0)
    df.is_offshore = (df.landtype .== 0)
    df.nturbines = df.nturbines .* df.inyears
    df.nturbines_on = df.nturbines .* df.is_onshore
    df.nturbines_off = df.nturbines .* df.is_offshore
    df.capac = df.turbinecapac/1000 .* df.inyears
    df.capac_on = df.capac .* df.is_onshore
    df.capac_off = df.capac .* df.is_offshore
    df.area_on = df.area .* df.is_onshore
    df.area_off = df.area .* df.is_offshore
    df.speed = usemiuu ? df.miuu : df.windspeed
    df.windspeed_range = CategoricalArrays.cut(df.speed, 0:.25:21, extend=true)
    gdf = groupby(df, :windspeed_range)
    cdf = combine(gdf, [:nturbines, :nturbines_on, :nturbines_off] .=> sum, 
        [:capac, :capac_on, :capac_off] .=> sum,
        nrow => :pixels, [:is_onshore, :is_offshore] .=> sum .=> [:pixels_on, :pixels_off],
        [:area, :area_on, :area_off] .=> sum, 
        renamecols = false)
    allranges = CategoricalArrays.cut(0:.25:20.75, 0:.25:21, extend=true)
    df_out = similar(cdf, length(allranges))
    df_out.windspeed_range = allranges
    df_out[:,2:end] .= 0
    indexes = [findfirst(allranges .== r) for r in cdf.windspeed_range]
    df_out[indexes, :] = cdf
    insertcols!(df_out, 2, :windspeed_low => 0:.25:20.75, :windspeed_high => .25:.25:21)

    countries = ["Sweden", "Denmark", "Germany", "Texas", "Iowa", "Oklahoma", "Kansas", "Illinois",
        "South Dakota", "North Dakota", "California", "Minnesota", "Colorado", "Oregon"]

    filename = "winddata 1978-2021.xlsx"
    XLSX.openxlsx(filename, mode=isfile(filename) ? "rw" : "w") do xf
        sheetname = "$(countries[country]) $(usemiuu ? "MIUU " : "")$firstyear-$lastyear"
        sheet = sheetname in XLSX.sheetnames(xf) ? xf[sheetname] : XLSX.addsheet!(xf, sheetname)
        XLSX.writetable!(sheet, [df_out[:,i] for i=2:size(df_out,2)], names(df_out)[2:end], anchor_cell=XLSX.CellRef("A1"))
    end
    # CSV.write(filename, df_out)
end

#=
# run in REPL with using Plots & plotly()
function plothists(region, windspeed, onshoreturbine, onshore, offshore, years)
    region = replace(region, "GADM" => "")
    region = replace(region, r"\d" => "")
    title = "$region $(years[1]) - $(years[2])"
    histopt = (legend=:none, size=(1200,550), tickfont=12, bins=2.8:.1:11.2,
                    xticks=3:1:11, normalize=true)
    barhist(windspeed[onshoreturbine], title=title, c=:lightgray; histopt...)
    stephist!(windspeed[onshore .& (windspeed.>0)], lw=1.5, c=:blue; histopt...)
    stephist!(windspeed[offshore .& (windspeed.>0)], lw=1.5, c=:red; histopt...)
end

    # groupedhist([windspeed[onshoreturbine]; windspeed[offshoreturbine]];
    #     group=[ones(sum(onshoreturbine)); 2*ones(sum(offshoreturbine))],
    #     bar_position = :stack, c=[:lightgray :red], histopt...)
=#

function exportGISturbinedata(; mincapac=0, minclass=0, firstyear=1978, plotmasks=false, optionlist...)
    cols = [:lon, :lat, :capac, :onshore, :elec2018]
    df_DK = DataFrame(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
    df_SE = DataFrame(CSV.File(in_datafolder("turbines_SE.csv")))[:, cols[1:4]]
    df_DE = DataFrame(CSV.File(in_datafolder("turbines_DE.csv")))[:, cols[1:3]]
    df_USA = DataFrame(CSV.File(in_datafolder("turbines_USA.csv")))[:, cols[1:3]]
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
    capac_ok, area, area_ok, maxcapac, maxcapac_ok, class, windspeed =
        zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg), zeros(nreg)

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
            windspeed[reg] += windatlas[flon, flat]
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
    class = class ./ n_onshore
    windspeed = windspeed ./ n_onshore
    offshore = round.(offshore ./ n_offshore * 100, digits=1)
    pop = round.(pop ./ n_onshore ./ capac, digits=1)
    masked = round.(maskedturbines ./ n_onshore * 100, digits=1)
    nogrid = round.(nogrid ./ maskedturbines * 100, digits=1)
    highpop = round.(highpop ./ maskedturbines * 100, digits=1)
    badland = round.(badland ./ maskedturbines * 100, digits=1)
    protected = round.(protected ./ maskedturbines * 100, digits=1)    
    println("Calculating results per municipality/county...")
    updateprogress = Progress(nreg, 1)
    for reg = 1:nreg
        area[reg] += sum((regions .== reg) .* cellarea')
        maxcapac[reg] = area[reg] * onshore_density * area_onshore
        area_ok[reg] += sum(((regions .== reg) .& mask_onshore) .* cellarea')
        maxcapac_ok[reg] = area_ok[reg] * onshore_density * area_onshore
        ss = sortslices([landcount[reg, :] collect(1:16)], dims=1, rev=true)
        commonland[reg,:] = getindex.(Ref(landtypes), ss[1:3, 2])
        freq[reg,:] = round.(ss[1:3, 1]/n_onshore[reg] * 100, digits=1)
        next!(updateprogress)
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
            exploit_tot=round.(area_onshore*capac./maxcapac*100, digits=1),
            exploit_ok=round.(area_onshore*capac_ok./maxcapac_ok*100, digits=1), 
            popdens=pop, land1=commonland[:,1], freq1=freq[:,1], land2=commonland[:,2], freq2=freq[:,2],
            land3=commonland[:,3], freq3=freq[:,3], area=round.(Int, area), class=class,
            windspeed=windspeed)
    outfile = "regionalwindGIS_$(gisregion)_minturbine=$(mincapac)_maxpopdens=$(Int(options.persons_per_km2)).csv"
    CSV.write(in_datafolder("output", outfile), df)
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

function myhist(data; normalize=false, args...)
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

    df_DK = DataFrame(CSV.File(in_datafolder("turbines_DK.csv")))[:, [:lon, :lat, :capac]]
    df_USA = DataFrame(CSV.File(in_datafolder("turbines_USA.csv")))[:, [:lon, :lat, :capac]]
    df_SE = DataFrame(CSV.File(in_datafolder("turbines_SE.csv")))[:, [:lon, :lat, :capac]]
    df_UK = DataFrame(CSV.File(in_datafolder("turbines_UK.csv")))[:, [:lon, :lat, :capac]]
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
    source = Projection("+proj=longlat +datum=WGS84")
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
    turbines = DataFrame(CSV.File(in_datafolder("turbines_DK.csv")))[:, cols]
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
