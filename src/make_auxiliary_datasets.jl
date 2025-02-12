export rasterize_datasets, create_scenario_datasets, cleanup_datasets, makeprotected, savelandcover,
        createGDP, creategridaccess, getpopulation, getwindatlas, similarity, closeststring

# cleanup options: :none, :limited, :all
function rasterize_datasets(; cleanup=:all)
    rasterize_GADM()
    rasterize_NUTS()
    rasterize_protected()
    downscale_landcover()
    savelandcover()
    upscale_topography()
    saveregions_global()
    rasterize_timezones()
    maketimezones()
    cleanup_datasets(cleanup=cleanup)
end

function create_scenario_datasets(scen, year)
    if !isfile(in_datafolder("population_$(scen)_$year.jld"))
        println("\nCreating population dataset for $scen $year...")
        downscale_population(scen, year)
    end
    if !isfile(in_datafolder("gdp_$(scen)_$year.jld"))
        println("\nCreating GDP dataset for $scen $year...")
        createGDP(scen, year)
    end
    if !isfile(in_datafolder("gridaccess_$(scen)_$year.jld"))
        println("\nCreating grid access dataset for $scen $year...")
        creategridaccess(scen, year)
    end
end

# cleanup options: :none, :limited, :all
function cleanup_datasets(; cleanup=:all)
    cleanup == :none && return
    for i = 0:2
        rm(in_datafolder("WDPA", "protected_raster$i.tif"), force=true)
        rm(in_datafolder("WDPA", "protectedfields$i.csv"), force=true)
        rm(in_datafolder("WDPA", "protected.jld$i"), force=true)
        rm(in_datafolder("WDPA", "WDPA-shapefile$i"), force=true, recursive=true)
    end
    rm(in_datafolder("landcover.tif"), force=true)
    rm(in_datafolder("topography.tif"), force=true)
    rm(in_datafolder("timezones.tif"), force=true)
    rm(in_datafolder("timezone_names.csv"), force=true)
    if cleanup == :all
        rm(in_datafolder("Landcover - USGS MODIS.tif"), force=true)
        rm(in_datafolder("ETOPO1_Ice_c_geotiff.tif"), force=true)
        rm(in_datafolder("gadm36"), force=true, recursive=true)
        rm(in_datafolder("nuts2016-level3"), force=true, recursive=true)
        rm(in_datafolder("WDPA"), force=true, recursive=true)
        rm(in_datafolder("timezones-with-oceans.shapefile"), force=true, recursive=true)
    end
end

function rasterize_GADM()
    println("\nRasterizing GADM shapefile for global administrative areas (1-10 minute run time)...")
    shapefile = in_datafolder("gadm36", "gadm36.shp")
    outfile = in_datafolder("gadm.tif")
    options = "-a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize(shapefile, outfile, split(options, ' '))
 
    println("Creating .csv file for regional index and name lookup...")
    sql = "select uid,name_0,name_1,name_2 from gadm36"
    # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    outfile = in_datafolder("gadmfields.csv")
    ogr2ogr_path() do ogr2ogr
        @time run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
    end
    nothing
end

function rasterize_NUTS()
    println("\nRasterizing NUTS shapefile for European administrative areas...")
    name = "NUTS_RG_01M_2016_4326_LEVL_3"
    shapefile = in_datafolder("nuts2016-level3", "$name.shp")
    outfile = in_datafolder("nuts.tif")
    options = "-a ROWID -ot Int16 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW -dialect SQlite"
    sql = "select ROWID+1 AS ROWID,* from $name"
    @time rasterize(shapefile, outfile, split(options, ' '), sql=sql)
 
    println("Creating .csv file for regional index and name lookup...")
    outfile = in_datafolder("nutsfields.csv")
    sql = "select ROWID+1 AS ROWID,* from $name"
    ogr2ogr_path() do ogr2ogr
        @time run(`$ogr2ogr -f CSV $outfile -dialect SQlite -sql $sql $shapefile`)
    end
    nothing
end

function read_gadm()
    println("Reading GADM rasters...")
    gadmfields = readdlm(in_datafolder("gadmfields.csv"), ',', header=true)[1]
    imax = maximum(gadmfields[:,1])
    subregionnames = fill("", (imax,3))
    subregionnames[gadmfields[:,1],:] = string.(gadmfields[:,2:4])
    gadm = readraster(in_datafolder("gadm.tif"))
    return gadm, subregionnames
end

function read_nuts()
    println("Reading NUTS rasters...")
    nutsfields = readdlm(in_datafolder("nutsfields.csv"), ',', header=true)[1]
    imax = maximum(nutsfields[:,1])
    subregionnames = nutsfields[:,3]    # indexes of NUTS regions are in order 1:2016, let's use that 
    nuts = readraster(in_datafolder("nuts.tif"))
    return nuts, subregionnames
end

function rasterize_protected()
    println("\nRasterizing three WDPA shapefiles for protected areas (total run time 6 minutes - 2 hours)...")

    for i = 0:2
        println("\nFile $(i+1)/3:")
        shapefile = in_datafolder("WDPA", "WDPA-shapefile_$i", "WDPA-shapefile-polygons.shp")

        println("Rasterizing...")
        gdal_rasterize_path() do gdal_rasterize
            outfile = in_datafolder("WDPA", "protected_raster$i.tif")
            sql = "select FID from \"WDPA-shapefile-polygons\""
            options = "-a FID -a_nodata -1 -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
            @time run(`$gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)
        end

        println("Creating .csv file for WDPA index and name lookup...")
        ogr2ogr_path() do ogr2ogr
            outfile = in_datafolder("WDPA", "protectedfields$i.csv")
            sql = "select FID,IUCN_CAT from \"WDPA-shapefile-polygons\""
            run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
        end

        makeprotected(i)
    end

    println("\nMerging the three rasters...")
    protected = max.(
        JLD.load(in_datafolder("WDPA", "protected0.jld"), "protected"),
        JLD.load(in_datafolder("WDPA", "protected1.jld"), "protected"),
        JLD.load(in_datafolder("WDPA", "protected2.jld"), "protected")
    )
    JLD.save(in_datafolder("protected.jld"), "protected", protected, compress=true)
    println("Done.")

    nothing
end

function makeprotected(n)
    println("Reading rasters...")
    protectedfields = readdlm(in_datafolder("WDPA", "protectedfields$n.csv"), ',', header=true)[1]
    IUCNcodes = ["Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Reported", "Not Applicable", "Not Assigned"]
    IUCNlookup = Dict(c => i for (i,c) in enumerate(IUCNcodes))
    protected0 = readraster(in_datafolder("WDPA", "protected_raster$n.tif"))

    println("Converting indexes to protected area types...")
    protected = similar(protected0, UInt8)
    # could replace loop with:  map!(p -> p == -1 ? 0 : IUCNlookup[protectedfields[p+1,2], protected, protected0)
    # alternatively             map!(p -> ifelse(p == -1, 0, IUCNlookup[protectedfields[p+1,2]), protected, protected0)
    @time for (i, p) in enumerate(protected0)
        protected[i] = (p == -1) ? 0 : IUCNlookup[protectedfields[p+1,2]]
    end
    println("Saving...")
    # JLD.save(in_datafolder("protected$n.jld"), "protected", protected, compress=true)
    JLD.save(in_datafolder("WDPA", "protected$n.jld"), "protected", protected)
end

function rasterize_timezones()
    println("Rasterizing shapefile of time zones...")
    shapefile = in_datafolder("timezones-with-oceans-now.shapefile", "combined-shapefile-with-oceans-now.shp")
    sql = "select FID+1 as FID from \"combined-shapefile-with-oceans-now\""
    outfile = in_datafolder("timezones.tif")
    options = "-a FID -a_nodata 0 -ot Int16 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    gdal_rasterize_path() do gdal_rasterize
        @time run(`$gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)
    end

    println("Creating .csv file for time zone index and name lookup...")
    sql = "select FID+1 as FID,tzid from \"combined-shapefile-with-oceans-now\""
    outfile = in_datafolder("timezone_names.csv")
    ogr2ogr_path() do ogr2ogr
        @time run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
    end
    nothing
end

function maketimezones()
    println("Reading time zone raster file...")
    tznames = string.(readdlm(in_datafolder("timezone_names.csv"), ',', header=true)[1][:,2])
    timezones = readraster(in_datafolder("timezones.tif"))

    f0 = findall(timezones.==0)     # find any "no data" pixels (should only be one pixel in Canada near Halifax)
    for f in f0
        timezones[f] = maximum(timezones[f .+ CartesianIndices((-1:1, -1:1))])  # replace with largest neighbor
    end

    println("Saving time zones dataset...")
    JLD.save(in_datafolder("timezones.jld"), "timezones", timezones, "tznames", tznames, compress=true)
end

function loadtimezones(lonrange, latrange)
    jldopen(in_datafolder("timezones.jld"), "r") do file
        return read(file, "timezones")[lonrange, latrange], read(file, "tznames")
    end
end

function resample(infile::String, outfile::String, options::Vector{<:AbstractString})
    gdal_translate_path() do gdal_translate
        @time run(`$gdal_translate $options -co COMPRESS=LZW $infile $outfile`)
    end
end

function downscale_landcover()
    println("\nDownscaling landcover dataset (2-10 minutes)...")
    infile = in_datafolder("Landcover - USGS MODIS.tif")
    options = "-r mode -ot Byte -tr 0.01 0.01"
    resample(infile, in_datafolder("landcover.tif"), split(options, ' '))
    nothing
end

function savelandcover()
    println("Converting landcover dataset from TIFF to JLD...")
    landcover = readraster(in_datafolder("landcover.tif"))
    landtypes = [
        "Water", "Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests", "Deciduous Needleleaf Forests", "Deciduous Broadleaf Forests", 
        "Mixed Forests", "Closed Shrublands", "Open Shrublands", "Woody Savannas", "Savannas", "Grasslands", "Permanent Wetlands",
        "Croplands", "Urban", "Cropland/Natural", "Snow/Ice", "Barren"
    ]
    landcolors = 1/255 * [
        190 247 255; 0 100 0; 77 167 86; 123 204 6; 104 229 104;
        55 200 133; 216 118 118; 255 236 163; 182 231 140; 255 228 18; 255 192 107; 40 136 213; 
        255 255 0; 255 0 0; 144 144 0; 255 218 209; 190 190 190; 
    ]
    println("Saving landcover dataset...")
    JLD.save(in_datafolder("landcover.jld"), "landcover", landcover, "landtypes", landtypes,
                "landcolors", landcolors, compress=true)
end

function upscale_topography()
    println("\nUpscaling topography dataset...")
    infile = in_datafolder("ETOPO1_Ice_c_geotiff.tif")
    options = "-r cubicspline -tr 0.01 0.01"
    outfile = in_datafolder("topography.tif")
    resample(infile, outfile, split(options, ' '))
    println("Reading new topography raster...")
    topography = readraster(outfile)
    println("Saving topography dataset...")
    JLD.save(in_datafolder("topography.jld"), "topography", topography, compress=true)
end

# gettopography() = readraster("topography.tif")

function downscale_population(scen, year)
    scen = lowercase(scen)
    println("Reading population dataset...")
    dataset = Dataset(in_datafolder("SSP_1km", "$(scen)_total_$year.nc4"))
    pop = replace(dataset["Band1"][:,:], missing => Float32(0))

    lat = dataset["lat"][:]
    res = 0.5/60    # source resolution 0.5 arcminutes

    println("Padding and saving intermediate dataset...")
    skiptop = round(Int, (90-(lat[end]+res/2)) / res)
    skipbottom = round(Int, (lat[1]-res/2-(-90)) / res)
    nlons = size(pop,1)
    # the factor (.01/res)^2 is needed to conserve total population
    pop = [zeros(Float32,nlons,skiptop) reverse(pop, dims=2)*Float32((.01/res)^2) zeros(Float32,nlons,skipbottom)]
    temptiff = "$(tempname()).tif"
    temptiff2 = "$(tempname()).tif"
    saveTIFF(pop, temptiff)

    println("Downscaling population dataset...")
    options = "-r cubicspline -tr 0.01 0.01"
    resample(temptiff, temptiff2, split(options, ' '))
    newpop = readraster(temptiff2)

    println("Saving population dataset...")
    JLD.save(in_datafolder("population_$(scen)_$year.jld"), "population", newpop, compress=true)

    rm(temptiff, force=true)
    rm(temptiff2, force=true)
end

getpopulation(scen, year) = JLD.load(in_datafolder("population_$(scen)_$year.jld"), "population")

function createGDP(scen, year)
    scen = lowercase(scen)
    scennum = scen[end]
    println("Reading low resolution population and GDP datasets...")
    pop, extent = readraster(in_datafolder("global_population_and_gdp", "p$(scennum)_$year.tif"), :getextent) # million people
    gdp = readraster(in_datafolder("global_population_and_gdp", "g$(scennum)_$year.tif"))    # billion USD(2005), PPP

    # Convert to USD 2010 using US consumer price index (CPI-U). CPI-U 2005: 195.3, CPI-U 2010: 218.056 
    # https://www.usinflationcalculator.com/inflation/consumer-price-index-and-annual-percent-changes-from-1913-to-2008/
    tempfile = tempname()
    gdp_per_capita = gdp./pop * 218.056/195.3 * 1000    # new unit: USD(2010)/person, PPP
    gdp_per_capita[pop.<=0] .= 0    # non land cells have pop & gdp set to -infinity, set to zero instead
    saveTIFF(gdp_per_capita, tempfile, extent)

    # println("Downscaling to high resolution and saving...")
    # options = "-r average -tr 0.01 0.01"
    # resample(tempfile, "gdp_per_capita_$(scen)_$year.tif", split(options, ' '))
    # rm(tempfile)

    @time gdphigh = downscale_lowres_gdp_per_capita(tempfile, scen, year)     # unit: USD(2010)/grid cell, PPP
    rm(tempfile, force=true)
    println("Saving high resolution GDP...")
    JLD.save(in_datafolder("gdp_$(scen)_$year.jld"), "gdp", gdphigh, compress=true)
end

function downscale_lowres_gdp_per_capita(tempfile, scen, year)
    println("Create high resolution GDP set using high resolution population and low resolution GDP per capita...")
    gpclow, extent = readraster(tempfile, :extend_to_full_globe)
    pop = getpopulation(scen, year)
    gdphigh = similar(pop, Float32)

    nrows, ncols = size(gpclow)
    sizemult = size(pop,1) ÷ nrows

    for c = 1:ncols
        cols = (c-1)*sizemult .+ (1:sizemult)
        for r = 1:nrows
            gpc = gpclow[r,c]
            rows = (r-1)*sizemult .+ (1:sizemult)
            gdphigh[rows,cols] = gpc * pop[rows,cols]
        end
    end
    return gdphigh
end

function creategridaccess(scen, year)
    println("Estimate high resolution grid access dataset by filtering gridded GDP...")
    gdp = JLD.load(in_datafolder("gdp_$(scen)_$year.jld"), "gdp")
    res = 360/size(gdp,1)

    disk = diskfilterkernel(1/6/res)                        # filter radius = 1/6 degrees
    gridaccess = gridsplit(gdp .> 100_000, x -> imfilter(x, disk), Float32)
    # gridaccess = Float32.(imfilter(gdp .> 100_000, disk))   # only "high" income cells included (100 kUSD/cell), cell size = 1x1 km          
    println("\nCompressing...")
    selfmap!(x -> ifelse(x<1e-6, 0, x), gridaccess)         # force small values to zero to reduce dataset size
    println("Saving high resolution grid access...")
    JLD.save(in_datafolder("gridaccess_$(scen)_$year.jld"), "gridaccess", gridaccess, compress=true)

    # maybe better:
    # loop through countries, index all pixels into vector, sort by GDP, use electrification to assign grid access
end

function getwindatlas(altitude=100)
    # filename = in_datafolder("gwa3_250_wind-speed_100m.tif") # v3.0 (lon extent [-180.3, 180.3], strangely)
    # filename = in_datafolder("global_ws.tif") # v2.3 (lon extent [-180.3, 180.3], strangely)
    # filename = in_datafolder("Global Wind Atlas v1 - 100m wind speed.tif")   # v1.0
    # windatlas = readraster(filename, :extend_to_full_globe)[1]
    filename = in_datafolder("Global Wind Atlas v3 - $(altitude)m wind speed.tif")   # v3.0
    windatlas = readraster(filename)
    clamp!(windatlas, 0, 25)
end

# Convert the Global Wind Atlas 3.0 dataset from 250 m to 1 km resolution. This reduces file size
# from 13 GB to 1 GB. Also change its weird lon-lat extents to standard [-180,-90] - [180, 90].
# Interpolate using bilinear, which is a simple and robust choice for GIS applications that avoids
# creating artifacts. See: 
# https://gis.stackexchange.com/questions/10931/what-is-lanczos-resampling-useful-for-in-a-spatial-context
function downsample_windatlas3(altitude=100)
    infile = in_datafolder("gwa3_250_wind-speed_$(altitude)m.tif")
    gdalinfo_path() do gdalinfo
        run(`$gdalinfo $infile`)
    end
    println("\n")
    outfile = in_datafolder("Global Wind Atlas v3 - $(altitude)m wind speed.tif")
    options = split("-r bilinear -te -180 -90 180 90 -tr 0.01 0.01", ' ')
    gdalwarp_path() do gdalwarp
        @time run(`$gdalwarp $options -co COMPRESS=LZW $infile $outfile`)
    end
end

function ogrinfo(file, options=["-al", "-so"])
    ogrinfo_path() do ogrinfo
        @time run(`$ogrinfo $options $file`)
    end
end

# GE.ogrinfo("D:/GISdata/Natura2000_end2019_Shapefile/Natura2000_end2019_epsg3035.shp")
# GE.ogrinfo("D:/GISdata/Natura2000_end2019_Shapefile/Natura2000_end2019_epsg3035.shp", ["-dialect", "sqlite",  "-sql", "select sitecode, sitename, release_da, ms, sitetype, inspire_id from Natura2000_end2019_epsg3035 limit 10"])

function rasterize_Natura2000()
    shapefile = in_datafolder("Natura2000_end2019_Shapefile", "Natura2000_end2019_epsg3035.shp")
    shapefile_proj = in_datafolder("Natura2000_end2019_Shapefile", "Natura2000_reprojected.shp")

    # Couldn't figure out how to make gdal_rasterize reproject on the fly, so...
    if !isfile(shapefile_proj)
        println("Reprojecting shapefile to EPSG:4326...")
        ogr2ogr_path() do ogr2ogr
            @time run(`$ogr2ogr -t_srs epsg:4326 -lco ENCODING=UTF-8 $shapefile_proj $shapefile`)
        end
    end

    println("Rasterizing...")
    gdal_rasterize_path() do gdal_rasterize
        outfile = in_datafolder("natura2000.tif")
        isfile(outfile) && rm(outfile)
        # https://sdi.eea.europa.eu/catalogue/copernicus9129929/api/records/e40ca403-b81a-4ecb-b484-cade980e9a2f
        # SITETYPE contains "A", "B" or "C":  
        #   A: SPAs (Special Protection Areas - sites designated under the Birds Directive); 
        #   B: SCIs and SACs (Sites of Community Importance and Special Areas of Conservation - sites designated under the Habitats Directive); 
        #   C: where SPAs and SCIs/SACs boundaries are identical (sites designated under both directives).
        sql = "select unicode(SITETYPE)-64 as CODE,* from Natura2000_reprojected"
        @time run(`$gdal_rasterize -a CODE -ot Byte -tr 0.01 0.01 -te -32 28 34 70
                -co COMPRESS=LZW -dialect sqlite -sql $sql $shapefile_proj $outfile`)
    end
    readraster(in_datafolder("natura2000.tif"))
end

function rasterize_MIUU()
    shapefile = in_datafolder("MIUU vindkartering-100m-sweref99", "vindkartering 2011_100m.shp")
    shapefile_proj = in_datafolder("MIUU vindkartering-100m-sweref99", "MIUU_reprojected.shp")

    # Couldn't figure out how to make gdal_rasterize reproject on the fly, so...
    if !isfile(shapefile_proj)
        println("Reprojecting shapefile to EPSG:4326...")
        ogr2ogr_path() do ogr2ogr
            @time run(`$ogr2ogr -t_srs epsg:4326 -lco ENCODING=UTF-8 $shapefile_proj $shapefile`)
        end
    end

    println("Rasterizing...")
    gdal_rasterize_path() do gdal_rasterize
        outfile = in_datafolder("miuu_windatlas.tif")
        isfile(outfile) && rm(outfile)
        @time run(`$gdal_rasterize -a Z -ot Float32 -tr 0.01 0.01 -te 10.65 55.07 24.17 69.07
                -co COMPRESS=LZW $shapefile_proj $outfile`)
    end
    readraster(in_datafolder("miuu_windatlas.tif"))
end

function readfarms()
    df = DataFrame(CSV.File(in_datafolder("Windfarms_World_20240407.csv"); quotechar='\'', missingstring=["#ND", ""]))
    # ["ID (#ND = no data)", "Continent", "ISO code (Code ISO 3166.1)", "Country", "State code", "Area", "City", "Name", "2nd name",
    #     "Latitude (WGS84)", "Longitude (WGS84)", "Altitude/Depth (m)", "Location accuracy (Yes = accurate location)", "Offshore - Shore distance (km)",
    #     "Manufacturer", "Turbine", "Hub height (m)", "Number of turbines", "Total power (kW)", "Developer", "Operator", "Owner",
    #     "Commissioning date (Format: yyyy or yyyymm)", "Status", "Decommissioning date (Format: yyyy or yyyymm)", "Link", "Update"]
    newnames = [:id, :continent, :iso, :country, :state, :area, :city, :name, :name2, :lat, :lon, :altitude, :accurate_location, :offshore,
        :manufacturer, :turbine, :hubheight, :num_turbines, :capac, :developer, :operator, :owner, :startdate, :status, :enddate, :link, :updated]
    rename!(df, newnames)
    df[!, [:lat, :lon]] = [ismissing(x) ? missing : parse(Float64, replace(x, "," => ".")) for x in Array(df[!, [:lat, :lon]])]
    df.year = [ismissing(x) ? missing : parse(Int, x[1:4]) for x in df.startdate]
    df.endyear = [ismissing(x) ? missing : parse(Int, x[1:4]) for x in df.enddate]
    df.month = [ismissing(x) ? missing : (m = match(r".*/(\d+)", x)) === nothing ? missing : parse(Int, m[1]) for x in df.startdate]
    df.endmonth = [ismissing(x) ? missing : (m = match(r".*/(\d+)", x)) === nothing ? missing : parse(Int, m[1]) for x in df.enddate]
    df.capac .= round.(df.capac / 1e6, digits=3)
    df.accurate_location = (df.accurate_location .== "Yes")     # Yes or No, no missing data
    df.onshore = .!startswith.(df.offshore, "Yes")              # all these also have df.area=="Offshore", no missing data
    df.altitude = [ismissing(x) ? missing : (m = match(r"(\d+)/(\d+)", x)) === nothing ? round(Int, parse(Float64, x)) :
                            round(Int, mean(parse.(Int, [m[1], m[2]]))) for x in df.altitude]

    replacecountries = ["United-Kingdom" => "United Kingdom", "New-Zealand" => "New Zealand", "North Macedonia" => "Macedonia"]
    replace!(df.country, replacecountries...)
    df.accurate_location[ismissing.(df.lat) .&& df.accurate_location] .= false      # one weird record in Poland
    filter!(row -> row.status == "Production", df)
    select!(df, [newnames[1:13]; :onshore; newnames[15:22]; :year; :month; :status])
    return df
end

function add_gisdata_to_farms(df; optionlist...)
    regions, regionsEU, regionlist, gadm, subregionnames, lons, lats, res, lonlim, latlim, lonrange, latrange = get_europe_datasets()
    reg54country = Dict(i => NUTScountries[string(reg)[1:2]] for (i, reg) in enumerate(regionlist))

    # assume wind_speed_altitude = wind_class_altitude = 1001!
    windatlas = getwindatlas(100)[lonrange,latrange]
    options = WindOptions(merge(windoptions(), optionlist))                 
    onshoreclass, offshoreclass = makewindclasses(options, windatlas)

    invest_onoffshore_per_region_class_yearcode = zeros(2, 5, length(regionlist), 11) 

    df.oldreg54 .= 0
    df.windclass .= 0
    df.yearcode .= 232323

    updateprogress = Progress(nrow(df), 1)
    for row in eachrow(df)
        lon, lat = row.lon, row.lat
        next!(updateprogress)
        (lon < lonlim[1] || lon > lonlim[2] || lat < latlim[1] || lat > latlim[2]) && continue
        rasterindex = lonlat_index(regionsEU, lon, lat)
        if get(reg54country, row.reg54, "") != row.country
            regindexes = [i for (i, reg) in reg54country if reg == row.country] |> sort
            if !isempty(regindexes)
                ii = regions.>=regindexes[1] .&& regions.<=regindexes[end]
                if sum(ii) > 0
                    rr = GeoArray(regions[feature_transform(ii)], res, lonlim, latlim)
                    rasterindex = lonlat_index(rr, lon, lat)
                    row.oldreg54 = row.reg54
                    row.reg54 = rr[rasterindex]
                    get(reg54country, row.reg54, "") != row.country && error("Country mismatch")
                end
            end
        end
        row.windclass = (row.onshore ? onshoreclass[rasterindex] : offshoreclass[rasterindex])
        row.yearcode = round_yearcode(row.year)             # 1=missing, 2=1980, 3=1985, 4=1990, ..., 10=2020, 11=2025
        if row.reg54 < 30000 && !ismissing(row.capac)
            onoff = 2 - row.onshore                         # 1=onshore, 2=offshore
            invest_onoffshore_per_region_class_yearcode[onoff, row.windclass, row.reg54, row.yearcode] += row.capac
        end
    end
    println()
    df.reg .= [rr < 999 ? string(regionlist[rr]) : "" for rr in df.reg54]
    df.reg_guess .= [rr < 999 ? string(regionlist[rr]) : "" for rr in df.reg54_guess]
    
    return df, invest_onoffshore_per_region_class_yearcode
end

function get_europe_datasets()
    println("\nEUROPE 56!!!!!")
    regions, _, regionlist, lonrange, latrange = loadregions("Europe54_SEfix")
    res = 0.01
    res2 = res/2
    lons = (-180+res/2:res:180-res/2)[lonrange]         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    lonlim = (lons[1]-res2, lons[end]+res2)
    latlim = (lats[end]-res2, lats[1]+res2)
    gadm, subregionnames = read_gadm()
    gadm = gadm[lonrange, latrange]
    # @show lonlim, latlim
    regions = regions[feature_transform(regions.>0)]        # ensure regions go offshore (maybe unnecessary?)
    regionsEU = GeoArray(regions, res, lonlim, latlim)
    return regions, regionsEU, regionlist, gadm, subregionnames, lons, lats, res, lonlim, latlim, lonrange, latrange
end

# run time 1.5-2 minutes
function guess_locations(df; skipguess=false)
    regions, regionsEU, regionlist, gadm, subregionnames, lons, lats, res, lonlim, latlim, lonrange, latrange = get_europe_datasets()

    df.reg54 .= 232323
    df.reg54_guess .= 232323
    if !skipguess
        df.lon_guess .= 232323.0
        df.lat_guess .= 232323.0
        df.guesslevel .= 232323
    end
    regioncache = Dict{String, Tuple{Float64, Float64}}()

    updateprogress = Progress(nrow(df), 1)
    for row in eachrow(df)
        lon, lat = row.lon, row.lat
        if !ismissing(lon) && lon >= lonlim[1] && lon <= lonlim[2] && lat >= latlim[1] && lat <= latlim[2]
            rasterindex = lonlat_index(regionsEU, lon, lat)
            row.reg54 = regionsEU[rasterindex]
        end
        skipguess && continue
        lon, lat = guess_lonlat_from_windfarm_regions(row, gadm, subregionnames, lons, lats, regioncache)
        next!(updateprogress)
        isnan(lon) && continue
        rasterindex_guess = lonlat_index(regionsEU, lon, lat)
        row.lon_guess, row.lat_guess = lon, lat
        row.reg54_guess = regionsEU[rasterindex_guess]
    end
    return df
end

# Distributes investments with missing years over the other years in proportion to the sum of investments in those years
# invest: 2×5×54×11, onshore/offshore x wind class x region x yearcode
function distribute_investments_with_missing_years!(invest)
    # [vec(invest[:,:,1,:]) vec(sum(invest[:,:,2:end,:], dims=3))]
    inv_yearmissing = invest[:,:,:,1:1]
    inv_sum = sum(invest, dims=4)
    invmult = 1.0 .+ inv_yearmissing ./ inv_sum
    invmult[isnan.(invmult)] .= 1.0
    invest[:,:,:,2:end] .*= invmult
    invest[:,:,:,1] .= 0.0
    # inv_sum2 = sum(invest, dims=4)
    invest = round.(invest, digits=3)
    return invest
end

round_year5(x) = ismissing(x) ? missing : round(Int, x / 5) * 5
round_yearcode(x) = ismissing(x) ? 1 : round(Int, (x - 1970)/5)
decodeyear(y) = (y == 1) ? 1111 : 1970 + 5*y

function guess_lonlat_from_windfarm_regions(farmrow, gadm, subregionnames, lons, lats, regioncache)
    regnames = subregionnames[subregionnames[:, 1] .== farmrow.country, :]
    gadmlevel2, gadmlevel3 = unique(regnames[:,2]), unique(regnames[:,3])
    level2, level3 = "", ""
    m = ismissing(farmrow.area) ? nothing : match(r"(.*)\s*\((.+)\)", farmrow.area)
    if m !== nothing
        level2 = closeststring(m.captures[2], gadmlevel2)
        level3 = closeststring(m.captures[1], gadmlevel3)
        if isempty(level2) && isempty(level3)
            level2 = closeststring(m.captures[1], gadmlevel2)
            level3 = closeststring(farmrow.city, gadmlevel3)
        end
    end
    if m === nothing || isempty(level2) && isempty(level3)
        level2 = closeststring(farmrow.area, gadmlevel2)
        level3 = closeststring(farmrow.city, gadmlevel3)
    end

    # println("$level3, $level2")
    if !isempty(level3) && !isempty(level2)
        key = "level3_$(level2)_$level3"
        lon, lat = get(regioncache, key, (Inf,Inf))
        farmrow.guesslevel = 23
        isfinite(lon) && return lon, lat
        gadmindexes = findall(subregionnames[:,3] .== level3 .&& subregionnames[:,2] .== level2 .&& subregionnames[:,1] .== farmrow.country)
        if !isempty(gadmindexes)
            lon, lat = loop_gadm(gadm, gadmindexes, lons, lats)
            regioncache[key] = lon, lat
            return lon, lat
        end
    end
    
    if !isempty(level3)
        key = "level3_$level3"
        farmrow.guesslevel = 3
        lon, lat = get(regioncache, key, (Inf,Inf))
        isfinite(lon) && return lon, lat
        gadmindexes = findall(subregionnames[:,3] .== level3 .&& subregionnames[:,1] .== farmrow.country)
    elseif !isempty(level2)
        key = "level2_$level2"
        farmrow.guesslevel = 2
        lon, lat = get(regioncache, key, (Inf,Inf))
        isfinite(lon) && return lon, lat
        gadmindexes = findall(subregionnames[:,2] .== level2 .&& subregionnames[:,1] .== farmrow.country)
    else
        key = "country_$(farmrow.country)"
        farmrow.guesslevel = 1
        lon, lat = get(regioncache, key, (Inf,Inf))
        isfinite(lon) && return lon, lat
        gadmindexes = findall(subregionnames[:,1] .== farmrow.country)
    end

    # isempty(gadmindexes) && println("\n\n$level3, $level2")
    lon, lat = loop_gadm2(gadm, gadmindexes, lons, lats)
    regioncache[key] = lon, lat
    return lon, lat
end

function clean_windfarm_database()
    println("Reading raw global wind farm database...")
    df0 = readfarms()
    println("Guessing coordinates for all European wind farms from region & area names (ETA 2 minutes)...")
    df1 = guess_locations(df0[df0.continent .== "Europe", :])

    println("Filling missing coordinates with guessed ones...")
    mm = ismissing.(df1.lon)
    df1.lon[mm] .= df1.lon_guess[mm]
    df1.lat[mm] .= df1.lat_guess[mm]
    df1.reg54[mm] .= df1.reg54_guess[mm]
    df1.lon, df1.lat = coalesce.(df1.lon), coalesce.(df1.lat)
    df1.dist = sqrt.((df1.lon - df1.lon_guess).^2 .+ (df1.lat - df1.lat_guess).^2)
    df1.dist[df1.lon_guess .> 1000] .= 232323

    CSV.write(in_datafolder("Windfarms_Europe_20240407_CLEANED.csv"), df1)
    return nothing
end

function GISdata_for_ELLI_model(; plotmasks=true)
    # fixSEinEurope54()
    # makedistances("Europe54_SEfix")
    # createmaps("Europe54_SEfix")
    # clean_windfarm_database()

    df0 = CSV.File(in_datafolder("Windfarms_Europe_20240407_CLEANED.csv")) |> DataFrame
    select!(df0, [:country, :lat, :lon, :onshore, :capac, :year, :reg54, :reg54_guess])
    dfvb = read_vindbrukskollen()
    guess_locations(dfvb, skipguess=true)   # adds reg54 to dataframe
    select!(dfvb, [:country, :lat, :lon, :onshore, :capac, :year, :reg54, :reg54_guess])
    delete!(df0, df0.country .== "Sweden" .&& df0.onshore)

    println("Add wind classes from GIS, correct offshore countries and aggregate capacity by class and investment year...")
    df, invest = add_gisdata_to_farms(vcat(df0, dfvb))     # Can add GIS options here

    filter!(row -> row.reg54 > 0 && !ismissing(row.capac) && row.capac > 0, df)     # remove farms with missing or zero capacity 
    df.year5 .= round_year5.(df.year)

    sort!(df, [:reg54, :year5])
    gdf = groupby(df, [:reg54, :year5])
    # gdf_tot = combine(gdf, :capac => sum)

    gisregion="Europe54_SEfix"

    GISsolar(; gisregion, era_year=1991, grid_everywhere=true, plant_area=1.0, pvroof_area=1.0, plotmasks)
    GISwind(; gisregion, era_year=1991, grid_everywhere=true, area_onshore=1.0, area_offshore=1.0, plotmasks)
    predictdemand(; gisregion, sspscenario="ssp2-26", sspyear=2020, era_year=1991)

    GISsolar(; gisregion, era_year=1992, grid_everywhere=true, plant_area=1.0, pvroof_area=1.0, plotmasks=false)
    GISwind(; gisregion, era_year=1992, grid_everywhere=true, area_onshore=1.0, area_offshore=1.0, plotmasks=false)
    predictdemand(; gisregion, sspscenario="ssp2-26", sspyear=2020, era_year=1992)

    distribute_investments_with_missing_years!(invest)
    matlab2multinode(invest; gisregion, year=1991)
    matlab2multinode(invest; gisregion, year=1992)

    open(in_datafolder("output", "README_GISparameters_$(gisregion).txt"), "w") do f
        commands = """
        The GIS data in this folder was created on $(Dates.now()) using GlobalEnergyGIS commit 4c2fd4f.
        Below are the commands and parameters used to create the data (see GISdata_for_ELLI_model()):

        gisregion="Europe54_SEfix"

        GISsolar(; gisregion, era_year=1991, grid_everywhere=true, plant_area=1.0, pvroof_area=1.0, plotmasks)
        GISwind(; gisregion, era_year=1991, grid_everywhere=true, area_onshore=1.0, area_offshore=1.0, plotmasks)
        predictdemand(gisregion, sspscenario="ssp2-26", sspyear=2020, era_year=1991)

        GISsolar(; gisregion, era_year=1992, grid_everywhere=true, plant_area=1.0, pvroof_area=1.0, plotmasks=false)
        GISwind(; gisregion, era_year=1992, grid_everywhere=true, area_onshore=1.0, area_offshore=1.0, plotmasks=false)
        predictdemand(gisregion, sspscenario="ssp2-26", sspyear=2020, era_year=1992)

        distribute_investments_with_missing_years!(invest)
        matlab2multinode(invest; gisregion, year=1991)
        matlab2multinode(invest; gisregion, year=1992)
        """
        println(f, commands)
    end

    nothing
end

function swedish_capacity_diagnostic()
    df0 = CSV.File(in_datafolder("Windfarms_Europe_20240407_CLEANED.csv")) |> DataFrame
    df, invest = add_gisdata_to_farms(df0) 
    # fix_investments(invest)
    # regions, _, regionlist, lonrange, latrange = loadregions("Europe54")
    # df.regname = [r < 30000 ? regionlist[r] : Symbol() for r in df.reg54]
    df[.!ismissing.(df.reg54 .+ df.capac) .&& df.reg54 .>= 48 .&& df.reg54 .<= 51 .&& df.iso .!= "SE", [4; 6:8; 10; 11; 13:19; 22:25; 32; 33]] |> display
    df[.!ismissing.(df.reg54 .+ df.capac) .&& (df.reg54 .< 48 .|| df.reg54 .> 51) .&& df.iso .== "SE", [4; 6:8; 10; 11; 13:19; 22:25; 32; 33]] |> display
    df[df.name .== "Kriegers Flak" .|| df.name .== "Lillgrund", [4; 6:8; 10; 11; 13:19; 22:25; 32; 33]] |> display
    df
end

function loop_gadm(gadm, gadmindexes, lons, lats)
    sumlons, sumlats, n = 0.0, 0.0, 0
    low, hi = extrema(gadmindexes)
    Threads.@threads for i in eachindex(IndexCartesian(), gadm)
        g = gadm[i]
        if g >= low && g <= hi && g in gadmindexes
            sumlons += lons[i[1]]
            sumlats += lats[i[2]]
            n += 1
        end
    end
    meanlon, meanlat = sumlons/n, sumlats/n
    return meanlon, meanlat
end

function loop_gadm2(gadm, gadmindexes, lons, lats)
    gg = similar(gadm, Bool)
    low, hi = extrema(gadmindexes)
    Threads.@threads for i in eachindex(IndexCartesian(), gadm)
        g = gadm[i]
        if g >= low && g <= hi && g in gadmindexes
            gg[i] = true
        end
    end
    ggf = findall(gg)
    isempty(ggf) && return NaN, NaN
    lon = [lons[cc[1]] for cc in ggf]
    lat = [lats[cc[2]] for cc in ggf]
    return median(lon), median(lat)
end

similarity(x,y) = Levenshtein()(x, y)

function closeststring(s, targetstrings)
    (ismissing(s) || isempty(s)) && return ""
    dist, ndx = findmin(similarity.(targetstrings, s))
    return (dist/length(s) <= 0.3 ? targetstrings[ndx] : "")
end

function savefarms(winddir = "C:/Users/niclas/Downloads/wind data")
    wf = readfarms(winddir)
    CSV.write("$winddir/windfarms.csv", wf)
end

function read_irena(winddir = "C:/Users/niclas/Downloads/wind data")
    df = DataFrame(CSV.File("$winddir/IRENA ELECCAP_20220411-085642.csv", header=3, missingstring=".."))
    # "Country/area", "Technology", "Grid connection", "Year", "Installed electricity capacity by country/area (MW)"
    newnames = [:country, :onshore, :gridconnected, :year, :capac]
    rename!(df, newnames)
    df.onshore = (df.onshore .== "Onshore wind energy")
    df.gridconnected = (df.gridconnected .== "On-grid")
    replacecountries = ["United Kingdom of Great Britain and Northern Ireland" => "United Kingdom", "United States of America" => "USA"]
    replace!(df.country, replacecountries...)
    # rows = df.gridconnected .&& df.year .== 2020
    # onshore2020 = Dict(row.country => row.capac for row in eachrow(df[rows .&& df.onshore, :]) if !ismissing(row.capac))
    # offshore2020 = Dict(row.country => row.capac for row in eachrow(df[rows .&& .!df.onshore, :]) if !ismissing(row.capac))
    return df
end

# Use bidding zone shapefiles to adjust borders of SE1-4
function fixSEinEurope54()
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions("Europe54")
    landcover = JLD.load(in_datafolder("landcover.jld"), "landcover")
    company, dfcompany = getcompanyraster()
    landcover = landcover[lonrange, latrange]
    company = company[lonrange, latrange]
    company = company[feature_transform(company.>0)]    # expand company map to fill in holes and artifacts in the shapefile
    iii = regions.>=48 .&& regions.<=51 .&& company.>0
    regions[iii] .= 51 .- company[iii] .+ 1
    territory = regions[feature_transform(regions.>0)]
    offshoreregions = territory .* (landcover .== 0)
    regionname = "Europe54_SEfix"    
    JLD.save(in_datafolder("regions_$regionname.jld"), "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

function getcompanydata()
    dfcompany = GDF.read("C:/Griddata/Elnätsområden Therese/omraden.shp")
    dfcompany.bolag[276] = "Hedemorahyttorna"
    dfcompany.snitt = parse.(Int, dfcompany.snitt)
    disallowmissing!(dfcompany)
    return dfcompany
end

"""Rasterize Swedish power company areas using center points of our 1 km grid squares.
    Return (company, dfcompany), where company is a Raster and dfcompany a GeoDataFrame of company areas."""
function getcompanyraster()
    # Shapefile of Swedish power companies (purchased by Therese)
    # Is this the source?   https://www.natomraden.se/
    dfcompany = getcompanydata()
    source = GeoFormatTypes.ProjString("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")    # SWEREF99 TM = EPSG:3006
    dest = GeoFormatTypes.ProjString("+proj=longlat +datum=WGS84 +no_defs")
    dfcompany = GDF.reproject(dfcompany, source, dest)
    # dfcompany = GDF.reproject(dfcompany, EPSG(3006), EPSG(4326))
    # company = rasterizeSWEREF99(dfcompany, :geometry)
    company = reverse(rasterize_global(dfcompany, :geometry, :snitt), dims=2)
    return company, dfcompany
end

"""Rasterize a SWEREF99 GeoDataFrame using center points of our 1 km grid squares.
    Return (swe, dfdeso), where swe is a Raster and dfdeso a GeoDataFrame of DeSO areas."""
function rasterizeSWEREF99(gdf, column)
    # Build a raster by sampling the vector data every 1000 m in the *center* of our 1 km grid squares.
    # The SWEREF99 cell ID ("Ruta") refers to the lower left corner of each cell, not the center.
    # (outer coordinate limits: x: 181896.33 - 1086312.94, y: 6090353.78 - 7689478.31)
    rasterize_geovector(gdf[!, column], (181500,1086500), (6089500,7689500), 1000, EPSG(3006))
end

# function rasterize_global(gdf, column)
#     println("\nRasterizing GADM shapefile for global administrative areas (1-10 minute run time)...")
#     # shapefile = in_datafolder("gadm36", "gadm36.shp")
#     outfile = in_datafolder("hyunkyo_gbg.tif")
#     options = "-a id -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
#     @time rasterize(geojson, outfile, split(options, ' '))
# end

function rasterize_global(gdf, geomcolumn, datacolumn)
    rasterize_geovector(gdf[!, geomcolumn], gdf[!, datacolumn], (-179.995,179.995), (-89.995,89.995), 0.01, EPSG(4326))
end

"""rasterize_geovector(gv, xlim, ylim, interval, crs; T=UInt32, missingval=UInt32(0))

    Rasterize a geovector `gv` (i.e. the geometry column containing polygons of a GeoDataFrame)
    using the values `vv` between coordinate limits `xlim` and `ylim` (i.e. Tuples in format (xmin, xmax)),
    spaced by `interval` and using coordinate system `crs`. Return a Raster."""
function rasterize_geovector(gv, vv, xlim, ylim, interval, crs; T=UInt32, missingval=UInt32(0))
    # x, y = values(ext)
    dimz = Rasters.X(xlim[1]:interval:xlim[2]; mode=Rasters.Projected(; sampling=Rasters.Points(), crs)),
        Rasters.Y(ylim[1]:interval:ylim[2]; mode=Rasters.Projected(; sampling=Rasters.Points(), crs))
    raster = Rasters.Raster(zeros(T, dimz); missingval)
    for (g, v) in zip(gv, vv)
        Rasters.rasterize!(raster, g, fill=v)
    end
    return raster
end
rasterize_geovector(gv, xlim, ylim, interval, crs; T=UInt32, missingval=UInt32(0)) =
    rasterize_geovector(gv, 1:length(gv), xlim, ylim, interval, crs; T, missingval)

function openmap(df::DataFrame, turbinenumber::Int)
    openmap(df, turbinenumber, :google)
    openmap(df, turbinenumber, :bing)
end
    
function openmap(df::DataFrame, turbinenumber::Int, source::Symbol, use_guess=false)
    if use_guess
        openmap(df.lon_guess[turbinenumber], df.lat_guess[turbinenumber], source)
    else
        openmap(df.lon[turbinenumber], df.lat[turbinenumber], source)
    end
    return df[turbinenumber, :]
end

function openmap(lon::Real, lat::Real, source=:google)
    # Extra quotes to avoid errors with special chars ? and &:
    # https://superuser.com/questions/36728/can-i-launch-urls-from-command-line-in-windows
    # https://discourse.julialang.org/t/quoting-special-characters-of-a-url-in-cmd-objects-on-windows/44324

    if source == :google
        # url = "https://www.google.com/maps/@?api=1\"&\"map_action=map\"&\"basemap=satellite\"&\"center=$lat%2C$lon\"&\"zoom=18"   # satellite map but no pin
        # url = "http://maps.google.com/maps?t=k\"&\"q=loc:$lat+$lon"   # used to do both but no longer places pin
        # https://stackoverflow.com/questions/47038116/google-maps-url-with-pushpin-and-satellite-basemap
        # https://stackoverflow.com/questions/60219254/show-location-marker-in-new-browser-window

        url = "https://www.google.com/maps/search/?api=1\"&\"query=$lat%2C$lon\"&\"basemap=satellite\"&\"zoom=18"   # map with pin
        # url = "https://www.google.com/maps/@?api=1\"&\"map_action=map\"&\"basemap=satellite\"&\"center=$lat%2C$lon\"&\"zoom=18" # no pin
    else
        url = "https://bing.com/maps/default.aspx?cp=$lat~$lon\"&\"lvl=18\"&\"style=a\"&\"sp=point.$(lat)_$(lon)_"
    end
    c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
    run(c)
end

# Use Google Map directions to show two sets of points
function openmapdir(lon1, lat1, lon2, lat2)
    url = "https://www.google.com/maps/dir/$lat1,$lon1/$lat2,$lon2/"
    c = Cmd(`cmd /c start \"\" $url`, windows_verbatim=true)
    run(c)
end

function openmapdir(df::DataFrame, turbinenumber::Int)
    openmapdir(df.lon[turbinenumber], df.lat[turbinenumber], df.lon_guess[turbinenumber], df.lat_guess[turbinenumber])
end

function read_vindbrukskollen()
    # Länsstyrelsen: Vindbrukskollen, https://vbk.lansstyrelsen.se/  (click "Excel-export")
    df = DataFrame(CSV.File(in_datafolder("Vindbrukskollen land 2025-01-28.csv")))
    select!(df, ["Status", "Placering", "E-Koordinat", "N-Koordinat", "Navhöjd (m)", "Rotordiameter (m)", "Maxeffekt (MW)", "Uppfört", "Fabrikat", "Modell"])
    rename!(df, [:status, :type, :lon, :lat, :hubheight, :rotordiam, :capac, :year, :brand, :model])
    delete!(df, df.status .!= "Uppfört")
    select!(df, Not(:status))
    df.onshore .= true      # startswith.(df.type, "Land")     # all turbines are marked "Land" or "Vatten"
    df.country .= "Sweden"
    select!(df, Not(:type))
    delete!(df, ismissing.(df.capac) .|| df.capac .== 0)
    df.capac .= coalesce.(df.capac / 1000)                  # 
    df.year .= max.(1980, year.(df.year))
    df.model = strip.(coalesce.(df.brand, "")) .* " " .* strip.(coalesce.(df.model, ""))
    df.model = [m == " " ? missing : m for m in df.model]
    select!(df, Not(:brand))

    # df2 = DataFrame(CSV.File(in_datafolder("Vindbrukskollen hav 2025-01-28.csv")))
    # select!(df2, ["Projektstatus", "Parken uppförd", "Uppfört antal verk", "Installerad effekt (MW)", "Elområde", "Län", "Kommun"])
    # rename!(df2, [:status, :year, :nturbines, :capac, :zone, :region, :munic])
    # delete!(df2, df2.status .!= "Uppförd")

    # https://www.lantmateriet.se/sv/Kartor-och-geografisk-information/gps-geodesi-och-swepos/referenssystem/tvadimensionella-system/sweref-99-projektioner/
    ec, nc = df.lon, df.lat
    source = "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"     # SWEREF99 TM = EPSG:3006
    dest = "+proj=longlat +datum=WGS84 +no_defs"
    trans = Proj.Transformation(source, dest)
    len = size(df,1)
    lon, lat = zeros(len), zeros(len)
    for i = 1:len
        lon[i], lat[i] = trans(ec[i], nc[i])
    end
    df.lon = lon
    df.lat = lat
    return df
end
