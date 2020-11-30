export rasterize_datasets, create_scenario_datasets, cleanup_datasets, makeprotected, savelandcover,
        createGDP, creategridaccess, getpopulation, getwindatlas

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
    println("\nRasterizing GADM shapefile for global administrative areas (2-10 minute run time)...")
    shapefile = in_datafolder("gadm36", "gadm36.shp")
    outfile = in_datafolder("gadm.tif")
    options = "-a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize(shapefile, outfile, split(options, ' '))
 
    println("Creating .csv file for regional index and name lookup...")
    sql = "select uid,name_0,name_1,name_2 from gadm36"
    # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    outfile = in_datafolder("gadmfields.csv")
    gdal_utility("ogr2ogr") do ogr2ogr
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
    gdal_utility("ogr2ogr") do ogr2ogr
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
    println("\nRasterizing three WDPA shapefiles for protected areas (total run time 10 minutes - 2 hours)...")

    for i = 0:2
        println("\nFile $(i+1)/3:")
        shapefile = in_datafolder("WDPA", "WDPA-shapefile$i", "WDPA-shapefile-polygons.shp")

        println("Rasterizing...")
        gdal_utility("gdal_rasterize") do gdal_rasterize
            outfile = in_datafolder("WDPA", "protected_raster$i.tif")
            sql = "select FID from \"WDPA-shapefile-polygons\""
            options = "-a FID -a_nodata -1 -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
            @time run(`$gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)
        end

        println("Creating .csv file for WDPA index and name lookup...")
        gdal_utility("ogr2ogr") do ogr2ogr
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
    for (i, p) in enumerate(protected0)
        protected[i] = (p == -1) ? 0 : IUCNlookup[protectedfields[p+1,2]]
    end
    println("Saving...")
    # JLD.save(in_datafolder("protected$n.jld"), "protected", protected, compress=true)
    JLD.save(in_datafolder("WDPA", "protected$n.jld"), "protected", protected)
end

function rasterize_timezones()
    println("Rasterizing shapefile of time zones...")
    shapefile = in_datafolder("timezones-with-oceans.shapefile", "dist", "combined-shapefile-with-oceans.shp")
    sql = "select FID+1 as FID from \"combined-shapefile-with-oceans\""
    outfile = in_datafolder("timezones.tif")
    options = "-a FID -a_nodata 0 -ot Int16 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    gdal_utility("gdal_rasterize") do gdal_rasterize
        @time run(`$gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)
    end

    println("Creating .csv file for time zone index and name lookup...")
    sql = "select FID+1 as FID,tzid from \"combined-shapefile-with-oceans\""
    outfile = in_datafolder("timezone_names.csv")
    gdal_utility("ogr2ogr") do ogr2ogr
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
    gdal_utility("gdal_translate") do gdal_translate
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

function getwindatlas()
    # filename = in_datafolder("gwa3_250_wind-speed_100m.tif") # v3.0 (lon extent [-180.3, 180.3], strangely)
    # filename = in_datafolder("global_ws.tif") # v2.3 (lon extent [-180.3, 180.3], strangely)
    # filename = in_datafolder("Global Wind Atlas v1 - 100m wind speed.tif")   # v1.0
    # windatlas = readraster(filename, :extend_to_full_globe)[1]
    filename = in_datafolder("Global Wind Atlas v3 - 100m wind speed.tif")   # v3.0
    windatlas = readraster(filename)
    clamp!(windatlas, 0, 25)
end

# Convert the Global Wind Atlas 3.0 dataset from 250 m to 1 km resolution.
# This reduces file size from 13 GB to 1 GB. Interpolate using cubic splines. 
# Also change its weird lon-lat extents to standard [-180,-90] - [180, 90].
function convert_windatlas3()
    infile = in_datafolder("gwa3_250_wind-speed_100m.tif")
    gdal_utility("gdalinfo") do gdalinfo
        run(`$gdalinfo $infile`)
    end
    println("\n")
    outfile = in_datafolder("Global Wind Atlas v3 - 100m wind speed.tif")
    options = split("-r cubicspline -te -180 -90 180 90 -tr 0.01 0.01", ' ')
    gdal_utility("gdalwarp") do gdalwarp
        @time run(`$gdalwarp $options -co COMPRESS=LZW $infile $outfile`)
    end
end
