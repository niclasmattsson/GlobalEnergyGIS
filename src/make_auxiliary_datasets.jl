export rasterize_datasets, create_scenario_datasets, cleanup_datasets, makeprotected, savelandcover,
        createGDP, creategridaccess, getpopulation, getwindatlas

# cleanup options: :none, :limited, :all
function rasterize_datasets(; cleanup=:all)
    rasterize_GADM()
    rasterize_NUTS()
    rasterize_protected()
    makeprotected()
    downscale_landcover()
    savelandcover()
    upscale_topography()
    saveregions_global()
    cleanup_datasets(cleanup=cleanup)
end

function create_scenario_datasets(scen, year)
    println("\nCreate population dataset...")
    downscale_population(scen, year)
    println("\nCreate GDP dataset...")
    createGDP(scen, year)
    println("\nCreate grid access dataset...")
    creategridaccess(scen, year)
end

# cleanup options: :none, :limited, :all
function cleanup_datasets(; cleanup=:all)
    cleanup == :none && return
    datafolder = getconfig("datafolder")
    rm(joinpath(datafolder, "protected_raster.tif"), force=true)
    rm(joinpath(datafolder, "protectedfields.csv"), force=true)
    rm(joinpath(datafolder, "landcover.tif"), force=true)
    rm(joinpath(datafolder, "topography.tif"), force=true)
    if cleanup == :all
        rm(joinpath(datafolder, "Landcover - USGS MODIS.tif"), force=true)
        rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"), force=true)
        rm(joinpath(datafolder, "gadm36"), force=true, recursive=true)
        rm(joinpath(datafolder, "nuts2016-level3"), force=true, recursive=true)
        rm(joinpath(datafolder, "WDPA_Feb2020"), force=true, recursive=true)
    end
end

function rasterize_GADM()
    println("Rasterizing GADM shapefile for global administrative areas (2-10 minute run time)...")
    ENV["PROJ_LIB"] = abspath(dirname(pathof(GDAL)), "../deps/usr/share/proj") # hack to fix GDAL/PROJ path
    oldpath = pwd()
    cd(Sys.BINDIR)  # Another hack to access libs in BINDIR, e.g. libstdc++-6.dll
    datafolder = getconfig("datafolder")
    shapefile = joinpath(datafolder, "gadm36", "gadm36.shp")
    outfile = joinpath(datafolder, "gadm.tif")
    options = "-a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize(shapefile, outfile, split(options, ' '))
 
    println("Creating .csv file for regional index and name lookup...")
    sql = "select uid,name_0,name_1,name_2 from gadm36"
    # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    outfile = joinpath(datafolder, "gadmfields.csv")
    ogr2ogr = GDAL.gdal.ogr2ogr_path
    @time run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
    cd(oldpath)
    nothing
end

function rasterize_NUTS()
    println("Rasterizing NUTS shapefile for European administrative areas...")
    ENV["PROJ_LIB"] = abspath(dirname(pathof(GDAL)), "../deps/usr/share/proj") # hack to fix GDAL/PROJ path
    oldpath = pwd()
    cd(Sys.BINDIR)  # Another hack to access libs in BINDIR, e.g. libstdc++-6.dll
    datafolder = getconfig("datafolder")
    name = "NUTS_RG_01M_2016_4326_LEVL_3"
    shapefile = joinpath(datafolder, "nuts2016-level3", "$name.shp")
    outfile = joinpath(datafolder, "nuts.tif")
    options = "-a ROWID -ot Int16 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW -dialect SQlite"
    sql = "select ROWID+1 AS ROWID,* from $name"
    @time rasterize(shapefile, outfile, split(options, ' '), sql=sql)
 
    println("Creating .csv file for regional index and name lookup...")
    outfile = joinpath(datafolder, "nutsfields.csv")
    sql = "select ROWID+1 AS ROWID,* from $name"
    ogr2ogr = GDAL.gdal.ogr2ogr_path
    @time run(`$ogr2ogr -f CSV $outfile -dialect SQlite -sql $sql $shapefile`)
    cd(oldpath)
    nothing
end

function read_gadm()
    println("Reading GADM rasters...")
    datafolder = getconfig("datafolder")
    gadmfields = readdlm(joinpath(datafolder, "gadmfields.csv"), ',', header=true)[1]
    imax = maximum(gadmfields[:,1])
    subregionnames = fill("", (imax,3))
    subregionnames[gadmfields[:,1],:] = string.(gadmfields[:,2:4])
    gadm = readraster(joinpath(datafolder, "gadm.tif"))
    return gadm, subregionnames
end

function read_nuts()
    println("Reading NUTS rasters...")
    datafolder = getconfig("datafolder")
    nutsfields = readdlm(joinpath(datafolder, "nutsfields.csv"), ',', header=true)[1]
    imax = maximum(nutsfields[:,1])
    subregionnames = nutsfields[:,3]    # indexes of NUTS regions are in order 1:2016, let's use that 
    nuts = readraster(joinpath(datafolder, "nuts.tif"))
    return nuts, subregionnames
end

function rasterize_protected()
    println("Rasterizing WDPA shapefile for protected areas (run time 10 minutes - 2 hours)...")
    ENV["PROJ_LIB"] = abspath(dirname(pathof(GDAL)), "../deps/usr/share/proj") # hack to fix GDAL/PROJ path
    oldpath = pwd()
    cd(Sys.BINDIR)  # Another hack to access libs in BINDIR, e.g. libstdc++-6.dll
    datafolder = getconfig("datafolder")
    shapefile = joinpath(datafolder, "WDPA_Feb2020", "WDPA_Feb2020-shapefile-polygons.shp")
    sql = "select FID from \"WDPA_Feb2020-shapefile-polygons\""
    outfile = joinpath(datafolder, "protected_raster.tif")
    options = "-a FID -a_nodata -1 -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    gdal_rasterize = GDAL.gdal.gdal_rasterize_path
    @time run(`$gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)

    println("Creating .csv file for WDPA index and name lookup...")
    sql = "select FID,IUCN_CAT from \"WDPA_Feb2020-shapefile-polygons\""
    outfile = joinpath(datafolder, "protectedfields.csv")
    ogr2ogr = GDAL.gdal.ogr2ogr_path
    @time run(`$ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
    cd(oldpath)
    nothing
end

function makeprotected()
    println("Reading protected area rasters...")
    datafolder = getconfig("datafolder")
    protectedfields = readdlm(joinpath(datafolder, "protectedfields.csv"), ',', header=true)[1]
    IUCNcodes = ["Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Reported", "Not Applicable", "Not Assigned"]
    IUCNlookup = Dict(c => i for (i,c) in enumerate(IUCNcodes))
    protected0 = readraster(joinpath(datafolder, "protected_raster.tif"))

    println("Converting indexes to protected area types (2-3 minutes)...")
    protected = similar(protected0, UInt8)
    # could replace loop with:  map!(p -> p == -1 ? 0 : IUCNlookup[protectedfields[p+1,2], protected, protected0)
    # alternatively             map!(p -> ifelse(p == -1, 0, IUCNlookup[protectedfields[p+1,2]), protected, protected0)
    for (i, p) in enumerate(protected0)
        protected[i] = (p == -1) ? 0 : IUCNlookup[protectedfields[p+1,2]]
    end

    println("Saving protected area dataset...")
    JLD.save(joinpath(datafolder, "protected.jld"), "protected", protected, compress=true)
end

function resample(infile::String, outfile::String, options::Vector{<:AbstractString})
    ENV["PROJ_LIB"] = abspath(dirname(pathof(GDAL)), "../deps/usr/share/proj") # hack to fix GDAL/PROJ path
    oldpath = pwd()
    cd(Sys.BINDIR)  # Another hack to access libs in BINDIR, e.g. libstdc++-6.dll
    gdal_translate = GDAL.gdal.gdal_translate_path
    @time run(`$gdal_translate $options -co COMPRESS=LZW $infile $outfile`)
    cd(oldpath)
end

function downscale_landcover()
    println("Downscaling landcover dataset (5-10 minutes)...")
    datafolder = getconfig("datafolder")
    infile = joinpath(datafolder, "Landcover - USGS MODIS.tif")
    options = "-r mode -ot Byte -tr 0.01 0.01"
    resample(infile, joinpath(datafolder, "landcover.tif"), split(options, ' '))
    nothing
end

function savelandcover()
    println("Converting landcover dataset from TIFF to JLD...")
    datafolder = getconfig("datafolder")
    landcover = readraster(joinpath(datafolder, "landcover.tif"))
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
    JLD.save(joinpath(datafolder, "landcover.jld"), "landcover", landcover, "landtypes", landtypes,
                "landcolors", landcolors, compress=true)
end

function upscale_topography()
    println("Upscaling topography dataset...")
    datafolder = getconfig("datafolder")
    infile = joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif")
    options = "-r cubicspline -tr 0.01 0.01"
    outfile = joinpath(datafolder, "topography.tif")
    resample(infile, outfile, split(options, ' '))
    println("Reading new topography raster...")
    topography = readraster(outfile)
    println("Saving topography dataset...")
    JLD.save(joinpath(datafolder, "topography.jld"), "topography", topography, compress=true)
end

# gettopography() = readraster("topography.tif")

function downscale_population(scen, year)
    scen = lowercase(scen)
    println("Reading population dataset...")
    datafolder = getconfig("datafolder")
    filename = joinpath(datafolder, "SSP2_1km", "$(scen)_total_$year.nc4")
    dataset = Dataset(filename)
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
    JLD.save(joinpath(datafolder, "population_$(scen)_$year.jld"), "population", newpop, compress=true)

    rm(temptiff, force=true)
    rm(temptiff2, force=true)
end

function getpopulation(scen, year)
    datafolder = getconfig("datafolder")
    JLD.load(joinpath(datafolder, "population_$(scen)_$year.jld"), "population")
end

function createGDP(scen, year)
    scen = lowercase(scen)
    scennum = scen[end]
    println("Reading low resolution population and GDP datasets...")
    datafolder = getconfig("datafolder")
    folder = joinpath(datafolder, "global_population_and_gdp")
    pop, extent = readraster(joinpath(folder, "p$(scennum)_$year.tif"), :getextent) # million people
    gdp = readraster(joinpath(folder, "g$(scennum)_$year.tif"))    # billion USD(2005), PPP

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
    JLD.save(joinpath(datafolder, "gdp_$(scen)_$year.jld"), "gdp", gdphigh, compress=true)
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
    datafolder = getconfig("datafolder")
    gdp = JLD.load(joinpath(datafolder, "gdp_$(scen)_$year.jld"), "gdp")
    res = 360/size(gdp,1)

    disk = diskfilterkernel(1/6/res)                        # filter radius = 1/6 degrees
    gridaccess = gridsplit(gdp .> 100_000, x -> imfilter(x, disk), Float32)
    # gridaccess = Float32.(imfilter(gdp .> 100_000, disk))   # only "high" income cells included (100 kUSD/cell), cell size = 1x1 km          
    println("\nCompressing...")
    selfmap!(x -> ifelse(x<1e-6, 0, x), gridaccess)         # force small values to zero to reduce dataset size
    println("Saving high resolution grid access...")
    JLD.save(joinpath(datafolder, "gridaccess_$(scen)_$year.jld"), "gridaccess", gridaccess, compress=true)

    # maybe better:
    # loop through countries, index all pixels into vector, sort by GDP, use electrification to assign grid access
end

function getwindatlas()
    # filename = "D:/datasets/Global Wind Atlas v3.0/gwa3_250_wind-speed_100m.tif"     # v3.0 (lon extent [-180.3, 180.3], strangely)
    # filename = "D:/datasets/Global Wind Atlas v2.3/global_ws.tif"     # v2.3 (lon extent [-180.3, 180.3], strangely)
    datafolder = getconfig("datafolder")
    filename = joinpath(datafolder, "Global Wind Atlas v1 - 100m wind speed.tif")   # v1.0
    windatlas = readraster(filename, :extend_to_full_globe)[1]
    clamp!(windatlas, 0, 25)
end

# @time q,cc = readraster("D:/datasets/Global Wind Atlas v3.0/gwa3_250_wind-speed_100m.tif", :none, 1);
# size(q,1)/(cc[3]-cc[1])
# cc[1]+132/400
