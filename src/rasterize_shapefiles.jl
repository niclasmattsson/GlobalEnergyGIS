import GDAL
using ArchGDAL
const AG = ArchGDAL

export rasterize, readraster, saveTIFF, GADM, makeregions, makeregions_nuts, makeoffshoreregions, makeprotected,
        savelandcover, getpopulation, createGDP, creategridaccess, getwindatlas, saveregions, loadregions, read3draster
        # Node, findchild, printnode, print_tree, Leaves

function rasterize_AG(infile::String, outfile::String, options::Vector{<:AbstractString})
    AG.registerdrivers() do
        AG.read(infile) do dataset
            GDAL.close(GDAL.rasterize(
                outfile,
                Ptr{GDAL.GDALDatasetH}(C_NULL),
                dataset.ptr,
                GDAL.rasterizeoptionsnew(options, C_NULL), C_NULL))
        end
    end
end

# works (creates the TIFF and saves it) but then crashes
# registerdrivers() should soon be unnecessary, see https://github.com/yeesian/ArchGDAL.jl/pull/76
function rasterize_AG2(infile::String, outfile::String, options::Vector{<:AbstractString})
    AG.registerdrivers() do
        AG.read(infile) do dataset
            AG.unsafe_gdalrasterize(dataset, options, dest=outfile)
        end
    end
end

# uses the command line version instead (gdal_rasterize)
# significantly faster for some reason, also give a simple progress indication
function rasterize(infile::String, outfile::String, options::Vector{<:AbstractString}; sql::String="")
    GDAL_BINPATH = joinpath(dirname(pathof(GDAL)), "../deps/usr/bin")
    println("Add GDAL_BINPATH to the gdal_rasterize commands before release!")
    if isempty(sql)
        run(`gdal_rasterize $options $infile $outfile`)
    else
        run(`gdal_rasterize $options -sql $sql $infile $outfile`)
    end        
end

function getextent(geotransform::Vector{Float64}, rastersize::Tuple{Int,Int})
    @assert length(geotransform) == 6 "A GeoTransform vector must have 6 elements."
    left, xres, _, top, _, yres = geotransform
    width, height = rastersize
    bottom, right = top+yres*height, left+xres*width
    return [left, bottom, right, top]
end

function read3draster(infile::String)
    ArchGDAL.registerdrivers() do
        ArchGDAL.read(infile) do dataset
            ArchGDAL.read(dataset)
        end
    end
end

function readraster(infile::String, extentflag::Symbol, dim::Int=1)
    local raster, geotransform
    ArchGDAL.registerdrivers() do
        raster = ArchGDAL.read(infile) do dataset
            # display(ArchGDAL.getproj(dataset))
            geotransform = ArchGDAL.getgeotransform(dataset)
            ArchGDAL.read(dataset)[:,:,dim]
        end
    end
    coordextent = getextent(geotransform, size(raster))
    if extentflag == :extend_to_full_globe
        left, bottom, right, top = coordextent
        xres, yres = geotransform[2], geotransform[6]
        newwidth, newheight = round.(Int, (360/xres, -180/yres))
        xindexes = 1+round(Int, (left-(-180))/xres):newwidth-round(Int, (right-180)/xres)
        yindexes = 1+round(Int, (top-90)/yres):newheight+round(Int, (bottom-(-90))/yres)
        adjusted = zeros(eltype(raster), (newwidth, newheight))
        adjusted[xindexes, yindexes] = raster
        return adjusted, coordextent
    else    # extentflag == :getextent
        return raster, coordextent
    end
end

readraster(infile::String, dim::Int=1) = readraster(infile, :none, dim)[1]

function saveTIFF(x::AbstractMatrix, filename::String, extent::Vector{Float64})
    wkt_string = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]"
    ArchGDAL.registerdrivers() do
        width, height = size(x)
        xres = (extent[3]-extent[1])/width
        yres = (extent[4]-extent[2])/height
        raster = AG.unsafe_create(
            filename,
            AG.getdriver("GTiff"),
            width = width,
            height = height,
            nbands = 1,
            dtype = eltype(x),
            options = ["COMPRESS=LZW"]
        )
        ## assign the projection and transformation parameters
        AG.setgeotransform!(raster, [extent[1], xres, 0, extent[4], 0, -yres])
        AG.setproj!(raster, wkt_string)
        
        ## write the raster    
        AG.write!(
            raster,
            x,      # image to "burn" into the raster
            1,      # update band 1
        )
        AG.destroy(raster)
    end
    nothing
end

saveTIFF(x::AbstractMatrix, filename::String) = saveTIFF(x, filename, [-180.0, -90.0, 180.0, 90.0])

function rasterize_GADM()
    println("Rasterizing global shapefile...")
    shapefile = "C:/Stuff/Datasets/gadm36/gadm36.shp"
    outfile = "gadm.tif"
    options = "-a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize(shapefile, outfile, split(options, ' '))
 
    println("Creating .csv file for regional index and name lookup...")
    sql = "select uid,name_0,name_1,name_2 from gadm36"
    # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    outfile = "gadmfields.csv"
    @time run(`ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
end

function rasterize_NUTS()
    println("Rasterizing global shapefile...")
    name = "NUTS_RG_60M_2016_4326_LEVL_3"
    shapefile = "C:/Stuff/Datasets/NUTS regions (nuts-2016-60m)/$name.shp"
    outfile = "nuts.tif"
    options = "-a ROWID -ot Int16 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW -dialect SQlite"
    sql = "select ROWID+1 AS ROWID,* from $name"
    @time rasterize(shapefile, outfile, split(options, ' '), sql=sql)
 
    println("Creating .csv file for regional index and name lookup...")
    outfile = "nutsfields.csv"
    sql = "select ROWID+1 AS ROWID,* from $name"
    @time run(`ogr2ogr -f CSV $outfile -dialect SQlite -sql $sql $shapefile`)
end

function saveregions(regionname, regiondefinitionarray; crop=true)
    land = JLD.load("landcover.jld", "landcover")
    saveregions(regionname, regiondefinitionarray, land, crop)
end

function saveregions(regionname, regiondefinitionarray, landcover, crop)
    regions = makeregions(regiondefinitionarray) .* (landcover.>0)
    if crop
        # get indexes of the bounding box containing onshore region data with 3 degrees of padding
        lonrange, latrange = getbboxranges(regions, round(Int, 3/0.01))
    else
        nlon, nlat = size(regions)
        lonrange, latrange = 1:nlon, 1:nlat
    end
    regions = regions[lonrange, latrange]
    offshoreregions = makeoffshoreregions(regions, landcover[lonrange, latrange])
    regionlist = Symbol.(regiondefinitionarray[:,1])

    println("\nSaving regions and offshoreregions...")
    JLD.save("regions_$regionname.jld", "regions", regions, "offshoreregions", offshoreregions,
                "regionlist", regionlist, "lonrange", lonrange, "latrange", latrange, compress=true)
end

function saveregions_global_gadm0()
    g = readdlm("gadmfields.csv", ',', skipstart=1)
    gadm0 = unique(string.(g[:,2]))
    regiondefinitionarray = [gadm0 GADM.(gadm0)]
    saveregions("Global_GADM0", regiondefinitionarray, crop=false)
end

function loadregions(regionname)
    jldopen("regions_$regionname.jld", "r") do file
        return read(file, "regions"), read(file, "offshoreregions"), read(file, "regionlist"),
                    read(file, "lonrange"), read(file, "latrange")
    end
end

function makeregions_nuts(regiondefinitionarray)
    println("Reading NUTS rasters...")
    nutsfields = readdlm("nutsfields.csv", ',', header=true)[1]
    imax = maximum(nutsfields[:,1])
    subregionnames = nutsfields[:,3]    # indexes of NUTS regions are in order 1:2016, let's use that 
    nuts = readraster("nuts.tif")

    regionnames = regiondefinitionarray[:,1]
    regiondefinitions = [regdef for regdef in regiondefinitionarray[:,2]]

    println("Making region index matrix...")
    return makeregions_main_nuts(nuts, subregionnames, regiondefinitions)
end

function makeregions_main_nuts(nuts, subregionnames, regiondefinitions)
    regionlookup = Dict(r => i for (i,tup) in enumerate(regiondefinitions) for r in tup)
    sz = size(nuts)
    region = zeros(Int16, sz)
    updateprogress = Progress(prod(sz), 1)
    for (i, g_uid) in enumerate(nuts)
        next!(updateprogress)
        g_uid == 0 && continue
        reg = subregionnames[g_uid]
        while length(reg) >= 2
            regid = get(regionlookup, reg, 0)
            if regid > 0
                region[i] = regid
                break
            end
            reg = reg[1:end-1]
        end
    end
    return region
end

function makeregions(regiondefinitionarray)
    println("Reading GADM rasters...")
    gadmfields = readdlm("gadmfields.csv", ',', header=true)[1]
    imax = maximum(gadmfields[:,1])
    subregionnames = fill("", (imax,3))
    subregionnames[gadmfields[:,1],:] = string.(gadmfields[:,2:4])
    gadm = readraster("gadm.tif")

    regionnames = regiondefinitionarray[:,1]
    regiondefinitions = [isa(regdef, Tuple) ? regdef : (regdef,) for regdef in regiondefinitionarray[:,2]]

    println("Making region index matrix...")
    return makeregions_main(gadm, subregionnames, regiondefinitions)
end

function makeregions_main(gadm, subregionnames, regiondefinitions)
    regionlookup = build_inverseregionlookup(regiondefinitions)
    sz = size(gadm)
    region = zeros(Int16, sz)
    updateprogress = Progress(prod(sz), 1)
    for (i, g_uid) in enumerate(gadm)
        next!(updateprogress)
        g_uid == 0 && continue
        reg0, reg1, reg2 = subregionnames[g_uid,:]
        region[i] = lookup_regionnames(regionlookup, reg0, reg1, reg2)
    end
    return region
end

function lookup_regionnames(regionlookup, reg0, reg1, reg2)
    v = get(regionlookup, (reg0, "*", "*"), 0)
    v > 0 && return v
    v = get(regionlookup, (reg0, reg1, "*"), 0)
    v > 0 && return v
    return get(regionlookup, (reg0, reg1, reg2), 0)
end

function build_inverseregionlookup(regiondefinitions)
    d = Dict{Tuple{String,String,String}, Int}()
    for reg = 1:length(regiondefinitions)
        for regdef in regiondefinitions[reg]
            parentregions, subregionnames = regdef.parentregions, regdef.subregionnames
            regions = ["*", "*", "*"]
            regions[1:length(parentregions)] = parentregions
            for s in subregionnames
                regions[length(parentregions)+1] = s
                d[regions...] = reg
            end
        end
    end
    return d
end

# Find the closest region pixel for each ocean pixel and major lake.
# Even VERY far offshore pixels will be allocated to whatever region is nearest, but
# those areas still won't be available for offshore wind power because of the
# requirement to be close enough to the electricity grid (or rather the grid proxy).
function makeoffshoreregions(regions, landcover)
    println("\nMaking offshore region index matrix...")
    closest_region = regions[feature_transform(regions.>0)]
    @time lakes = majorlakes(regions)
    return closest_region .* ((regions .== 0) .| lakes) .* (landcover.==0)
end

# Use ImageSegmentation.jl to identify large lakes (large enough for offshore wind).
# This will incorrectly trigger for some very long and wide rivers, but that's not a
# big problem since near-shore areas are not allowed for offshore wind. 
function majorlakes(regions)
    println("\nIdentifying major lakes (>1000 km2)...")
    println("...segmenting regions...")
    seg = fast_scanning(regions.>0, 0.1)
    println("...removing small segments and land areas...")
    large_water_segments = prune_segments(seg, i -> (segment_pixel_count(seg,i)<1000 || segment_mean(seg,i)>0.8),
                                            (i,j) -> -segment_pixel_count(seg,j))
    lakes = labels_map(large_water_segments)
    println("...removing the largest land segment...")
    pixelcount = countmap(lakes, alg=:dict)
    # the most common index will (almost certainly) be the sole non-lake segment
    mostcommonindex = sort(collect(pixelcount), by=x->x[2], rev=true)[1][1]
    println("...lake identification complete.")    
    return lakes .!= mostcommonindex
end

struct GADM{T}
    parentregions::Vector{T}
    subregionnames::NTuple{N,T} where N
end
GADM(regionnames::T...) where T = GADM(T[], regionnames)
GADM(parentregions::Vector{T}, subregionnames::T...) where T = GADM(parentregions, subregionnames)

# struct NUTS{T}
#     subregionnames::NTuple{N,T} where N
# end
NUTS(regionnames...) = regionnames


# ArchGDAL tutorial: http://www.acgeospatial.co.uk/julia-prt3/

# ogrinfo -al -so C:/Stuff/Datasets/gadm36/gadm36.shp

# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "testtest.tif", "-a ID_0 -ts 4000 2000 -ot Byte")
# shapefile2tif("C:/Stuff/Datasets/gadm36/gadm36.shp", "Europe", "ID_0", 4300, [-11, 34, 32, 72], ')

# ogr2ogr -f CSV C:/Stuff/Julia/gadmfields012.csv -sql "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36" C:/Stuff/Datasets/gadm36/gadm36.shp
# gdal_rasterize -a UID -ot Int32 -ts 5000 2500 C:\Stuff\Datasets\gadm36\gadm36.shp C:/Stuff/Julia/globtest.tif
# gdal_rasterize -a UID -ot Int32 -ts 36000 18000 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif
# gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif

# run(`gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif`)
# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))
# rasterize("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))

# timemem-1.0 gdal_translate -r mode -tr 0.1 0.1 -co COMPRESS=LZW gadm.tif gadmsmall.tif


function rasterize_protected()
    println("Rasterizing global shapefile (10+ minute run time)...")
    shapefile = "C:/Stuff/Datasets/WDPA - Protected areas/WDPA_Feb2019-shapefile-polygons.shp"
    sql = "select FID from \"WDPA_Feb2019-shapefile-polygons\""
    outfile = "protected_raster.tif"
    options = "-a FID -a_nodata -1 -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time run(`gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile`)

    println("Creating .csv file for WDPA index and name lookup...")
    sql = "select FID,IUCN_CAT from \"WDPA_Feb2019-shapefile-polygons\""
    outfile = "protectedfields.csv"
    @time run(`ogr2ogr -f CSV $outfile -sql $sql $shapefile`)
end

function makeprotected()
    println("Reading protected area rasters...")
    protectedfields = readdlm("protectedfields.csv", ',', header=true)[1]
    IUCNcodes = ["Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Reported", "Not Applicable", "Not Assigned"]
    IUCNlookup = Dict(c => i for (i,c) in enumerate(IUCNcodes))
    protected0 = readraster("protected_raster.tif")

    println("Converting indexes to protected area types...")
    protected = similar(protected0, UInt8)
    # maybe replace loop with:  map!(p -> p == -1 ? 0 : IUCNlookup[protectedfields[p+1,2], protected, protected0)
    # alternatively             map!(p -> ifelse(p == -1, 0, IUCNlookup[protectedfields[p+1,2]), protected, protected0)
    for (i, p) in enumerate(protected0)
        protected[i] = p == -1 ? 0 : IUCNlookup[protectedfields[p+1,2]]
    end

    println("Saving protected area dataset...")
    JLD.save("protected.jld", "protected", protected, compress=true)
end


function resample(infile::String, outfile::String, options::Vector{<:AbstractString})
    @time run(`gdal_translate $options -co COMPRESS=LZW $infile $outfile`)
end

function downscale_landcover()
    println("Downscaling landcover dataset...")
    infile = "C:/Stuff/GET.GIS/datasets/landcover/Landcover - USGS MODIS.tif"
    options = "-r mode -ot Byte -tr 0.01 0.01"
    resample(infile, "landcover.tif", split(options, ' '))
end

function savelandcover()
    println("Reading landcover dataset (TIFF)...")
    landcover = readraster("landcover.tif")
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
    JLD.save("landcover.jld", "landcover", landcover, "landtypes", landtypes, "landcolors", landcolors, compress=true)
end

function upscale_topography()
    println("Upscaling topography dataset...")
    infile = "C:/Stuff/GET.GIS/datasets/topography/ETOPO1_Ice_c_geotiff.tif"
    options = "-r cubicspline -tr 0.01 0.01"
    resample(infile, "topography.tif", split(options, ' '))
    topography = readraster("topography.tif")
    println("Saving topography dataset...")
    JLD.save("topography.jld", "topography", topography, compress=true)
    rm("topography.tif")
end

# gettopography() = readraster("topography.tif")

function downscale_population(scen, year)
    scen = lowercase(scen)
    println("Reading population dataset...")
    filename = "D:/datasets/population/Gao SSP 1km/$(scen)_total_$year.nc4"
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
    JLD.save("population_$(scen)_$year.jld", "population", newpop, compress=true)

    rm(temptiff)
    rm(temptiff2)
end

getpopulation(scen, year) = JLD.load("population_$(scen)_$year.jld", "population")

function createGDP(scen, year)
    scen = lowercase(scen)
    println("Reading low resolution population and GDP datasets...")
    pop, extent = readraster("C:/Stuff/Datasets/Population & GDP/Murakami & Yamagata/pop_$(scen)_$year.tif", :getextent) # million people
    gdp = readraster("C:/Stuff/Datasets/Population & GDP/Murakami & Yamagata/gdp_$(scen)_$year.tif")    # billion USD(2005), PPP

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
    rm(tempfile)
    println("Saving high resolution GDP...")
    JLD.save("gdp_$(scen)_$year.jld", "gdp", gdphigh, compress=true)
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
    gdp = JLD.load("gdp_$(scen)_$year.jld", "gdp")
    res = 360/size(gdp,1)

    disk = diskfilterkernel(1/6/res)                        # filter radius = 1/6 degrees
    gridaccess = Float32.(imfilter(gdp .> 100_000, disk))   # only "high" income cells included (100 kUSD/cell), cell size = 1x1 km          
    selfmap!(x -> ifelse(x<1e-6, 0, x), gridaccess)         # force small values to zero to reduce dataset size
    println("Saving high resolution grid access...")
    JLD.save("gridaccess_$(scen)_$year.jld", "gridaccess", gridaccess, compress=true)

    # better:
    # loop through countries, index all pixels into vector, sort by GDP, use electrification to assign grid access
end

function getwindatlas()
    # filename = "D:/datasets/Global Wind Atlas v2.3/global_ws.tif"     # v2.3 (lon extent [-180.3, 180.3], strangely)
    filename = "C:/Stuff/GET.GIS/datasets/Global Wind Atlas/Global Wind Atlas - 100m mean wind speed.tif"   # v1.0
    windatlas = readraster(filename, :extend_to_full_globe)[1]
    clamp!(windatlas, 0, 25)
end


# grid12 + electrification
# turbine curves
# SSP
# makehydro
# listregions
