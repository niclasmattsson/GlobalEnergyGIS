import GDAL
using ArchGDAL
const AG = ArchGDAL

export rasterize, readraster, saveTIFF, GADM, makeregions, eurasia38, makeoffshoreregions, makeprotected, makelandcover
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
function rasterize(infile::String, outfile::String, options::Vector{<:AbstractString})
    GDAL_BINPATH = joinpath(dirname(pathof(GDAL)), "../deps/usr/bin")
    run(` gdal_rasterize $options $infile $outfile` )
end

function readraster(infile::String)
    ArchGDAL.registerdrivers() do
        ArchGDAL.read(infile) do dataset
            # display(ArchGDAL.getproj(dataset))
            # display(ArchGDAL.getgeotransform(dataset))
            # display(ArchGDAL.importEPSG(4326))
            dropdims(ArchGDAL.read(dataset), dims=3)
        end
    end
end

function saveTIFF(x::AbstractArray, filename::String)
    wkt_string = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]]"
    ArchGDAL.registerdrivers() do
        width, height = size(x)
        res = 360/width
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
        AG.setgeotransform!(raster, [-180, res, 0, 90, 0, -res])
        AG.setproj!(raster, wkt_string)
        
        ## write the raster    
        AG.write!(
            raster,
            x,      # image to "burn" into the raster
            1,      # update band 1
        )
        AG.destroy(raster)
    end
end

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
    @time run(` ogr2ogr -f CSV $outfile -sql $sql $shapefile` )
end

function makeregions(regiondefinitionarray)
    println("Reading GADM rasters...")
    gadmfields = readdlm("gadmfields.csv", ',', header=true)[1]
    imax = maximum(gadmfields[:,1])
    subregionnames = fill("", (imax,3))
    subregionnames[gadmfields[:,1],:] = string.(gadmfields[:,2:4])
    gadm = readraster("gadm.tif")

    regionnames = regiondefinitionarray[:,1]
    regiondefinitions = [isa(regdef, GADM) ? (regdef,) : regdef for regdef in regiondefinitionarray[:,2]]

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

# Integer version (commented out) is somewhat faster but less clear
function lookup_regionnames(regionlookup, reg0, reg1, reg2)
    v = get(regionlookup, (reg0, "*", "*"), 0)
    # v = get(regionlookup, (reg0, 0, 0), 0)
    v > 0 && return v
    v = get(regionlookup, (reg0, reg1, "*"), 0)
    # v = get(regionlookup, (reg0, reg1, 0), 0)
    v > 0 && return v
    return get(regionlookup, (reg0, reg1, reg2), 0)
end

# Integer version (commented out) is somewhat faster but less clear
function build_inverseregionlookup(regiondefinitions)
    d = Dict{Tuple{String,String,String}, Int}()
    # d = Dict{Tuple{Int,Int,Int}, Int}()
    for reg = 1:length(regiondefinitions)
        for regdef in regiondefinitions[reg]
            parentregions, subregionnames = regdef.parentregions, regdef.subregionnames
            regions = ["*", "*", "*"]
            # regions = zeros(3)
            regions[1:length(parentregions)] = parentregions
            for s in subregionnames
                regions[length(parentregions)+1] = s
                d[regions...] = reg
            end
        end
    end
    return d
end

function makeoffshoreregions(regions)
    println("Making offshore region index matrix...")
    regions[feature_transform(regions.>0)]
end

struct GADM{T}
    parentregions::Vector{T}
    subregionnames::NTuple{N,T} where N
end
GADM(regionnames::T...) where T = GADM(T[], regionnames)
GADM(parentregions::Vector{T}, subregionnames::T...) where T = GADM(parentregions, subregionnames)

eurasia38 = [
    # Europe: "NOR","FRA","GER","UK","MED","BAL","SPA","CEN"
    "NOR"   GADM("Sweden","Norway","Denmark","Finland","Åland","Faroe Islands")
    "FRA"   GADM("France","Monaco")
    "GER"   GADM("Germany","Netherlands","Belgium","Luxembourg")
    "UK"    GADM("United Kingdom","Ireland","Guernsey","Isle of Man","Jersey")
    "MED"   GADM("Greece","Bulgaria","Italy","San Marino","Vatican City","Slovenia","Croatia","Bosnia and Herzegovina","Serbia","Montenegro","Kosovo","Albania","Macedonia","Malta")
    "BAL"   GADM("Poland","Estonia","Latvia","Lithuania")
    "SPA"   GADM("Spain","Portugal","Andorra","Gibraltar")
    "CEN"   GADM("Austria","Switzerland","Czech Republic","Hungary","Slovakia","Romania","Liechtenstein")

    # Eastern Europe & Middle East: "Belarus/Ukraine/Moldova","Turkey/Caucasus","Middle East","Iran","Arabian peninsula"
    "BUK"   GADM("Belarus","Ukraine","Moldova")
    "TCC"   GADM("Turkey","Georgia","Armenia","Azerbaijan")
    "MEA"   GADM("Israel","Palestina","Lebanon","Jordan","Syria","Iraq","Kuwait")
    "IRN"   GADM("Iran")
    "ARB"   GADM("Saudi Arabia","Yemen","Oman","United Arab Emirates","Bahrain","Qatar")
    
    # Central & South-Central Asia: "Central Asia","South-Central Asia"
    "KZK"   GADM("Kazakhstan")
    "CAS"   GADM("Uzbekistan","Turkmenistan","Tajikistan","Kyrgyzstan")
    "SCA"   GADM("Afghanistan","Pakistan")
    
    # India: "North","West","Central","South","East","Northeast"
    # https://en.wikipedia.org/wiki/Administrative_divisions_of_India
    # Excluding island groups: "Andaman and Nicobar Islands","Lakshadweep"
    # South also includes Sri Lanka
    # Northeast also includes Nepal,Bhutan and Bangladesh
    "IN_N"  GADM(["India"], "Jammu and Kashmir","Himachal Pradesh","Punjab","Rajasthan","Chandigarh","Haryana","NCT of Delhi","Uttarakhand","Uttar Pradesh")
    "IN_W"  GADM(["India"], "Gujarat","Goa","Maharashtra","Dadra and Nagar Haveli","Daman and Diu")
    "IN_C"  GADM(["India"], "Madhya Pradesh","Chhattisgarh")
    "IN_S"  (GADM(["India"], "Karnataka","Kerala","Tamil Nadu","Andhra Pradesh","Telangana","Puducherry"), GADM("Sri Lanka"))
    "IN_E"  GADM(["India"], "Odisha","Jharkhand","West Bengal","Bihar")
    "IN_NE" (GADM(["India"], "Sikkim","Assam","Meghalaya","Tripura","Mizoram","Manipur","Nagaland","Arunachal Pradesh"),
                GADM("Nepal","Bhutan","Bangladesh"))

    # Southeast Asia: "Central Asia","South-Central Asia"
    # also includes Peninsular Malaysia,https://en.wikipedia.org/wiki/Peninsular_Malaysia
    "SEA"   (GADM("Myanmar","Thailand","Laos","Vietnam","Cambodia","Singapore"),
                GADM(["Malaysia"], "Perlis","Kedah","Pulau Pinang","Perak","Kelantan","Trengganu","Pahang","Selangor","Kuala Lumpur","Putrajaya","Negeri Sembilan","Melaka","Johor"))
    
    # Russia: "Northwest","Central","Southwest","Volga","Ural","Siberia","East"
    # https://en.wikipedia.org/wiki/Federal_districts_of_Russia
    # questionable: Novgorod,Altay,Yevrey,Maga Buryatdan/Magadan
    "RU_NW"  GADM(["Russia"], "Arkhangel'sk","Vologda","Kaliningrad","Karelia","Komi","Leningrad","Murmansk","Nenets","Novgorod","Pskov","City of St. Petersburg")
    "RU_C"   GADM(["Russia"], "Belgorod","Bryansk","Vladimir","Voronezh","Ivanovo","Kaluga","Kostroma","Kursk","Lipetsk","Moscow City","Moskva","Orel","Ryazan'","Smolensk","Tambov","Tver'","Tula","Yaroslavl'")
    "RU_SW"  GADM(["Russia"], "Adygey","Astrakhan'","Volgograd","Kalmyk","Krasnodar","Rostov","Dagestan","Ingush","Kabardin-Balkar","Karachay-Cherkess","North Ossetia","Stavropol'","Chechnya")
    "RU_VL"  GADM(["Russia"], "Bashkortostan","Kirov","Mariy-El","Mordovia","Nizhegorod","Orenburg","Penza","Perm'","Samara","Saratov","Tatarstan","Udmurt","Ul'yanovsk","Chuvash")
    "RU_UR"  GADM(["Russia"], "Kurgan","Sverdlovsk","Tyumen'","Khanty-Mansiy","Chelyabinsk","Yamal-Nenets")
    "RU_SB"  GADM(["Russia"], "Altay","Gorno-Altay","Irkutsk","Kemerovo","Krasnoyarsk","Novosibirsk","Omsk","Tomsk","Tuva","Khakass")
    "RU_E"   GADM(["Russia"], "Amur","Buryat","Yevrey","Zabaykal'ye","Kamchatka","Maga Buryatdan","Primor'ye","Sakha","Sakhalin","Khabarovsk","Chukot")
    
    # China: "North" (Huáběi),"Northeast" (Dōngběi),"East" (Huádōng),"South Central" (Zhōngnán),"Southwest" (Xīnán),"Northwest" (Xīběi)
    # https://en.wikipedia.org/wiki/List_of_regions_of_China
    # East also includes Taiwan
    # South Central also includes Hong Kong and Macao
    "CH_N"   GADM(["China"], "Beijing","Tianjin","Hebei","Shanxi","Nei Mongol")
    "CH_NE"  GADM(["China"], "Liaoning","Jilin","Heilongjiang")
    "CH_E"   (GADM(["China"], "Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi","Shandong"), GADM("Taiwan"))
    "CH_SC"  (GADM(["China"], "Henan","Hubei","Hunan","Guangdong","Guangxi","Hainan"), GADM("Hong Kong","Macao"))
    "CH_SW"  GADM(["China"], "Chongqing","Sichuan","Guizhou","Yunnan","Xizang")
    "CH_NW"  GADM(["China"], "Shaanxi","Gansu","Qinghai","Ningxia Hui","Xinjiang Uygur")

    # Other
    "MON"    GADM("Mongolia")
    "JKR"    GADM("Japan","South Korea","North Korea")
]



# ArchGDAL tutorial: http://www.acgeospatial.co.uk/julia-prt3/

# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "testtest.tif", "-a ID_0 -ts 4000 2000 -ot Byte")
# shapefile2tif("C:/Stuff/Datasets/gadm36/gadm36.shp", "Europe", "ID_0", 4300, [-11, 34, 32, 72], ')

# ogr2ogr -f CSV C:/Stuff/Julia/gadmfields012.csv -sql "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36" C:/Stuff/Datasets/gadm36/gadm36.shp
# gdal_rasterize -a UID -ot Int32 -ts 5000 2500 C:\Stuff\Datasets\gadm36\gadm36.shp C:/Stuff/Julia/globtest.tif
# gdal_rasterize -a UID -ot Int32 -ts 36000 18000 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif
# gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif

# run(` gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif` )
# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))
# rasterize("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))

# timemem-1.0 gdal_translate -r mode -tr 0.1 0.1 -co COMPRESS=LZW gadm.tif gadmsmall.tif


function rasterize_protected()
    println("Rasterizing global shapefile (10+ minute run time)...")
    shapefile = "C:/Stuff/Datasets/WDPA - Protected areas/WDPA_Feb2019-shapefile-polygons.shp"
    sql = "select FID from \"WDPA_Feb2019-shapefile-polygons\""
    outfile = "protected_raster.tif"
    options = "-a FID -a_nodata -1 -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time run(` gdal_rasterize $(split(options, ' ')) -sql $sql $shapefile $outfile` )

    println("Creating .csv file for WDPA index and name lookup...")
    sql = "select FID,IUCN_CAT from \"WDPA_Feb2019-shapefile-polygons\""
    outfile = "protectedfields.csv"
    @time run(` ogr2ogr -f CSV $outfile -sql $sql $shapefile` )
end

function makeprotected()
    println("Reading protected area rasters...")
    protectedfields = readdlm("protectedfields.csv", ',', header=true)[1]
    IUCNcodes = ["Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Reported", "Not Applicable", "Not Assigned"]
    IUCNlookup = Dict(c => i for (i,c) in enumerate(IUCNcodes))
    protected0 = readraster("protected_raster.tif")

    println("Converting indexes to protected area types...")
    protected = similar(protected0, UInt8)
    for (i, p) in enumerate(protected0)
        protected[i] = p == -1 ? 0 : IUCNlookup[protectedfields[p+1,2]]
    end

    println("Saving protected area dataset...")
    saveTIFF(protected, "protected.tif")
end


function resample(infile::String, outfile::String, options::Vector{<:AbstractString})
    @time run(` gdal_translate $options -co COMPRESS=LZW $infile $outfile` )
end

function downscale_landcover()
    println("Downscaling landcover dataset...")
    infile = "C:/Stuff/GET.GIS/datasets/landcover/Landcover - USGS MODIS.tif"
    options = "-r mode -ot Byte -tr 0.01 0.01"
    resample(infile, "landcover.tif", split(options, ' '))
end

function getlandcover()
    landcover = readraster("landcover.tif")
    landtypes = [
        "Barren", "Snow/Ice", "Cropland/Natural", "Urban", "Croplands",
        "Permanent Wetlands", "Grasslands", "Savannas", "Woody Savannas", "Open Shrublands",
        "Closed Shrublands", "Mixed Forests", "Deciduous Broadleaf Forests", "Deciduous Needleleaf Forests", "Evergreen Broadleaf Forests",
        "Evergreen Needleleaf Forests", "Water"
    ]
    # landcolors = 1/255 * [
    #     190 190 190; 255 218 209; 144 144 0; 255 0 0; 255 255 0;
    #     40 136 213; 255 192 107; 255 228 18; 182 231 140; 255 236 163;
    #     216 118 118; 55 200 133; 104 229 104; 123 204 6; 77 167 86;
    #     0 100 0; 190 247 255
    # ]
    return landcover, landtypes
end

function upscale_topography()
    println("Upscaling topography dataset...")
    infile = "C:/Stuff/GET.GIS/datasets/topography/ETOPO1_Ice_c_geotiff.tif"
    options = "-r cubicspline -tr 0.01 0.01"
    resample(infile, "topography.tif", split(options, ' '))
end

gettopography() = readraster("topography.tif")

function downscale_population(year)
    println("Reading population dataset...")
    filename = "D:/datasets/population/Gao SSP2_1km/ssp2_total_$year.nc4"
    pop = ncread(filename, "Band1")
    lat = ncread(filename, "lat")
    res = 0.5/60    # source resolution 0.5 arcminutes

    println("Padding and saving intermediate dataset...")
    skiptop = round(Int, (90-(lat[end]+res/2)) / res)
    skipbottom = round(Int, (lat[1]-res/2-(-90)) / res)
    pop = pop'
    pop[pop.<0] .= 0
    lons = size(pop,2)
    # the factor (.01/res)^2 is needed to conserve total population
    pop = collect([zeros(Float32,skiptop,lons); reverse(pop, dims=1)*(.01/res)^2; zeros(Float32,skipbottom,lons)]')
    temptiff = tempname()
    saveTIFF(pop, temptiff)

    println("Downscaling population dataset...")
    options = "-r cubicspline -tr 0.01 0.01"
    resample(temptiff, "population.tif", split(options, ' '))
    rm(temptiff)
end

getpopulation() = readraster("population.tif")

# grid12 + electrification
# windatlas
# turbine curves
# SSP
# GDP (maybe)
# makehydro

