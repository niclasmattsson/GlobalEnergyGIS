import GDAL
using ArchGDAL
const AG = ArchGDAL

export rasterize_CLI, readraster, globalGADMtiff, GADM, makeregions, eurasia38, makeoffshoreregions
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

# uses the command line version instead (gdal_rasterize)
# significantly faster for some reason, also give a simple progress indication
function rasterize_CLI(infile::String, outfile::String, options::Vector{<:AbstractString})
    GDAL_BINPATH = joinpath(dirname(pathof(GDAL)), "../deps/usr/bin")
    run(` gdal_rasterize $options $infile $outfile` )
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

function readraster(infile::String)
    ArchGDAL.registerdrivers() do
        ArchGDAL.read(infile) do dataset
            dropdims(ArchGDAL.read(dataset), dims=3)
        end
    end
end

function shapefile2tif(shapefile::String, outfile::String, layer::String, width::Int, bbox::Vector{<:Real}, options::String)
	GDAL_BINPATH = joinpath(dirname(pathof(GDAL)), "../deps/usr/bin")
	bboxstr = join(bbox, " ")
	height = round(Int, width*(bbox[4]-bbox[2])/(bbox[3]-bbox[1]))
    path = mkpath(joinpath(tempdir(), "temp_shapefiles"))
    cmd = `$GDAL_BINPATH/ogr2ogr -lco ENCODING=UTF-8 -clipsrc $bbox $path/$outfile.shp $shapefile`
    # println("$bboxstr $tempfile $shapefile")
    # dump(cmd)
    println("Clipping shapefile to bounding box...")
	run(cmd)
    println("Rasterizing clipped shapefile...")
	rasterize_AG("$path/$outfile.shp", "$outfile.tif", split("-a $layer -ts $height $width -ot Byte"))
end

function globalGADMtiff()
    println("Rasterizing global shapefile...")
    shapefile = "C:/Stuff/Datasets/gadm36/gadm36.shp"
    outfile = "gadm.tif"
    options = "-a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"
    # options = "-a UID -ot Int32 -tr 0.02 0.02 -te -180 -90 180 90 -co COMPRESS=LZW"
    @time rasterize_CLI(shapefile, outfile, split(options, ' '))
 
    println("Creating .csv file for regional index and name lookup...")
    sql = "select uid,name_0,name_1,name_2 from gadm36"
    # sql = "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36"
    outfile = "gadmfields.csv"
    @time run(` ogr2ogr -f CSV $outfile -sql $sql $shapefile` )
end

function read_gadmfields()
    gadmfields = readdlm("gadmfields.csv", ',', header=true)[1]
    gadm = readraster("gadm.tif")
    # getmax(x) = maximum([isempty(e) ? typemin(Int64) : e for e in x])
    # ncountries = getmax(g[:,2])
    gadmregions = Dict(0 => unique(g[:,3]), 1 => unique(g[:,5]), 2 => unique(g[:,7]), 3 => g[:,1])
    gadmregionindex = Dict(level => Dict(gadmregions[level][i] => i for i=1:length(gadmregions[level])) for level=0:3)
    index = Dict{Int,Int64}()
    id_multiplier = Dict(0 => 10^(6+5+4), 1 => 10^(6+5), 2 => 10^6, 3 => 1)
    for row in eachrow(gadmfields)
        index::Int64 = 0
        for level = 0:3
            column = level == 3 ? 1 : 3 + level*2
            index += id_multiplier[level] * gadmregionindex[level][row[column]]
        end

    end
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
    return regiondictloop(gadm, subregionnames, regiondefinitions)
end

makeoffshoreregions(regions) = regions[feature_transform(region.>0)]

function regionloop(gadm, subregionnames, regiondefinitions)
    sz = size(gadm)
    region = zeros(Int16, sz)
    updateprogress = Progress(prod(sz), 1)
    for (i, g_uid) in enumerate(gadm)
        next!(updateprogress)
        g_uid == 0 && continue
        reg0, reg1, reg2 = subregionnames[g_uid,:]
        region[i] = regionlookup(regiondefinitions, reg0, reg1, reg2)
    end
    return region
end

function regiontreeloop(gadm, subregionnames, regiondefinitions)
    regiontree = buildregiontree(regiondefinitions)
    sz = size(gadm)
    region = zeros(Int16, sz)
    updateprogress = Progress(prod(sz), 1)
    for (i, g_uid) in enumerate(gadm)
        next!(updateprogress)
        g_uid == 0 && continue
        reg0, reg1, reg2 = subregionnames[g_uid,:]
        region[i] = regiontreelookup(regiontree, reg0, reg1, reg2)
    end
    return region
end

function regiondictloop(gadm, subregionnames, regiondefinitions)
    regiondict = buildregiondict(regiondefinitions)
    sz = size(gadm)
    region = zeros(Int16, sz)
    updateprogress = Progress(prod(sz), 1)
    for (i, g_uid) in enumerate(gadm)
        next!(updateprogress)
        g_uid == 0 && continue
        reg0, reg1, reg2 = subregionnames[g_uid,:]
        region[i] = regiondictlookup(regiondict, reg0, reg1, reg2)
    end
    return region
end

# function regiondictloop_threads(gadm, subregionnames, regiondefinitions)
#     regiondict = buildregiondict(regiondefinitions)
#     sz = size(gadm)
#     region = [zeros(Int16, sz) for i = 1:Threads.nthreads()]
#     updateprogress = Progress(prod(sz), 1)
#     Threads.@threads for i = 1:prod(sz)
#         next!(updateprogress)
#         g_uid = gadm[i]
#         g_uid == 0 && continue
#         reg0, reg1, reg2 = subregionnames[g_uid,:]
#         region[Threads.threadid()][i] = regiondictlookup(regiondefinitions, reg0, reg1, reg2)
#     end
#     return maximum(region, dims=3)
# end

# function regionloop2(gadm, subregionnames, regiondefinitions)
#     sz = size(gadm)
#     region = [zeros(Int16, sz) for i = 1:Threads.nthreads()]
#     updateprogress = Progress(prod(sz), 1)
#     Threads.@threads for i = 1:prod(sz)
#         next!(updateprogress)
#         g_uid = gadm[i]
#         g_uid == 0 && continue
#         reg0, reg1, reg2 = subregionnames[g_uid,:]
#         region[Threads.threadid()][i] = regionlookup(regiondefinitions, reg0, reg1, reg2)
#     end
#     return maximum(region, dims=3)
# end

# function testtest()
#     lookup = rand(10)
#     rmat = rand(1:10, (5,6))
#     Threads.@threads for r in rmat
#         # r == 0 && continue
#         val = lookup[r]
#     end
#     return lookup, rmat
# end

function makeregions_faster(regiondefinitionarray)
    println("Reading GADM rasters...")
    gadmfields = readdlm("gadmfields.csv", ',', header=true)[1]
    imax = maximum(gadmfields[:,1])
    subregionnames = Dict(level => unique(gadmfields[:, level+2]) for level = 0:2)
    subregionindexes = Dict(level => Dict(subregionnames[level][i] => i for i=1:length(subregionnames[level])) for level = 0:2)

    subregioncodes = zeros(Int,imax,3)
    for i=1:imax, j=1:3
        subregioncodes[gadmfields[i,1],j] = subregionindexes[j-1][gadmfields[i,j+1]]
    end
    gadm = readraster("gadm.tif")

    regionnames = regiondefinitionarray[:,1]
    regiondefinitions = [isa(regdef, GADM) ? (regdef,) : regdef for regdef in regiondefinitionarray[:,2]]
    regiondefinitioncodes = [Tuple(
        GADM([subregionindexes[length(regdef.parentregions)-1][r] for r in regdef.parentregions],
            Tuple(subregionindexes[length(regdef.parentregions)][t] for t in regdef.subregionnames))
                                for regdef in regiondefinitions[reg]) for reg = 1:length(regiondefinitions)]

    println("Making region index matrix...")
    return regiondictloop(gadm, subregioncodes, regiondefinitioncodes)
end

# regiontreelookup(rootnode, reg0, reg1, reg2)
function regiontreelookup(node, reglist...)
    isempty(reglist) && return node.regionindex
    child = findfirst(x -> x.name == reglist[1], node.children)
    child == nothing ? node.regionindex : regiontreelookup(node.children[child], reglist[2:end]...)
end

function regionlookup(regiondefinitions, reg0, reg1, reg2)
    for reg = 1:length(regiondefinitions)
        for regdef in regiondefinitions[reg]
            parentregions, subregionnames = regdef.parentregions, regdef.subregionnames
            level = length(parentregions)
            level == 0 && reg0 in subregionnames && return reg
            level == 1 && reg0 == parentregions[1] && reg1 in subregionnames && return reg
            level == 2 && reg0 == parentregions[1] && reg1 == parentregions[2] && reg2 in subregionnames && return reg
        end
    end
    return 0
end

function regiondictlookup(regiondict, reg0, reg1, reg2)
    v = get(regiondict, (reg0, "*", "*"), 0)
    # v = get(regiondict, (reg0, 0, 0), 0)
    v > 0 && return v
    v = get(regiondict, (reg0, reg1, "*"), 0)
    # v = get(regiondict, (reg0, reg1, 0), 0)
    v > 0 && return v
    return get(regiondict, (reg0, reg1, reg2), 0)
end

# struct Node{T}
#     name::T
#     regionindex::Int64
#     children::Vector{Node{T}}
# end
# Node(name::T, regionindex::Int64) where T = Node(name, regionindex, Node{T}[])

# import AbstractTrees: children, printnode, print_tree, Leaves
# children(n::Node) = n.children
# inchildren(s, n::Node) = any(c -> (c.name == s), n.children)
# function findchild(s, n::Node)
#     index = findfirst(c -> (c.name == s), n.children)
#     if index == nothing
#         return nothing
#     else
#         return n.children[index]
#     end
# end

# function insertparentregions(tree, parentregions)
#     isempty(parentregions) && return tree
#     parent = findchild(parentregions[1], tree)
#     if parent == nothing
#         parent = Node(parentregions[1], 0)
#         push!(tree.children, parent)
#     end
#     insertparentregions(parent, parentregions[2:end])
# end

# function buildregiontree(regiondefinitions)
#     lookuptype = eltype(regiondefinitions[1][1].parentregions)
#     if lookuptype == String
#         tree = Node("root", 0)
#     elseif lookuptype == Int
#         tree = Node(0, 0)
#     else
#         error("Unexpected type $lookuptype")
#     end
#     for reg = 1:length(regiondefinitions)
#         for regdef in regiondefinitions[reg]
#             parentregions, subregionnames = regdef.parentregions, regdef.subregionnames
#             parent = insertparentregions(tree, parentregions)
#             for s in subregionnames
#                 push!(parent.children, Node(s, reg))
#             end
#         end
#     end
#     return tree
# end

function buildregiondict(regiondefinitions)
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

# const gadmfields = read_gadmfields()

# function regionlookup(regiondict)
#     regionnames = keys(regiondict)
#     inversedict = Dict()
#     lookupfunction = (r0,r1,r2) -> lookup(inversedict, r0, r1, r2)
#     return regionnames, lookupfunction
# end

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



# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "testtest.tif", "-a ID_0 -ts 4000 2000 -ot Byte")
# shapefile2tif("C:/Stuff/Datasets/gadm36/gadm36.shp", "Europe", "ID_0", 4300, [-11, 34, 32, 72], ')

# ogr2ogr -f CSV C:/Stuff/Julia/gadmfields012.csv -sql "select uid,id_0,name_0,id_1,name_1,id_2,name_2 from gadm36" C:/Stuff/Datasets/gadm36/gadm36.shp
# gdal_rasterize -a UID -ot Int32 -ts 5000 2500 C:\Stuff\Datasets\gadm36\gadm36.shp C:/Stuff/Julia/globtest.tif
# gdal_rasterize -a UID -ot Int32 -ts 36000 18000 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif
# gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif

# run(` gdal_rasterize -a UID -ot Int32 -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW C:/Stuff/Datasets/gadm36/gadm36.shp C:/Users/niclas/Downloads/globtest.tif` )
# rasterize_AG("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))
# rasterize_CLI("C:/Stuff/Datasets/gadm36/gadm36.shp", "C:/Users/niclas/Downloads/globtest3.tif", split("-a ID_0 -ot Byte -tr 0.01 0.01 -te -180 -90 180 90 -co COMPRESS=LZW"))

# timemem-1.0 gdal_translate -r mode -tr 0.1 0.1 -co COMPRESS=LZW gadm.tif gadmsmall.tif
