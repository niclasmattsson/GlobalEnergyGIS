import GDAL
using ArchGDAL, GDAL_jll

export rasterize, readraster, saveTIFF

function rasterize_AG(infile::String, outfile::String, options::Vector{<:AbstractString})
    ArchGDAL.read(infile) do dataset
        GDAL.close(GDAL.rasterize(
            outfile,
            Ptr{GDAL.GDALDatasetH}(C_NULL),
            dataset.ptr,
            GDAL.rasterizeoptionsnew(options, C_NULL), C_NULL))
    end
end

# works (creates the TIFF and saves it) but then crashes
function rasterize_AG2(infile::String, outfile::String, options::Vector{<:AbstractString})
    ArchGDAL.read(infile) do dataset
        ArchGDAL.unsafe_gdalrasterize(dataset, options, dest=outfile)
    end
end

# uses the command line version instead (gdal_rasterize)
# significantly faster for some reason, also gives a simple progress indication
function rasterize(infile::String, outfile::String, options::Vector{<:AbstractString}; sql::String="")
    gdal_rasterize_path() do gdal_rasterize
        if isempty(sql)
            run(`$gdal_rasterize $options $infile $outfile`)
        else
            run(`$gdal_rasterize $options -sql $sql $infile $outfile`)
        end
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
    ArchGDAL.read(infile) do dataset
        ArchGDAL.read(dataset)
    end
end

function readraster(infile::String, extentflag::Symbol, dim::Int=1)
    local raster, geotransform
    raster = ArchGDAL.read(infile) do dataset
        # display(ArchGDAL.getproj(dataset))
        geotransform = ArchGDAL.getgeotransform(dataset)
        ArchGDAL.read(dataset)[:,:,dim]
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
    width, height = size(x)
    xres = (extent[3]-extent[1])/width
    yres = (extent[4]-extent[2])/height
    raster = ArchGDAL.create(
        filename,
        driver = ArchGDAL.getdriver("GTiff"),
        width = width,
        height = height,
        nbands = 1,
        dtype = eltype(x),
        options = ["COMPRESS=LZW"]
    )
    ## assign the projection and transformation parameters
    ArchGDAL.setgeotransform!(raster, [extent[1], xres, 0, extent[4], 0, -yres])
    ArchGDAL.setproj!(raster, wkt_string)
    
    ## write the raster    
    ArchGDAL.write!(
        raster,
        x,      # image to "burn" into the raster
        1,      # update band 1
    )
    ArchGDAL.destroy(raster)
    nothing
end

saveTIFF(x::AbstractMatrix, filename::String) = saveTIFF(x, filename, [-180.0, -90.0, 180.0, 90.0])

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
