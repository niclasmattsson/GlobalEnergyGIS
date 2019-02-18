import GDAL
using ArchGDAL; const AG = ArchGDAL

export shapefile2tif

function agraster2tif(infile::String, outfile::String, options::Vector{<:AbstractString})
    AG.registerdrivers() do
        AG.read(infile) do source_ds
            GDAL.close(GDAL.rasterize(
                outfile,
                Ptr{GDAL.GDALDatasetH}(C_NULL),
                source_ds.ptr,
                GDAL.rasterizeoptionsnew(options, C_NULL), C_NULL))
        end
    end
end

function shapefile2tif(shapefile::String, outfile::String, width::Int, bbox::Vector{<:Real})
	GDAL_BINPATH = "$(DEPOT_PATH[1])\\packages\\GDAL\\vec6Y\\deps\\usr\\bin"
	bboxstr = join(bbox, " ")
	height = round(Int, width*(bbox[4]-bbox[2])/(bbox[3]-bbox[1]))
	tempfile = "temp_shapefiles/$outfile.shp"
	run(`$GDAL_BINPATH\\ogr2ogr -clipsrc $bboxstr $tempfile $shapefile`)
	agraster2tif("temp_shapefiles/$outfile.shp", "$outfile.tif", split("-a ID_0 -ts $height $width -ot Byte"))
end

gadmlevel = 1
shapefile2tif("C:\\Stuff\\GET.GIS\\datasets\\gadm36\\gadm36_$gadmlevel.shp", "Europe", 4300, [-11, 34, 32, 72])
