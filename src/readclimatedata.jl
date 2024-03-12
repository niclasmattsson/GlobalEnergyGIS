export savewind

include("coordinatedescent.jl")

const MODELDATA_CORDEX = Dict(
    "CNRM-CM5" =>   ("CNRM-CERFACS-CNRM-CM5", "r1i1p1", "v2"),
    "MPI-ESM" =>    ("MPI-M-MPI-ESM-LR", "r1i1p1", "v1"),
    "EC-EARTH" =>   ("ICHEC-EC-EARTH", "r12i1p1", "v1"),
    "NorESM1-M" =>  ("NCC-NorESM1-M", "r1i1p1", "v1"),
    "HadGEM2-ES" => ("MOHC-HadGEM2-ES", "r1i1p1", "v1")
)

const MODELDATA_HCLIM = Dict(
    "EC-EARTH" => ("EC-Earth_driven", "ICHEC-EC-EARTH")
)

const ALLSIMS = [
    ("CORDEX", "ictp", "CNRM-CM5", 85),
    ("CORDEX", "cnrm", "CNRM-CM5", 85),
    ("CORDEX", "cnrm", "CNRM-CM5", 26),
    ("CORDEX", "ictp", "MPI-ESM", 85),
    ("CORDEX", "ictp", "MPI-ESM", 26),
    ("CORDEX", "ictp", "EC-EARTH", 85),
    ("CORDEX", "ictp", "NorESM1-M", 85),
    ("CORDEX", "ictp", "NorESM1-M", 26),
    ("CORDEX", "cnrm", "NorESM1-M", 85),
    ("CORDEX", "cnrm", "HadGEM2-ES", 85),
    ("HCLIM", "", "EC-EARTH", 85),
]

function simname_hclim(model, rcp, variable, altitude, year)
    driven, modelname = MODELDATA_HCLIM[model]
    varalt = variable in ["u", "v"] ? "$(variable)a$(altitude)m" : variable
    yearrange = year <= 2060 ? "2040_2060" : "2080_2100"
    rcpfolder = year > 2005 ? "RCP$(rcp)_$yearrange" : "historical"
    rcpname = year > 2005 ? "rcp$(rcp)" : "historical"
    climatefolder = getconfig("climatefolder")
    path = "$climatefolder/hclim3km/$driven/$rcpfolder/3hr/$varalt"
    time = variable in ["rsns", "rsdsdir"] ? "3hr_$(year)01010130-$(year)12312230" :
                                             "3hr_$(year)01010000-$(year)12312100"
    file = "$(varalt)_NEU-3_$(modelname)_$(rcpname)_r12i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_$(time).nc"
    return "$path/$file"
end

function simname_cordex(model, org, rcp, variable, altitude, year)
    modelname, rip, v1v2 = MODELDATA_CORDEX[model]
    org = uppercase(org)
    variant = org == "ICTP" ? "RegCM4-6" : "ALADIN63"
    varalt = variable in ["u", "v"] ? "$(variable)a$(altitude)m" : variable
    rcpfolder = year > 2005 ? "rcp$rcp" : "historical"
    climatefolder = getconfig("climatefolder")
    path = "$climatefolder/CORDEX/$org/$modelname/$rcpfolder/$rip/$variant/$v1v2/3hr/$varalt/latest"
    time = variable == "rsds" ? "3hr_$(year)01010130-$(year)12312230" :
                                "3hr_$(year)01010300-$(year+1)01010000"
    file = "$(varalt)_EUR-11_$(modelname)_$(rcpfolder)_$(rip)_$(org)-$(variant)_$(v1v2)_$(time).nc"
    return "$path/$file"
end

function nearest_lonlat(lons, lats, xylon, xylat)
    indices = zeros(CartesianIndex{2}, length(lons), length(lats))
    distances = zeros(length(lons), length(lats))
    @time for (i, lon) in enumerate(lons)
        lastindex = i > 1 ? indices[i-1, 1] : CartesianIndex(200,200)
        for (j, lat) in enumerate(lats)
            f = ndx -> checkbounds(Bool, xylon, ndx) ? (xylon[ndx] - lon)^2 + (xylat[ndx] - lat)^2 : Inf
            dmin, newindex = coordinate_descent(f, lastindex)
            distances[i, j] = dmin
            indices[i, j] = newindex
        end
    end
    return distances, indices
end

extent2range(extent, res) =
    range(extent[1]+res/2, extent[3]-res/2, step=res), range(extent[4]-res/2, extent[2]+res/2, step=-res)

function filename_wind(datasource, org, model, rcp, altitude, year)
    rcpname = year > 2005 ? "rcp$(rcp)" : "historical"
    orgname = isempty(org) ? "" : "_$org"
    climatefolder = getconfig("climatefolder")
    filename = "$climatefolder/wind_$(datasource)$(orgname)_$(model)_$(altitude)m_$(rcpname)_$(year).h5"
end

function savewind(datasource, org, model, rcp, altitude, year)
    datasource, org = uppercase(datasource), lowercase(org)
    filename = filename_wind(datasource, org, model, rcp, altitude, year)
    println("\nProcessing $filename...")
    println("Reprojecting u variable...")
    u, extent = reproject(datasource, org, model, rcp, "u", altitude, year)
    println("Reprojecting v variable...")
    v, _ = reproject(datasource, org, model, rcp, "v", altitude, year)
    println("Calculating absolute and annual average wind speeds...")
    @time begin
        wind = sqrt.(u.^2 + v.^2)
        meanwind = meandrop(wind, dims=3)
    end
    sz = size(wind)
    println("Saving as $filename...")
    @time h5open(filename, "w") do file 
        group = file["/"]
        dataset_data = create_dataset(group, "wind", datatype(Float32), dataspace(sz), chunk=(16,16,sz[3]), blosc=3)
        dataset_mean = create_dataset(group, "meanwind", datatype(Float32), dataspace(sz[1:2]), chunk=(16,16), blosc=3)
        dataset_extent = create_dataset(group, "extent", datatype(Float64), dataspace(size(extent)))
        dataset_data[:,:,:] = wind
        dataset_mean[:,:] = meanwind
        dataset_extent[:] = extent
    end
    nothing
end

function savesolartemp(datavar, datasource, org, model, rcp, altitude, year)
    datasource, org = uppercase(datasource), lowercase(org)
    rcpname = year > 2005 ? "rcp$(rcp)" : "historical"
    orgname = isempty(org) ? "" : "_$org"
    climatefolder = getconfig("climatefolder")
    filename = "$climatefolder/$(datavar)_$(datasource)$(orgname)_$(model)_$(altitude)m_$(rcpname)_$(year).h5"
    println("\nProcessing $filename...")
    println("Reprojecting...")
    if datavar == "temp"
        var = datasource == "CORDEX" ? "tas" : "tasLand"
    else
        var = datasource == "CORDEX" ? "rsds" : "rsns"
    end
    data, extent = reproject(datasource, org, model, rcp, var, altitude, year)
    @time meandata = meandrop(data, dims=3)
    sz = size(data)
    println("Saving as $filename...")
    @time h5open(filename, "w") do file 
        group = file["/"]
        dataset_data = create_dataset(group, datavar, datatype(Float32), dataspace(sz), chunk=(16,16,sz[3]), blosc=3)
        dataset_mean = create_dataset(group, "mean$datavar", datatype(Float32), dataspace(sz[1:2]), chunk=(16,16), blosc=3)
        dataset_extent = create_dataset(group, "extent", datatype(Float64), dataspace(size(extent)))
        dataset_data[:,:,:] = data
        dataset_mean[:,:] = meandata
        dataset_extent[:] = extent
    end
    nothing
end

function save_climate_data()
    for (datasource, org, model, rcp) in ALLSIMS
        @time savewind(datasource, org, model, rcp, 100, 2050)
        model == "EC-EARTH" && @time savewind(datasource, org, model, rcp, 100, 2005)
        for st in ["solar", "temp"]
            @time savesolartemp(st, datasource, org, model, rcp, 0, 2050)
            model == "EC-EARTH" && datasource == "HCLIM" && @time savesolartemp(st, datasource, org, model, rcp, 0, 2005)
        end
    end
end

function savewind_10year()
    for year in [1996:2004; 2046:2049; 2051:2055; 2091:2100]
        savewind("CORDEX", "ictp", "EC-EARTH", 85, 100, year)
        savewind("HCLIM", "", "EC-EARTH", 85, 100, year)
    end
end

# create_10year_meanwind("CORDEX", "ictp", "EC-EARTH", 85, 100, 1996:2005)
# create_10year_meanwind("CORDEX", "ictp", "EC-EARTH", 85, 100, 2046:2055)
# create_10year_meanwind("CORDEX", "ictp", "EC-EARTH", 85, 100, 2091:2100)
# create_10year_meanwind("HCLIM", "", "EC-EARTH", 85, 100, 1996:2005)
# create_10year_meanwind("HCLIM", "", "EC-EARTH", 85, 100, 2046:2055)
# create_10year_meanwind("HCLIM", "", "EC-EARTH", 85, 100, 2091:2100)
function create_10year_meanwind(datasource, org, model, rcp, altitude, years)
    meanwindall = h5read(filename_wind(datasource, org, model, rcp, altitude, years[1]), "/meanwind")
    extent = h5read(filename_wind(datasource, org, model, rcp, altitude, years[1]), "/extent")
    for year in years[2:end]
        meanwindall += h5read(filename_wind(datasource, org, model, rcp, altitude, year), "/meanwind")
    end
    climatefolder = getconfig("climatefolder")
    @time h5open("$climatefolder/meanwind_$(years[1])_$(years[end])_$(datasource)_$(model)_$(altitude)m.h5", "w") do file 
        group = file["/"]
        dataset_mean = create_dataset(group, "meanwind", datatype(Float32), dataspace(size(meanwindall)), chunk=(16,16), blosc=3)
        dataset_extent = create_dataset(group, "extent", datatype(Float64), dataspace(size(extent)))
        dataset_mean[:,:] = meanwindall/length(years)
        dataset_extent[:] = extent
    end
end

function create_10year_meanwind_era5(years)
    meanwindall = h5read(in_datafolder("era5wind$(years[1]).h5"), "/meanwind")
    for year in years[2:end]
        meanwindall += h5read(in_datafolder("era5wind$(year).h5"), "/meanwind")
    end
    climatefolder = getconfig("climatefolder")
    @time h5open("$climatefolder/meanwind_$(years[1])_$(years[end])_ERA5.h5", "w") do file 
        group = file["/"]
        dataset_mean = create_dataset(group, "meanwind", datatype(Float32), dataspace(size(meanwindall)), chunk=(16,16), blosc=3)
        dataset_mean[:,:] = meanwindall/length(years)
    end
end

function annual_wind_deviations(datasource, org, model, rcp, altitude)
    years = [1996:2005; 2046:2055; 2091:2100]
    nyears = length(years)
    dev = zeros(nyears, nyears)
    gisregion = (datasource == "HCLIM") ? "NEurope8" : "CORDEXEurope8"
    regions, _, regionlist, lonrange, latrange = loadregions(gisregion);
    erares = (datasource == "ERA5") ? 0.28125 : (datasource == "HCLIM") ? 0.03 : 0.10
    nlon, nlat = length(lonrange), length(latrange)
    ilons, ilats = (datasource == "HCLIM") ? (2:nlon-1,2:nlat-1) : (1:nlon,1:nlat)
    smallregions = resize_categorical(regions[ilons,ilats], regionlist, lonrange[ilons], latrange[ilats], erares; skipNOREGION=true)
    ok = (datasource == "HCLIM") ? (smallregions .>= 1 .&& smallregions .<= 4) : (smallregions .== 1)
    for (i, year1) in enumerate(years)
        meanwind1 = h5read(filename_wind(datasource, org, model, rcp, altitude, year1), "/meanwind")[ok]
        for (j, year2) in enumerate(years)
            year1 == year2 && continue
            meanwind2 = h5read(filename_wind(datasource, org, model, rcp, altitude, year2), "/meanwind")[ok]
            dev[j,i] = mean(abs.(meanwind1 - meanwind2))
        end
    end
    return dev
end

function plot_meandata(var)
    dir = "E:/clim"
    files = readdir(dir, join=true)
    for file in files
        base, ext = splitext(file)
        if startswith(basename(base), var) && ext == ".h5"
            println(file)
            meandata = h5read(file, "/mean$var")
            cr = var == "wind" ? (2,11) : var == "solar" ? (60,250) : (263,295)
            println(extrema(meandata[isfinite.(meandata)]), " ", cr)
            s = heatmap(reverse(meandata, dims=2); colorrange=cr)
            Makie.save("$base.png", s, resolution=(800,600).*4)
        end
    end
end

remove_singleton_dims(a) = dropdims(a, dims = (findall(size(a) .== 1)...,))

function reproject(datasource, org, model, rcp, variable, altitude, year)
    if datasource == "HCLIM"
        simname = simname_hclim(model, rcp, variable, altitude, year)
        res = 0.03
        extent = [2, 50.5, 32, 71.5]
        maxdist = 2.1 * res^2
    else
        simname = simname_cordex(model, org, rcp, variable, altitude, year)
        res = 0.1
        extent = [-10.5, 34.5, 32, 71.5]
        maxdist = 2.7 * res^2
    end
    println("Loading NetCDF file...")
    @time begin
        nc = Dataset(simname)
        xylon = nc["lon"][:,:]
        xylat = nc["lat"][:,:]
        varalt = variable in ["u", "v"] ? "$(variable)a$(altitude)m" : variable
        ncdata = remove_singleton_dims(Float32.(replace(nc[varalt][:,:,:], missing => Inf32)))
        nhours = size(ncdata, 3)
        lons, lats = extent2range(extent, res)
    end
    println("Calculating nearest neighbors...")
    distances, indices = nearest_lonlat(lons, lats, xylon, xylat)
    println("Reading time series data from nearest neighbors...")
    @time begin
        newdata = zeros(Float32, size(indices)..., nhours)
        infdata = fill(Inf32, nhours)
        ci = CartesianIndices(indices)
        for (i,j) in enumerate(indices)
            index = ci[i]
            newdata[index, :] = distances[index] < maxdist ? ncdata[j, :] : infdata
        end
    end
    return newdata, extent
end

function reproject_hclim_gdal(altitude, year)
    println("Reprojecting u variable...")
    filename1 = reproject_hclim_gdal("u", altitude, year)
    println("Reprojecting v variable...")
    filename2 = reproject_hclim_gdal("v", altitude, year)
    println("Read & save absolute wind speed...")
    saveSMHIwind(altitude, year)
    rm(filename1)
    rm(filename1 * ".aux.xml")
    rm(filename2)
    rm(filename2 * ".aux.xml")
end

function reproject_hclim_gdal(variable, altitude, year)
    varalt = "$(variable)a$(altitude)m"
    bulkname = "NEU-3_ECMWF-ERAINT_evaluation_r1i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_3hr"
    ncfile = "D:/SMHI/$(varalt)_$(bulkname)_$(year)01010000-$(year)12312100.nc"
    outfile = "D:/SMHI/TEST$(varalt)_$(year)_03d.tif"
    # using nearest neighbor, retains more variation than bilinear  (see -r option)
    targetproj = split("-te 1 50.5 32 71.5 -tr 0.03 0.03 -t_srs EPSG:4326", ' ')
    options = split("--config GDAL_CACHEMAX 9999 -wm 9999 -co COMPRESS=ZSTD -co BIGTIFF=YES", ' ')
    @time run(`gdalwarp $(targetproj) $(options) NETCDF:"$(ncfile)"://$(varalt) $(outfile)`)
    return outfile
end

# not working because of incomplete metadata in many CORDEX files
function reproject_cordex_gdal(variable, altitude, year, hours=3)
    varalt = "$(variable)a$(altitude)m"
    bulkname = "EUR-11_MOHC-HadGEM2-ES_rcp85_r1i1p1_CNRM-ALADIN63_v1_$(hours)hr"
    ncfile = "D:/SMHI/$(varalt)_$(bulkname)_$(year)01010$(hours)00-$(year+1)01010000.nc"  # note not same hours of day, but no matter since it's the future
    outfile = "D:/SMHI/TESTcordex_$(varalt)_$(year)_10d3.tif"
    # using nearest neighbor, retains more variation than bilinear  (see -r option)
    proj4 = "+proj=lcca +lat_1=49.5 +lat_0=49.5 +lon_0=10.5 +k_0=1.0 +units=km" # +a=6371.220 +b=6371.220"
    targetproj = split("-te -10.5 34.5 32 71.5 -tr 0.10 0.10 -t_srs EPSG:4326", ' ')
    options = split("--config GDAL_CACHEMAX 9999 -wm 9999 -co COMPRESS=ZSTD -co BIGTIFF=YES", ' ')
    @time run(`gdalwarp -s_srs $(proj4) $(targetproj) $(options) NETCDF:"$(ncfile)"://$(varalt) $(outfile)`)
    return outfile
end

# gdalsrsinfo -o proj4 NETCDF:ua100m_NEU-3_ECMWF-ERAINT_evaluation_r1i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_3hr_200601010000-200612312100.nc://ua100m
# gdalsrsinfo -o proj4 NETCDF:va100m_EUR-11_MOHC-HadGEM2-ES_rcp85_r1i1p1_CNRM-ALADIN63_v1_3hr_205001010300-205101010000.nc://va100m

function readprojected(variable="v", altitude=100, year=2018)
    filename = "D:/SMHI/$(variable)a$(altitude)m_$(year)_03d.tif"
    ArchGDAL.readraster(filename)
end

AGread(filename) = ArchGDAL.readraster(filename)[:,:,:]
filesizeMB(filename) = round(filesize(filename)/1024^2, digits=1)

function readalldata(filename)
    ArchGDAL.readraster(filename) do dataset
        return dataset[:,:,:]
    end
end

function getGWArangesfromSMHI(variable="v", altitude=100, year=2018, method="")
    suffix = isempty(method) ? "" : "_$method"
    filename = "D:/SMHI/$(variable)a$(altitude)m_$(year)_03d$(suffix).tif"
    ArchGDAL.read(filename) do dataset
        geotransform = ArchGDAL.getgeotransform(dataset)
        size = convert.(Int, (ArchGDAL.width(dataset), ArchGDAL.height(dataset)))
        coordextent = getextent(geotransform, size)
        # display(size)
        # display(coordextent)
        latrange_gwa, lonrange_gwa = bbox2ranges(extent2bbox(coordextent), 100)
        return lonrange_gwa, latrange_gwa
    end
end

function plotSMHImaps(altitude=100, year=2018; method="")
    suffix = "$(altitude)m$year$method"
    # smhi = getSMHIwind(altitude, year, method)  # reads mean wind
    drawmap(smhi; scalefactor=(1.0,1.8), colorrange=(0,12), save="smhi_$suffix.png")
    lonrange, latrange = getGWArangesfromSMHI("v", altitude, year, method)
    gwa_full = getwindatlas(altitude)[lonrange,latrange]
    gwa = gwa_full[2:3:end, 2:3:end]
    drawmap(gwa; scalefactor=(1.0,1.8), colorrange=(0,12), save="gwa_$suffix.png")
    drawmap(gwa_full; scalefactor=(1.0,1.8), colorrange=(0,12), save="gwa_full_$suffix.png")
    return gwa, smhi
end

function saveSMHIwind(altitude=100, year=2018; compress=4)
    u = readprojected("u", altitude, year)[:,:,:]
    v = readprojected("v", altitude, year)[:,:,:]
    wind = sqrt.(u.^2 + v.^2)
    filename = "D:/SMHI/windSMHI_$(altitude)m$(year).h5"
    gridsize = size(u)[1:2]
    hours = size(u, 3)
    # saveTIFF(wind, filename, [1.0, 50.5, 31.99, 71.5])
    # similar_dataset(origdataset, wind, filename; compressmethod="ZSTD")
    h5open(filename, "w") do file 
        group = file["/"]
        dataset_wind = create_dataset(group, "wind", datatype(Float32), dataspace(size(u)...); chunk=(16,16,hours), compress)
        dataset_meanwind = create_dataset(group, "meanwind", datatype(Float32), dataspace(gridsize...); chunk=gridsize, compress)
        dataset_wind[:,:,:] = wind
        dataset_meanwind[:,:] = meandrop(wind, dims=3)
    end
end

testfilename_hclim(variable, altitude, year) =
    "D:/SMHI/$(variable)a$(altitude)m_NEU-3_ECMWF-ERAINT_evaluation_r1i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_3hr_$(year)01010000-$(year)12312100.nc"

function metadataSMHI(variable="v", altitude=100, year=2018)
    fn = testfilename_hclim(variable, altitude, year)
    nc = Dataset(fn)
    display(nc)
    # gdalinfo_path() do gdalinfo
    #     # run(`$gdalinfo NETCDF:$fn`)
    #     run(`$gdalinfo --formats`)
    # end
    proj = nc["Lambert_Conformal"].attrib["proj4"]
    # nc["lon"][1,1], nc["lon"][1,end], nc["lon"][end,1], nc["lon"][end,end]
    nc["lon"][:,:], nc["lat"][:,:]
end

function get_minimum_lonlat_spacing()
    lon, lat = metadataSMHI()
    get_minimum_lonlat_spacing(lon, lat)
end

function get_minimum_lonlat_spacing(lon, lat)
    dist = similar(lon)
    rows, cols = size(lon)
    for r=1:rows, c=1:cols
        dist[r,c] = sqrt(minimum((lon[rr,cc] - lon[r,c])^2 + (lat[rr,cc] - lat[r,c])^2
                                for rr = max(1,r-1):min(rows,r+1) for cc = max(1,c-1):min(cols,c+1)
                                    if rr != r || cc != c))
    end
    return dist
end

function testSMHI(variable="v", altitude=100, year=2018)
    fn = testfilename_hclim(variable, altitude, year)
    nc = Dataset(fn)
    proj = nc["Lambert_Conformal"].attrib["proj4"]
    var = "$(variable)a$(altitude)m"

    # infile = in_datafolder("gwa3_250_wind-speed_$(altitude)m.tif")
    # gdalinfo_path() do gdalinfo
    #     run(`$gdalinfo $infile`)
    # end
    # println("\n")
    # outfile = in_datafolder("Global Wind Atlas v3 - $(altitude)m wind speed.tif")
    # options = split("-r bilinear -te -180 -90 180 90 -tr 0.01 0.01", ' ')
    # gdalwarp_path() do gdalwarp
    #     @time run(`$gdalwarp $options -co COMPRESS=LZW $infile $outfile`)
    # end

    # options = split("-r bilinear -te -180 -90 180 90 -tr 0.01 0.01", ' ')
    # ds = ArchGDAL.gdalwarp(Float32.(nc[var][:,:,23]), options) do warped
    #     band = ArchGDAL.getband(warped, 1)
    #     ArchGDAL.read(band)
    # end

    # TimeMem-1.0 gdalwarp -te 1 50 36 70 -tr 0.03 0.03 -t_srs EPSG:4326 NETCDF:"ua200m_NEU-3_ECMWF-ERAINT_evaluation_r1i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_3hr_201801010000-201812312100.nc"://ua200m --config GDAL_CACHEMAX 9999 -wm 9999 -co COMPRESS=ZSTD -co BIGTIFF=YES ua200m_2018_03d.tif
    # gdal_translate -b 1 -b 2 -b 3 --config GDAL_CACHEMAX 9999 -co COMPRESS=ZSTD -co BIGTIFF=YES testag.tif testag_small.tif 

    # gdalwarp -of netCDF -te 1 51 32 70  -s_srs "+proj=lcca +lat_1=62.200000 +lat_0=62.200000 +lon_0=10.000000 +k_0=1.0 +x_0=669102.401911 +y_0=1265895.127345 +a=6371220.000000 +b=6371220.000000" -t_srs EPSG:4326 HDF5:"va200m_NEU-3_ECMWF-ERAINT_evaluation_r1i1p1_HCLIMcom-HCLIM38-AROME_x2yn2v1_3hr_201801010000-201812312100.nc"://va200m -to SRC_METHOD=NO_GEOTRANSFORM test.nc

    # return ds
    Float32.(nc[var][:,:,:])
end

# Can optionally zero cells that are zero in the Global Wind Atlas to save a lot of disk space.
function readSMHI(; year=2018, windatlas_only=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")
    
    filename = joinpath(datafolder, "era5wind$year.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    windatlas = reshape(imresize(getwindatlas(), gridsize), (1,gridsize...))

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        dataset_wind = create_dataset(group, "wind", datatype(Float32), dataspace(hours,gridsize...), chunk=(hours,16,16), blosc=3)
        dataset_meanwind = create_dataset(group, "meanwind", datatype(Float32), dataspace(gridsize...), chunk=gridsize, blosc=3)

        totalwind = zeros(gridsize)
        hour = 1

        count = 0
        for month = 1:12, monthhalf = 1:2
            if monthhalf == 1
                firstday, lastday = "01", "15"
            else
                firstday = "16"
                lastday = Dates.daysinmonth(Date("$year-$month"))
            end
            monthstr = lpad(month,2,'0')
            date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
            erafile = joinpath(downloadsfolder, "wind$year-$monthstr$firstday-$monthstr$lastday.nc")

            count += 1
            println("\nFile $count of 24:")
            println("Reading wind components from $erafile...")
            # Permute dimensions to get hours as dimension 1 (for efficient iteration in GISwind())
            ncdataset = Dataset(erafile)
            u100 = permutedims(ncdataset["u100"][:,:,:], [3,1,2])
            v100 = permutedims(ncdataset["v100"][:,:,:], [3,1,2])

            println("Calculating absolute speed...")
            wind = replace(sqrt.(u100.^2 + v100.^2), missing => 0.0) .* (windatlas .> 0)

            totalwind = totalwind + sumdrop(wind, dims=1)
            len = size(wind,1)
            println("Writing to $filename...")
            dataset_wind[hour:hour+len-1,:,:] = wind
            hour += len
        end
        println("\nWriting meanwind to $filename...")
        dataset_meanwind[:,:] = totalwind/hours
    end
    nothing
end

function copy_gtiff(infile, outfile; compressmethod="ZSTD")
    ArchGDAL.readraster(infile) do origdataset
        nbands = ArchGDAL.nraster(origdataset)
        band1 = ArchGDAL.getband(origdataset, 1)
        nodata = ArchGDAL.getnodatavalue(band1)
        isfile(outfile) && rm(outfile)
        origdata = origdataset[:,:,:]
        ArchGDAL.create(outfile,
            driver = ArchGDAL.getdriver("GTiff"),
            width = ArchGDAL.width(origdataset),
            height = ArchGDAL.height(origdataset),
            nbands = nbands,
            dtype = ArchGDAL.pixeltype(band1),
            options = ["BIGTIFF=YES", "COMPRESS=$compressmethod"]
        ) do dataset
            geotransform = ArchGDAL.getgeotransform(origdataset)
            proj = ArchGDAL.getproj(origdataset)
            ArchGDAL.setgeotransform!(dataset, geotransform)
            ArchGDAL.setproj!(dataset, proj)
            for b = 1:nbands
                band = ArchGDAL.getband(dataset, b)
                ArchGDAL.setnodatavalue!(band, nodata)
                ArchGDAL.write!(band, origdata[:,:,b])
            end
        end
        return origdata
    end
end

function recompress_gtiff(infile, outfile; compressmethod="ZSTD")
    gdal_translate_path() do gdal_translate
        run(`$gdal_translate -q -co BIGTIFF=YES -co COMPRESS=$compressmethod $infile $outfile`)
    end
end

function save_dataset_with_new_data(origdataset, newdata, outfile; compressmethod="ZSTD")
    nbands = ArchGDAL.nraster(origdataset)
    ArchGDAL.create(outfile,
        driver = ArchGDAL.getdriver("GTiff"),
        width = ArchGDAL.width(origdataset),
        height = ArchGDAL.height(origdataset),
        nbands = nbands,
        dtype = Float32,
        options = ["BIGTIFF=YES", "COMPRESS=$compressmethod"]
    ) do dataset
        geotransform = ArchGDAL.getgeotransform(origdataset)
        proj = ArchGDAL.getproj(origdataset)
        ArchGDAL.setgeotransform!(dataset, geotransform)
        ArchGDAL.setproj!(dataset, proj)
        nodata = ArchGDAL.getnodatavalue(ArchGDAL.getband(origdataset, 1))
        for b = 1:nbands
            band = ArchGDAL.getband(dataset, b)
            ArchGDAL.setnodatavalue!(band, nodata)
            ArchGDAL.write!(band, newdata[:,:,b])
        end
    end
    nothing
end

function save_and_recompress(data, filename, extent)
    ext = splitext(filename)[2]
    !in(ext, [".tif", ".h5"]) && error("Extension $ext not recognized.")
    println("Saving as $filename...")
    if ext == ".h5"
        @time h5open(filename, "w") do file 
            group = file["/"]
            dataset_data = create_dataset(group, "wind", datatype(Float32), dataspace(size(data)), chunk=(16,16,size(data,3)), blosc=3)
            dataset_extent = create_dataset(group, "extent", datatype(Float64), dataspace(size(extent)))
            dataset_data[:,:,:] = data
            dataset_extent[:] = extent
        end
    else
        @time saveTIFF(data, filename, extent, compressmethod="ZSTD")
        if filesize(filename) > 2^32
            println("Recompressing...")
            temp = tempname(dirname(filename))
            @time recompress_gtiff(filename, temp)
            mv(temp, filename, force=true)
        end
    end
    nothing
end

function read_gtiff(infile)
    ArchGDAL.readraster(infile) do origdataset
        nbands = ArchGDAL.nraster(origdataset)
        band1 = ArchGDAL.getband(origdataset, 1)
        nodata = ArchGDAL.getnodatavalue(band1)
        origdata = origdataset[:,:,:]
        width = ArchGDAL.width(origdataset)
        height = ArchGDAL.height(origdataset)
        dtype = ArchGDAL.pixeltype(band1)
        geotransform = ArchGDAL.getgeotransform(origdataset)
        proj = ArchGDAL.getproj(origdataset)
        @show(nbands, band1, nodata, width, height, dtype, geotransform, proj)
        return origdata
    end
end

# function XY_Lambert_Conformal_Conical_2_Lat_Lon(x,y,λ0,ϕ0,ϕ1)
#     R = 6371 # Radius of Earth
#     λ0_rad = λ0*π/180
#     ϕ0_rad = ϕ0*π/180
#     ϕ1_rad = ϕ1*π/180
#     x_skaliert = x/R
#     y_skaliert = y/R
#     n = sin(ϕ1_rad)
#     F = (cos(ϕ1_rad)*(tan(π/4 + ϕ1_rad/2))^n)/n
#     ρ0 = F*(cot(π/4 + ϕ0_rad/2))^n
#     θ = atan(x_skaliert/(ρ0 - y_skaliert))
#     ρ = sign(n)*sqrt(x_skaliert^2 + (ρ0 - y_skaliert)^2)
#     ϕ_rad = 2atan((F/ρ)^(1/n)) - π/2
#     λ_rad = λ0_rad + θ/n
#     ϕ = ϕ_rad*180/π
#     λ = λ_rad*180/π
#     return λ,ϕ
# end
