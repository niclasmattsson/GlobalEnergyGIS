using DataDeps

export download_datasets

function download_datasets()
    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")

    println("\nDownloading Global Wind Atlas:")
    download_progressbar("https://chalmersuniversity.box.com/shared/static/25cjcah213sk4wdkqwi8t8kxxoy3hh6k.tif",
                        joinpath(datafolder, "Global Wind Atlas v1 - 100m wind speed.tif"))

    println("\nDownloading WDPA (protected areas):")
    download_progressbar("https://www.protectedplanet.net/downloads/WDPA_Nov2019?type=shapefile",
                        joinpath(downloadsfolder, "WDPA_Nov2019.shp"))

    println("\nDownloading GADM (global administrative areas):")
    zipfile = joinpath(downloadsfolder, "gadm36.zip")
    download_progressbar("https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_shp.zip",
                        zipfile)
    println("\nUnzipping zip archive...")
    unpack(zipfile, joinpath(downloadsfolder, "GADM36"), ".zip")

    println("\nDownloading NUTS (administrative areas in Europe):")
    zipfile = joinpath(downloadsfolder, "nuts-2016-01m.shp.zip")
    download_progressbar("https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2016-01m.shp.zip",
                        zipfile)
    println("\nUnzipping zip archive...")
    unpack(zipfile, joinpath(downloadsfolder, "NUTS2016"), ".zip")

    println("\nDownloading Land Cover:")
    download_progressbar("https://chalmersuniversity.box.com/shared/static/zm8zdkv1mz0wns0u77afkna95wbo7ggo.tif",
                        joinpath(datafolder, "Landcover - USGS MODIS.tif"))

    println("\nDownloading ETOPO1 Topography:")
    zipfile = joinpath(downloadsfolder, "ETOPO1_Ice_c_geotiff.zip")
    download_progressbar("https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip",
                        zipfile)
    println("\nUnzipping zip archive...")
    unpack(zipfile, joinpath(downloadsfolder, "ETOPO1"), ".zip")

    println("\nDownloading 1-km downscaled spatial population scenarios (Gao et al):")
    tarfile = joinpath(downloadsfolder, "SSP2_1km_netcdf.tar")
    download_progressbar("http://www.cgd.ucar.edu/iam/modeling/data/SSP2_1km_netcdf.tar",
                        tarfile)
    println("\nExtracting tar archive...")
    unpack(tarfile, joinpath(downloadsfolder, "SSP2_1km"), ".tar")

    println("\nDownloading global population & GDP:")
    zipfile = joinpath(downloadsfolder, "global_population_and_gdp.zip")
    download_progressbar("https://osf.io/7jv3n/download",
                        zipfile)
    println("\nUnzipping zip archive...")
    unpack(zipfile, joinpath(downloadsfolder, "global_population_and_gdp"), ".zip")

    println("\nDownloads complete.")
end

function download_progressbar(url::AbstractString, filename::AbstractString)
    println("Downloading to $filename, please wait...")
    curl = Base.find_curl()
    if curl == nothing
        # if curl is not on the system, fall back to a download without progress bar
        download(url, filename)
    else
        download_curl(curl, url, filename)
    end
end

function download_curl(curl_exe::AbstractString, url::AbstractString, filename::AbstractString)
    run(`$curl_exe --progress-bar -g -L -f -o $filename $url`, wait=true)
    return filename
end

unpack(inputfilename, outputpath, extension) = run(DataDeps.unpack_cmd(inputfilename, outputpath, extension, ""))
