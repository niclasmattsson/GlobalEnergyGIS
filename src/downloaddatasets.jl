using DataDeps

export download_datasets

function download_datasets(startfile=1)
    datafolder = getconfig("datafolder")

    # tuples of (dataset_name, filename, url)
    datasets = [
        ("Global Wind Atlas", "Global Wind Atlas v1 - 100m wind speed.tif",
            "https://chalmersuniversity.box.com/shared/static/25cjcah213sk4wdkqwi8t8kxxoy3hh6k.tif"),
        ("WDPA (protected areas):   [add code to check for monthly updates]", "WDPA_Dec2019.shp",
            "https://www.protectedplanet.net/downloads/WDPA_Dec2019?type=shapefile"),
        ("GADM (global administrative areas)", "gadm36.zip",
            "https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_shp.zip"),
        ("NUTS (administrative areas in Europe)", "nuts-2016-01m.shp.zip",
            "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2016-01m.shp.zip"),
        ("Land Cover", "Landcover - USGS MODIS.tif",
            "https://chalmersuniversity.box.com/shared/static/zm8zdkv1mz0wns0u77afkna95wbo7ggo.tif"),
        ("ETOPO1 Topography", "ETOPO1_Ice_c_geotiff.zip",
            "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip"),
        ("1-km downscaled spatial population scenarios (Gao et al)", "temp_ssp2.tar",
            "http://www.cgd.ucar.edu/iam/modeling/data/SSP2_1km_netcdf.tar"),
        ("Global population & GDP", "temp_popgdp.zip",
            "https://osf.io/7jv3n/download")
    ]


    for (i, datasetinfo) in enumerate(datasets)
        i < startfile && continue
        name, filename, url = datasetinfo

        println("\nDownloading dataset $i: $name")
        download_progressbar(url, joinpath(datafolder, filename))
    end

    println("\nDownloads complete.")


    for datasetinfo in datasets
        name, filename, url = datasetinfo
        foldername, extension = splitext(filename)

        if extension in [".zip", ".tar"]
            println("\nUnpacking archive: $filename")
            unpack(joinpath(datafolder, filename), joinpath(datafolder, foldername), extension)
        end
    end

    filename = "NUTS_RG_01M_2016_4326_LEVL_3.shp.zip"
    println("\nUnpacking archive: $filename")
    unpack(joinpath(datafolder, "nuts-2016-01m.shp", filename), joinpath(datafolder, "nuts2016-level3"), ".zip")

    println("\n\n\nCleaning up...")

    mv(joinpath(datafolder, "ETOPO1_Ice_c_geotiff", "ETOPO1_Ice_c_geotiff.tif"), joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"))
    rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff"))
    mv(joinpath(datafolder, "temp_popgdp", "data"), joinpath(datafolder, "global_population_and_gdp"))
    rm(joinpath(datafolder, "temp_popgdp"))
    mv(joinpath(datafolder, "temp_ssp2", "SSP2_1km"), joinpath(datafolder, "SSP2_1km"))
    rm(joinpath(datafolder, "temp_ssp2"), force=true, recursive=true)
    rm(joinpath(datafolder, "SSP2_1km", "PaxHeaders.62069"), force=true, recursive=true)
    rm(joinpath(datafolder, "nuts-2016-01m.shp"), force=true, recursive=true)

    for datasetinfo in datasets
        name, filename, url = datasetinfo
        foldername, extension = splitext(filename)

        if extension in [".zip", ".tar"]
            rm(joinpath(datafolder, filename))
        end
    end
    println("\n\nDone.")
end

function download_progressbar(url::AbstractString, filename::AbstractString)
    println("Downloading to $filename...")
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
