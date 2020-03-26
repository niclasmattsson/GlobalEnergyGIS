using BinDeps

export download_datasets

function download_datasets(startfile=1)
    datafolder = getconfig("datafolder")
    if !isfile(datafolder)
        mkpath(datafolder)
    end

    # tuples of (dataset_name, filename, url)
    datasets = [
        ("Global Wind Atlas", "Global Wind Atlas v1 - 100m wind speed.tif",
            "https://chalmersuniversity.box.com/shared/static/25cjcah213sk4wdkqwi8t8kxxoy3hh6k.tif"),
        ("WDPA (protected areas):   [add code to check for monthly updates]", "WDPA_Mar2020.zip",
            "https://www.protectedplanet.net/downloads/WDPA_Mar2020?type=shapefile"),
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
            "https://osf.io/7jv3n/download"),
        ("WRI Global Power Plant Database", "WRI - Global Power Plant Database v1.10.zip",
            "https://chalmersuniversity.box.com/shared/static/ss6gycw7hf10e1fiicbxl5rgk08q5xr9.zip"),
            # switch to the official v1.2 link later (some plant cleanup is hardcoded for v1.1) 
            # "http://datasets.wri.org/dataset/540dcf46-f287-47ac-985d-269b04bea4c6/resource/c240ed2e-1190-4d7e-b1da-c66b72e08858/download/globalpowerplantdatabasev120"),
        ("Hydropower - future potential (Gernaat et al)", "Hydro database (Gernaat) - potential.csv",
            "https://chalmersuniversity.box.com/shared/static/uxsa2wfj7o6s6upknjagnhwfa88aw3cx.csv"),
        ("Hydropower - existing capacity (Gernaat et al)", "Hydro database (Gernaat) - existing (GRanD).csv",
            "https://chalmersuniversity.box.com/shared/static/cnbqvt1oxnzyky0gb5k978ki9xp1ea23.csv"),
        ("Country-level hydropower capacity (WEC)", "WEC hydro capacity 2015.csv",
            "https://chalmersuniversity.box.com/shared/static/0kmxcatjzt8ivo4ktw65yoybts2uht76.csv"),
        ("Country-level hydropower potential (WEC)", "WEC hydro potentials.csv",
            "https://chalmersuniversity.box.com/shared/static/czkpqr0b6572bzdslyn76gsrxvc87vvc.csv"),
        ("IEA electricity demand 2016", "ieademand_2016.csv",
            "https://chalmersuniversity.box.com/shared/static/1gz9nkhgteflunbinom7u8pc1r5qjglu.csv"),
        ("SSP v2 Final Energy/Electricity from IAM scenarios", "SSP v2 Final Energy - Electricity.csv",
            "https://chalmersuniversity.box.com/shared/static/hfrt1uq8slgucj3m22rv1j3yhe20d2sa.csv"),
        ("Synthetic demand training data", "syntheticdemand_trainingdata.csv",
            "https://chalmersuniversity.box.com/shared/static/e1kfo9qioe87r34aq6693w1ed69rd2lu.csv"),
        ("Synthetic demand demand data", "syntheticdemand_demanddata.csv",
            "https://chalmersuniversity.box.com/shared/static/yuguu9c083pqjdf00ljw7hh9ekofg8t7.csv"),
        ("Synthetic demand time zone offsets", "syntheticdemand_timezoneoffsets.csv",
            "https://chalmersuniversity.box.com/shared/static/ei4v2453k90qguuwjggz1sw5tl2kavau.csv"),
        ("Time zone shape file", "timezones-with-oceans.shapefile.zip",
            "https://github.com/evansiroky/timezone-boundary-builder/releases/download/2019b/timezones-with-oceans.shapefile.zip")
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
        fullpath = joinpath(datafolder, filename)

        if extension in [".zip", ".tar"] && isfile(fullpath)
            println("\nUnpacking archive: $filename")
            unpack(fullpath, joinpath(datafolder, foldername), extension)
        end
    end

    if startfile <= 4
        filename = "NUTS_RG_01M_2016_4326_LEVL_3.shp.zip"
        println("\nUnpacking archive: $filename")
        unpack(joinpath(datafolder, "nuts-2016-01m.shp", filename), joinpath(datafolder, "nuts2016-level3"), ".zip")
        rm(joinpath(datafolder, "nuts-2016-01m.shp"), force=true, recursive=true)
    end

    println("\n\n\nCleaning up...")

    if startfile <= 6
        mv(joinpath(datafolder, "ETOPO1_Ice_c_geotiff", "ETOPO1_Ice_c_geotiff.tif"), joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"))
        rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff"))
    end
    if startfile <= 7
        mv(joinpath(datafolder, "temp_ssp2", "SSP2_1km"), joinpath(datafolder, "SSP2_1km"))
        rm(joinpath(datafolder, "temp_ssp2"), force=true, recursive=true)
        rm(joinpath(datafolder, "SSP2_1km", "PaxHeaders.62069"), force=true, recursive=true)
    end
    if startfile <= 8
        mv(joinpath(datafolder, "temp_popgdp", "data"), joinpath(datafolder, "global_population_and_gdp"))
        rm(joinpath(datafolder, "temp_popgdp"))
    end
    if startfile <= 9
        mv(joinpath(datafolder, "WRI - Global Power Plant Database v1.10"), joinpath(datafolder, "tempWRI"))
        mv(joinpath(datafolder, "tempWRI", "WRI - Global Power Plant Database v1.10"),
            joinpath(datafolder, "WRI - Global Power Plant Database v1.10"))
        rm(joinpath(datafolder, "tempWRI"))
    end

    for datasetinfo in datasets
        name, filename, url = datasetinfo
        foldername, extension = splitext(filename)

        if extension in [".zip", ".tar"]
            fullpath = joinpath(datafolder, filename)
            isfile(fullpath) && rm(fullpath)
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

unpack(inputfilename, outputpath, extension) = run(BinDeps.unpack_cmd(inputfilename, outputpath, extension, ""))
