using DataDeps

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
        ("WDPA (protected areas):   [add code to check for monthly updates]", "WDPA_Feb2020.zip",
            "https://www.protectedplanet.net/downloads/WDPA_Feb2020?type=shapefile"),
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
        ("Synthetic demand input data", "synthetic_demand_data.zip",
            "https://chalmersuniversity.box.com/shared/static/2xe2tdo4ylmldv7e2mc8p0fd1x15gcis.zip"),
        ("IEA electricity demand 2016", "ieademand_2016.csv",
            "https://chalmersuniversity.box.com/shared/static/1gz9nkhgteflunbinom7u8pc1r5qjglu.csv"),
        ("SSP v2 Final Energy/Electricity from IAM scenarios", "SSP v2 Final Energy - Electricity.csv",
            "https://chalmersuniversity.box.com/shared/static/hfrt1uq8slgucj3m22rv1j3yhe20d2sa.csv")
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

    syntheticdemandfolder = joinpath(datafolder, "syntheticdemand")
    mkpath(syntheticdemandfolder)
    mkpath(joinpath(syntheticdemandfolder, "output"))
    cp(joinpath(dirname(@__FILE__), "synthetic_demand_demo.py"), joinpath(syntheticdemandfolder, "synthetic_demand_demo.py"))

    println("\n\n\nCleaning up...")

    mv(joinpath(datafolder, "ETOPO1_Ice_c_geotiff", "ETOPO1_Ice_c_geotiff.tif"), joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"))
    rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff"))
    mv(joinpath(datafolder, "temp_popgdp", "data"), joinpath(datafolder, "global_population_and_gdp"))
    rm(joinpath(datafolder, "temp_popgdp"))
    mv(joinpath(datafolder, "temp_ssp2", "SSP2_1km"), joinpath(datafolder, "SSP2_1km"))
    rm(joinpath(datafolder, "temp_ssp2"), force=true, recursive=true)
    rm(joinpath(datafolder, "SSP2_1km", "PaxHeaders.62069"), force=true, recursive=true)
    rm(joinpath(datafolder, "nuts-2016-01m.shp"), force=true, recursive=true)
    mv(joinpath(datafolder, "WRI - Global Power Plant Database v1.10"), joinpath(datafolder, "tempWRI"))
    mv(joinpath(datafolder, "tempWRI", "WRI - Global Power Plant Database v1.10"),
        joinpath(datafolder, "WRI - Global Power Plant Database v1.10"))
    rm(joinpath(datafolder, "tempWRI"))
    mv(joinpath(datafolder, "synthetic_demand_data", "data"), joinpath(syntheticdemandfolder, "data"))
    rm(joinpath(datafolder, "synthetic_demand_data"))

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
