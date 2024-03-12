using DataDeps, UrlDownload

export download_datasets

# tuples of (dataset_name, filename, url)
# not a const to avoid errors when updating urls
function get_dataset_info()
    [
    # https://globalwindatlas.info/api/gis/global/wind-speed/10
    # https://globalwindatlas.info/api/gis/global/wind-speed/50
    # https://globalwindatlas.info/api/gis/global/wind-speed/100
    # https://globalwindatlas.info/api/gis/global/wind-speed/150
    # https://globalwindatlas.info/api/gis/global/wind-speed/200
    "GWA100"       ("Global Wind Atlas", "Global Wind Atlas v3 - 100m wind speed.tif",
        "https://chalmersuniversity.box.com/shared/static/wfr6dm9bcmj0mcqtdn0uimhg0otd4ht1.tif")
    "GWA150"       ("Global Wind Atlas", "Global Wind Atlas v3 - 150m wind speed.tif",
        "https://chalmersuniversity.box.com/shared/static/ghexnwa7crukl58nkwric8v9kmlymvj2.tif")
    "GWA200"       ("Global Wind Atlas", "Global Wind Atlas v3 - 200m wind speed.tif",
        "https://chalmersuniversity.box.com/shared/static/7ib2jdni6hu9uwe1hp1mqoeeyqtg5601.tif")
    "WDPA"         ("WDPA (protected areas):", "WDPA.zip",
        "https://d1gam3xoknrgr2.cloudfront.net/current/$(filenameWDPA()).zip")
        # "https://chalmersuniversity.box.com/shared/static/wn1kznvy7qh1issqcxdlsq64kgtkaayi.zip")
    "GADM"         ("GADM (global administrative areas)", "gadm36.zip",
        "https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_shp.zip")
    "NUTS"         ("NUTS (administrative areas in Europe)", "nuts-2016-01m.shp.zip",
        "https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2016-01m.shp.zip")
    "landcover"    ("Land Cover", "Landcover - USGS MODIS.tif",
        "https://chalmersuniversity.box.com/shared/static/zm8zdkv1mz0wns0u77afkna95wbo7ggo.tif")
    "topography"   ("ETOPO1 Topography", "ETOPO1_Ice_c_geotiff.zip",
        "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/georeferenced_tiff/ETOPO1_Ice_c_geotiff.zip")
    "ssppop1"      ("1-km downscaled spatial population scenarios (Gao et al), file 1/3", "temp_ssp1.tar",
        "http://www.cgd.ucar.edu/iam/modeling/data/SSP1_1km_netcdf.tar")
    "ssppop2"      ("1-km downscaled spatial population scenarios (Gao et al), file 2/3", "temp_ssp2.tar",
        "http://www.cgd.ucar.edu/iam/modeling/data/SSP2_1km_netcdf.tar")
    "ssppop3"      ("1-km downscaled spatial population scenarios (Gao et al), file 3/3", "temp_ssp3.tar",
        "http://www.cgd.ucar.edu/iam/modeling/data/SSP3_1km_netcdf.tar")
    "gdppop"       ("Global population & GDP", "temp_popgdp.zip",
        "https://osf.io/7jv3n/download")
    "powerplants"  ("WRI Global Power Plant Database", "WRI - Global Power Plant Database v1.10.zip",
        "https://chalmersuniversity.box.com/shared/static/ss6gycw7hf10e1fiicbxl5rgk08q5xr9.zip")
        # switch to the official v1.2 link later (some plant cleanup is hardcoded for v1.1) 
        # "http://datasets.wri.org/dataset/540dcf46-f287-47ac-985d-269b04bea4c6/resource/c240ed2e-1190-4d7e-b1da-c66b72e08858/download/globalpowerplantdatabasev120")
    "timezones"    ("Time zone shape file", "timezones-with-oceans.shapefile.zip",
        "https://github.com/evansiroky/timezone-boundary-builder/releases/download/2019b/timezones-with-oceans.shapefile.zip")
    "monthlywind"  ("Average monthly wind speeds 1979-2019", "era5monthlywind.h5",
        "https://chalmersuniversity.box.com/shared/static/otvb5nq0lz5ntqx0e65kocos2afo7gas.h5")
    "monthlysolar" ("Average monthly solar insolation 1979-2019", "era5monthlysolar.h5",
        "https://chalmersuniversity.box.com/shared/static/jlpspmp9ou96hk7xno46rf79wg3d5hqc.h5")
    "various"      ("Various smaller datasets", "Various_smaller_datasets.zip",
        "https://chalmersuniversity.box.com/shared/static/w3pmx4xhgorgd6jejv23gn4ycsnza8s6.zip")
]
end

filenameWDPA(monthyear=Dates.format(now(), "uyyyy")) = "WDPA_WDOECM_$(monthyear)_Public_all_shp"

download_datasets(startfrom::Int) = download_datasets(get_dataset_info()[startfrom:end, 1]...)
download_datasets(datasetindices::AbstractVector) = download_datasets(get_dataset_info()[datasetindices, 1]...)

function download_datasets(shortnames::String...)
    datafolder = getconfig("datafolder")
    if !isfile(datafolder)
        mkpath(datafolder)
    end

    dataset_info = get_dataset_info()
    datasets = Dict(r[1] => r[2] for r in eachrow(dataset_info))
    shortnames = isempty(shortnames) ? dataset_info[:, 1] : shortnames

    for (i, shortname) in enumerate(shortnames)
        fullname, filename, url = datasets[shortname]

        println("\nDownloading dataset $i: $fullname")
        download_progressbar(url, joinpath(datafolder, filename))
        unpack_and_cleanup(shortname, filename, datafolder, dataset_info)
    end
    println("\nDownloads complete.")
end

function unpack_and_cleanup(shortname, filename, datafolder, dataset_info)
    foldername, extension = splitext(filename)
    fullpath = joinpath(datafolder, filename)

    if extension in [".zip", ".tar"] && isfile(fullpath)
        println("\nUnpacking archive: $filename")
        unpack(fullpath, joinpath(datafolder, foldername), extension)
    end

    function renameWDPAfiles(WDPAfolder)
        for filename in readdir(WDPAfolder)
            newname = replace(filename, filenameWDPA() => "WDPA-shapefile")
            if newname != filename 
                mv(joinpath(WDPAfolder, filename), joinpath(WDPAfolder, newname))
            end
        end
    end

    if shortname == "WDPA"
        renameWDPAfiles(joinpath(datafolder, "WDPA"))
        for i = 0:2
            foldername = "WDPA-shapefile_$i"
            println("\nUnpacking archive: $foldername.zip")
            unpack(joinpath(datafolder, "WDPA", "$foldername.zip"),
                    joinpath(datafolder, "WDPA", foldername), ".zip")
            renameWDPAfiles(joinpath(datafolder, "WDPA", foldername))
        end
    elseif shortname == "NUTS"
        filename = "NUTS_RG_01M_2016_4326_LEVL_3.shp.zip"
        println("\nUnpacking archive: $filename")
        unpack(joinpath(datafolder, "nuts-2016-01m.shp", filename), joinpath(datafolder, "nuts2016-level3"), ".zip")
        rm(joinpath(datafolder, "nuts-2016-01m.shp"), force=true, recursive=true)
    elseif shortname == "topography"
        mv(joinpath(datafolder, "ETOPO1_Ice_c_geotiff", "ETOPO1_Ice_c_geotiff.tif"), joinpath(datafolder, "ETOPO1_Ice_c_geotiff.tif"))
        rm(joinpath(datafolder, "ETOPO1_Ice_c_geotiff"))
    elseif shortname[1:end-1] == "ssppop"
        sspfolder = joinpath(datafolder, "SSP_1km")
        !isdir(sspfolder) && mkdir(sspfolder)
        for ssp = 1:3
            for y = 2010:10:2100
                sspfile = joinpath(datafolder, "temp_ssp$ssp", "SSP$(ssp)_1km", "ssp$(ssp)_total_$y.nc4")
                isfile(sspfile) && mv(sspfile, joinpath(sspfolder, "ssp$(ssp)_total_$y.nc4"))
            end
            rm(joinpath(datafolder, "temp_ssp$ssp"), force=true, recursive=true)
        end
    elseif shortname == "gdppop"
        mv(joinpath(datafolder, "temp_popgdp", "data"), joinpath(datafolder, "global_population_and_gdp"), force=true)
        rm(joinpath(datafolder, "temp_popgdp"))
    elseif shortname == "powerplants"
        mv(joinpath(datafolder, "WRI - Global Power Plant Database v1.10"), joinpath(datafolder, "tempWRI"), force=true)
        mv(joinpath(datafolder, "tempWRI", "WRI - Global Power Plant Database v1.10"),
            joinpath(datafolder, "WRI - Global Power Plant Database v1.10"))
        rm(joinpath(datafolder, "tempWRI"), force=true, recursive=true)
    elseif shortname == "various"
        mixeddir = joinpath(datafolder, "Various_smaller_datasets")
        for file in readdir(mixeddir)
            mv(joinpath(mixeddir, file), joinpath(datafolder, file), force=true)
        end
        rm(mixeddir)
    end

    if extension in [".zip", ".tar"] && isfile(fullpath)
        rm(fullpath)
    end
end

function download_progressbar(url::AbstractString, filename::AbstractString)
    println("Downloading to $filename...")
    # UrlDownload can take care of unpacking too, look into this.
    urldownload(url, true, compress=:none, parser=identity, save_raw=filename)
end

function unpack(inputfilename, outputpath, extension)
    !isdir(outputpath) && mkdir(outputpath)
    run(DataDeps.unpack_cmd(inputfilename, outputpath, extension, ""))
end
