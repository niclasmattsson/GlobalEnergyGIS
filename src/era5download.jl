using PyCall, Pkg.TOML

export era5download, monthlyera5download, saveconfig, download_and_convert_era5

getconfig(key) = getconfig()[key]

function getconfig()
    configfile = joinpath(homedir(), ".GlobalEnergyGIS_config")
    if !isfile(configfile)
        error("Configuration file missing, please run saveconfig(datafolder, uid, api_key) first. See GlobalEnergyGIS README.")
    end
    return TOML.parsefile(configfile)
end

function saveconfig(datafolder::AbstractString, uid::Int, api_key::AbstractString; agree_terms=false)
    !agree_terms && error("You must agree to the terms of use of all datasets to proceed. See GlobalEnergyGIS README.")
    downloadsfolder = joinpath(datafolder, "downloads")
    mkpath(downloadsfolder)
    configfile = joinpath(homedir(), ".GlobalEnergyGIS_config")
    open(configfile, "w") do io
        d = Dict("datafolder"=>datafolder, "agree_terms"=>agree_terms)
        TOML.print(io, d)
        println("Configuration file written to $configfile.")
    end
    cds_id(uid, api_key)
end

function cds_id(uid::Int, api_key::AbstractString)
    filename = joinpath(homedir(), ".cdsapirc")
    isfile(filename) && error("$filename already exists, no changes made. Please check its contents manually.")
    open(filename, "w") do file
        write(file, "url: https://cds.climate.copernicus.eu/api/v2\n")
        write(file, "key: $uid:$api_key\n")
        println("Copernicus credentials written to $filename.")
    end
end

function download_and_convert_era5(year=2018; datasets=["wind", "solar", "temp"])
    for dataset in datasets
        println("\nDownloading ERA5 $dataset data from Copernicus...")
        era5download(year; datasets=[dataset])
        println("\nConverting downloaded $dataset data to HDF5 and recompressing...")
        if dataset == "solar"
            makesolarera5(; year)
        elseif dataset == "wind"
            makewindera5(; year)
        elseif dataset == "temp"
            maketempera5(; year)
        end
        println("\nCleanup: deleting downloaded $dataset data...")
        clearvars_era5(; year, datasets=[dataset])
    end
    println("\nERA5 datasets $datasets downloaded and converted. Temporary files cleaned up.")
end

function era5download(year=2018; datasets=["wind", "solar", "temp"])
    mkpath(in_datafolder("downloads"))
    count = 0
    for dataset in datasets, month = 1:12, monthhalf = 1:2
        if dataset == "wind"
            vars = ["100m_u_component_of_wind", "100m_v_component_of_wind"]
        elseif dataset == "solar"
            vars = ["surface_solar_radiation_downwards", "total_sky_direct_solar_radiation_at_surface"]
        else
            vars = ["2m_temperature"]
        end
        if monthhalf == 1
            firstday, lastday = "01", "15"
        else
            firstday = "16"
            lastday = Dates.daysinmonth(Date("$year-$month"))
        end
        monthstr = lpad(month,2,'0')
        date1, date2 = "$year-$monthstr-$firstday", "$year-$monthstr-$lastday"
        outfile = in_datafolder("downloads", "$dataset$year-$monthstr$firstday-$monthstr$lastday.nc")
        count += 1
        println("\nFile $count of $(24*length(datasets)):")
        request_era5_vars(outfile, vars, date1, date2)
    end
end

function monthlyera5download(; datasets=["wind", "solar", "temp"])
    mkpath(in_datafolder("downloads"))
    for dataset in datasets
        if dataset == "wind"
            vars = ["100m_u_component_of_wind", "100m_v_component_of_wind"]
        elseif dataset == "solar"
            vars = ["surface_solar_radiation_downwards", "total_sky_direct_solar_radiation_at_surface"]
        else
            vars = ["2m_temperature"]
        end
        # avoid 2020 to avoid partial-year problems and near real-time ERA5T data
        # (which includes an additional dimension to indicate ERA5/ERA5T)
        years = 1979:2019
        outfile = in_datafolder("downloads", "monthly$(dataset)_$(years[1])-$(years[end]).nc")
        request_monthly_era5_vars(outfile, vars, collect(years))
    end
end

function request_era5_vars(outfile::String, vars::Vector{String}, firstdate::String, lastdate::String)
    datestring = "$firstdate/$lastdate"
    py"""
    import cdsapi

    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': $vars,
            'grid': '0.28125/0.28125',
            'area': '89.859375/-179.859375/-89.859375/179.859375',
            'date': $datestring,
            'time': '00/to/23/by/1'
        },
        $outfile)
    """
end

function request_monthly_era5_vars(outfile::String, vars::Vector{String}, years::Vector{Int})
    stryears = string.(years)
    months = string.(1:12)
    py"""
    import cdsapi

    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'format': 'netcdf',
            'variable': $vars,
            'grid': '0.28125/0.28125',
            'area': '89.859375/-179.859375/-89.859375/179.859375',
            'year': $stryears,
            'month': $months,
            'time': '00:00'
        },
        $outfile)
    """
end

# For some reason the delivered NetCDF files are unreadable unless they are limited to 15-16 days each.

# default resolution:  'grid': '0.25/0.25'  (but the max ERA5 resolution is 0.28125 degrees, so use that instead)
# default area:        'area': '90/0/-90/179.75'  (note max latitude is 90)

# The default grid size using a resolution of 0.25 is 721x1440 pixels.
# The default grid size using a resolution of 0.28125 is 641x1280 pixels.

# 'area': '89.859375/-179.859375/-89.859375/179.859375': get correct center points for a 640x1280 grid (north/west/south/east)

# The CDS server doesn't use high enough precision to evaluate the area coordinates exactly. Each lon & lat is off
# by about 0.000625 degrees, but this error is very small compared to the resolution (0.28125 degrees). Otherwise
# using the default grid we would have to interpolate to get correct center points, and this would also lose some
# precision. This way we avoid doing all that computational work for a 641x1280x8760 grid.

# Note that using the default grid with a resolution of 0.28125 also makes longitude errors of up to 0.000244 degrees
# (the numbers come from inspecting latitude & longitude coordinates of ERA5 test data downloads).

# https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows
# https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key
# https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5
# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
# https://cds.climate.copernicus.eu/cdsapp#!/yourrequests