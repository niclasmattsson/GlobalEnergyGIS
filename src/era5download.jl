using PyCall

export era5download, cds_id

function cds_id(uid::Int, api_key::String)
    filename = joinpath(homedir(), ".cdsapirc")
    isfile(filename) && error("$filename already exists, no changes made. Please check its contents manually.")
    open(filename, "w") do file
        write(file, "url: https://cds.climate.copernicus.eu/api/v2\n")
        write(file, "key: $uid:$api_key\n")
        println("Credentials written to $filename.")
    end
end

function era5download(year)
    for dataset in ["wind", "solar"], month = 1:12, monthhalf = 1:2       # ["wind", "solar"]
        if dataset == "wind"
            var1, var2 = "100m_u_component_of_wind", "100m_v_component_of_wind"
        else
            var1, var2 = "surface_solar_radiation_downwards", "total_sky_direct_solar_radiation_at_surface"
        end
        if monthhalf == 1
            firstday, lastday = "01", "15"
        else
            firstday = "16"
            lastday = Dates.daysinmonth(Date("$year-$month"))
        end
        monthstr = lpad(month,2,'0')
        date1, date2 = "$year-$monthstr-$firstday", "$year-$monthstr-$lastday"
        outfile = "D:/testera5/$dataset$year-$monthstr$firstday-$monthstr$lastday.nc"
        request_era5_vars(outfile, [var1, var2], date1, date2)
    end
end

function request_era5_vars(outfile::String, vars::Vector{String}, firstdate::String, lastdate::String)
    println("Something seems wrong with outfile in Python code...")
    # vars = ["'$v'" for v in vars]
    # varstring = join(vars, ", ")
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
        '$outfile')
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