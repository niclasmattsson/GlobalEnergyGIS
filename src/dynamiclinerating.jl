function dynamic_line_rating(year=2019)
    println("Reading line data and weather data for $year...")
    lines = read_line_data()
    weatherdata = read_weatherdata_DLR(year)

    println("Calculating line ratings...")
    ampacity = calculate_line_ratings(lines, weatherdata)
    return ampacity
    # add Global Wind Atlas data
end

getlinecoords(line) = (line.start_lon, line.start_lat), (line.end_lon, line.end_lat)

"""
Finds all grid cells intersected by a great circle path between two points (lon, lat) on a regular
lon-lat grid with specified resolution (in degrees). Returns path length and bearings within each cell.
"""
function greatcircle_waypoints(point1::Tuple, point2::Tuple, grid_resolution::Float64; step=0.0001)
    angular_dist = calc_angular_distance(point1, point2)

    # Determine number of steps for interpolation.
    # Step size is roughly half the grid resolution to ensure no grid cells are missed.
    num_steps = ceil(Int, angular_dist / step)

    path_points = interpolate_greatcircle_path(point1, point2, angular_dist, num_steps)
    cell_segments = group_segments_by_cell(path_points, grid_resolution)
    waypoints = calculate_lengths_and_bearings(cell_segments, point1)

    if !all_cells_adjacent(waypoints)
        @warn "Cells not adjacent, try a finer step size."
    end

    return waypoints
end

"Calculate angular distance between two points in degrees using Haversine formula"
function calc_angular_distance(point1::Tuple, point2::Tuple)
    lon1, lat1 = point1
    lon2, lat2 = point2

    Δlon = lon2 - lon1
    Δlat = lat2 - lat1
    a = sind(Δlat / 2)^2 + cosd(lat1) * cosd(lat2) * sind(Δlon / 2)^2
    angular_dist = 2 * asind(sqrt(a))

    return angular_dist
end

"Interpolate points along the great circle path"
function interpolate_greatcircle_path(point1, point2, angular_dist, num_steps)
    lon1, lat1 = point1
    lon2, lat2 = point2
    path_points = Vector{Tuple{Float64, Float64}}(undef, num_steps + 1)
    path_points[1] = point1

    for i in 1:num_steps
        f = i / num_steps
        
        A = sind((1 - f) * angular_dist) / sind(angular_dist)
        B = sind(f * angular_dist) / sind(angular_dist)

        x = A * cosd(lat1) * cosd(lon1) + B * cosd(lat2) * cosd(lon2)
        y = A * cosd(lat1) * sind(lon1) + B * cosd(lat2) * sind(lon2)
        z = A * sind(lat1) + B * sind(lat2)

        lon = atand(y, x)
        lat = atand(z, sqrt(x^2 + y^2))

        path_points[i+1] = (lon, lat)
    end
    return path_points
end

# Round a value to the nearest multiple of the given resolution, offset by resolution/2, treating -0.0 as 0.0.
function round_res(value, resolution)
    res2 = resolution / 2
    rounded = round((value - res2) / resolution) * resolution + res2
    return rounded == -0.0 ? 0.0 : rounded
end

"Group the path into segments belonging to each grid cell"
function group_segments_by_cell(path_points, grid_resolution)
    num_segments = length(path_points) - 1
    cell_segments = Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64}}}()

    for i in 1:num_segments
        p_start = path_points[i]
        p_end = path_points[i+1]
        
        # Determine grid cell for the midpoint of the segment
        mid_lon = (p_start[1] + p_end[1]) / 2
        mid_lat = (p_start[2] + p_end[2]) / 2

        cell_key = (round_res(mid_lon, grid_resolution), round_res(mid_lat, grid_resolution))
        if !haskey(cell_segments, cell_key)
            cell_segments[cell_key] = []
        end
        # Store the start and end points of the small segment
        push!(cell_segments[cell_key], p_start, p_end)
    end
    return cell_segments
end

"Calculate length and bearings for the path within each cell"
function calculate_lengths_and_bearings(cell_segments, start_point)
    results = @NamedTuple{cell::Tuple{Float64, Float64}, len::Float64, mean_bearing::Float64, bearing_error::Float64}[]
    for (cell, segments) in cell_segments
        entry_point, exit_point = segments[1], segments[end]

        len = greatcircledistance(entry_point, exit_point)
        bearings = greatcirclebearings(entry_point, exit_point)

        mean_bearing = mean(bearings)
        bearing_error = maximum(abs.(bearings .- mean_bearing))

        push!(results, (; cell, len, mean_bearing, bearing_error))
    end
    sort!(results, by=x->greatcircledistance(start_point, cell_segments[x.cell][1]))
    return results
end

"Calculate length and bearings for the path within each cell"
function calculate_lengths_and_bearings_alt(cell_segments)
    results = map(collect(cell_segments)) do (cell, segments)
        entry_point, exit_point = segments[1], segments[end]

        len = greatcircledistance(entry_point, exit_point)
        bearings = greatcirclebearings(entry_point, exit_point)

        mean_bearing = mean(bearings)
        bearing_error = maximum(abs.(bearings .- mean_bearing))

        (; cell, len, mean_bearing, bearing_error)
    end
    sort!(results, by=x->greatcircledistance(start_point, cell_segments[x.cell][1]))
    return results
end

adjacentcells(cell1, cell2) = xor(cell1[1] == cell2[1], cell1[2] == cell2[2])
all_cells_adjacent(waypoints) = all(adjacentcells(waypoints[i].cell, waypoints[i+1].cell) for i = 1:length(waypoints)-1)

function download_era5_DLR(year)
    date1, date2 = "$year-01-01", "$year-12-31"

    # Split into two requests, one for instantaneous variables (wind/temp)
    # and one for accumulated (solar) (otherwise Copernicus returns a zip file)
    windvars = ["100m_u_component_of_wind", "100m_v_component_of_wind", "2m_temperature"]
    outfile = in_datafolder("DLR", "ehubDLR_windtemp_$year.nc")
    request_era5_vars(outfile, windvars, date1, date2; res=0.25, bbox=(3.5, 33.0, 53.5, 72.5))

    solarvars = ["surface_solar_radiation_downwards", "total_sky_direct_solar_radiation_at_surface"]
    outfile = in_datafolder("DLR", "ehubDLR_solar_$year.nc")
    request_era5_vars(outfile, solarvars, date1, date2; res=0.25, bbox=(3.5, 33.0, 53.5, 72.5))
end

function read_weatherdata_DLR(year)
    nc_wt = Dataset(in_datafolder("DLR", "ehubDLR_windtemp_$year.nc"))
    nc_s = Dataset(in_datafolder("DLR", "ehubDLR_solar_$year.nc"))
    sz = size(nc_wt["u100"])

    res = 0.25
    bbox = (3.5, 33.0, 53.5, 72.5)

    geo = GeoArray(zeros(sz[1:2]), res, bbox)
    u100 = permutedims(nc_wt["u100"][:,:,:], [3,1,2])   # instantaneous, eastward [m/s]
    v100 = permutedims(nc_wt["v100"][:,:,:], [3,1,2])   # instantaneous, northward [m/s]
    t2m = permutedims(nc_wt["t2m"][:,:,:], [3,1,2])     # instantaneous [K]
    ssrd = permutedims(nc_s["ssrd"][:,:,:], [3,1,2])    # accumulated [J/m2/period]
    fdir = permutedims(nc_s["fdir"][:,:,:], [3,1,2])    # accumulated [J/m2/period]

    lons = bbox[1] + res/2 : res : bbox[2] - res/2
    lats = bbox[4] - res/2 : -res : bbox[3] + res/2     # lats in descending order
    return (; geo, u100, v100, t2m, ssrd, fdir, lons, lats, res, year)
end

function get_cell_weather!(cell_weather, cell, line_azimuth, line_diameter, weatherdata)
    lon, lat = cell
    (; geo, u100, v100, t2m, ssrd, fdir, year) = weatherdata   # u is eastward, v is northward
    (; temp_air, wind_speed, wind_angle, insolation, wind_u, wind_v, SSRD, FDIR) = cell_weather

    # Note that while the solar position calculations are instantaneous positions, ERA5 radiation variables
    # represent accumulated radiation over the hour *ending* at the indicated time. Therefore, solar positions
    # must be shifted 30 minutes BACK to correspond to the midpoint time of the ERA5 accumulations.
    # Source: ERA5 "accumulations are over the hour ending at the forecast step"
    # https://confluence.ecmwf.int//display/CKB/ERA5+data+documentation#ERA5datadocumentation-Meanratesandaccumulations
    datetime = DateTime(year,1,1) - Minute(30):Hour(1):DateTime(year,12,31,23) - Minute(30)
    time = 1:8760
    index = lonlat_index(geo, lon, lat)
    temp_air .= t2m[time, index] .- 273.15   # [°C]
    wind_u .= u100[time, index]
    wind_v .= v100[time, index]
    wind_speed .= sqrt.(wind_u.^2 + wind_v.^2)           # [m/s]
    wind_angle .= mod.(atand.(wind_v, wind_u), 360)      # angle from North, clockwise

    almostzero = eps(Float32)

    # ERA5 radiations are in J/m2/period, so for hourly data divide by 3600 to get W/m2
    SSRD .= ssrd[time, index] ./ 3600  # total (global) horizontal insolation [W/m2]
    # FDIR .= fdir[time, index] ./ 3600  # direct insolation on a horizontal surface [W/m2]
    insolation .= 0.0
    for (i, dt) in enumerate(datetime)
        TSI = solarinsolation(dt)                       # Total Solar Irradiance (top of atmosphere, perpendicular to sun) [W/m2]
        δ, H = solarposition(dt, lon)                   # absolute solar position (declination, hour angle)
        solarpos = sines_and_cosines(δ, H)
        zen, az = zenith_azimuth(lat, solarpos...)      # relative solar position (radians)
        cos_zen = max(almostzero, cos(zen))
        
        # When the solar elevation is close to 0, both FDIR and cos(zenith) will also be near 0, and
        # calculated DNI will approach "0/0". So we'll clamp DNI to avoid artifacts.
        # That wasn't enough, so we'll add an artificial term to increase the denominator near the horizon.
        # Also, we'll use SSRD instead of FDIR to capture total insolation, so we have GNI instead of DNI.
        GNI = clamp(SSRD[i] / (cos_zen + horizoncorrection(zen)), 0, TSI)  # Global Normal Irradiance [W/m2]

        zenith, azimuth = rad2deg(zen), rad2deg(az)
        zenith > 90 && continue                         # sun below horizon

        Δaz = azimuth - line_azimuth
        cosθ = sind(zenith) * cosd(Δaz)
        sinθ = sqrt(1 - cosθ^2)
        insolation[i] = GNI * line_diameter*1e-3 * sinθ   # direct component only [W/m] (per unit length of line)
    end
end

function hourly_weather(hour, cell_weather)
    return (;
        temp_air = cell_weather.temp_air[hour],
        wind_speed = cell_weather.wind_speed[hour],
        wind_angle = cell_weather.wind_angle[hour],
        insolation = cell_weather.insolation[hour]
    )
end

"Collect all the line and conductor data we need from Excel files into a DataFrame."
function read_line_data()
    xlsx_buslines = in_datafolder("DLR", "Bus_and_line_data_EHUB400_future_data_v1_11_2.xlsx")
    lines = XLSX.readtable(xlsx_buslines, "lines"; infer_eltypes=true) |> DataFrame
    buses = XLSX.readtable(xlsx_buslines, "buses"; infer_eltypes=true) |> DataFrame
    rename!(lines, ["line_id", "start_node", "end_node", "resistance", "reactance", "voltage", "transformer",
                    "c_rating", "length", "geometry", "conductor_count", "conductor_type"])
    select!(buses, ["bus_id", "x-coordinate", "y-coordinate"])
    rename!(buses, "x-coordinate"=>"lon", "y-coordinate"=>"lat")

    xlsx_conductors = in_datafolder("DLR", "Conductor_data.xlsx")
    df_cond = XLSX.readtable(xlsx_conductors, "Line_data"; infer_eltypes=true) |> DataFrame
    select!(df_cond, ["Name", "Diam", "DC-res"])
    rename!(df_cond, ["name", "diameter", "dc_resistance"])
    df_cond.name .= uppercase.(strip.(df_cond.name))

    leftjoin!(lines, df_cond, on=:conductor_type=>:name, matchmissing=:notequal)

    leftjoin!(lines, buses, on=:start_node=>:bus_id)
    rename!(lines, :lon=>:start_lon, :lat=>:start_lat)
    leftjoin!(lines, buses, on=:end_node=>:bus_id)
    rename!(lines, :lon=>:end_lon, :lat=>:end_lat)

    distance = zeros(nrow(lines)) 
    for (i, row) in enumerate(eachrow(lines))
        p1 = (row.start_lon, row.start_lat)
        p2 = (row.end_lon, row.end_lat)
       distance[i] = greatcircledistance(p1, p2)
    end
    insertcols!(lines, 10, :distance => distance)

    # lines = lines[.!lines.transformer, :]
    # select!(lines, Not(:transformer, :geometry))

    lines.conductor_count .= [ismissing(c) ? 1 : c == "duplex" ? 2 : c == "triplex" ? 3 : 1 for c in lines.conductor_count]

    # lines.length .= coalesce.(lines.length)
    # lines.diameter .= coalesce.(lines.diameter)
    # lines.dc_resistance .= coalesce.(lines.dc_resistance)
    # lines.conductor_type .= coalesce.(lines.conductor_type)

    return lines
end

function write_csv(basename, var, columns)
    df = DataFrame(var, string.(columns))
    CSV.write(in_datafolder("DLR", "$basename.csv"), df)
end

function read_csv(basename)
    df = CSV.read(in_datafolder("DLR", "$basename.csv"), DataFrame)
    return df
end
