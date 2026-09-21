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
lon-lat grid with specified resolution (in degrees). Returns path length and bearings within each cell,
ordered from point1 to point2.

Crossings of the path with grid meridians and parallels are calculated analytically, so no cells are
missed even when the path only clips the corner of a cell. A path that passes (numerically) exactly
through a grid corner continues directly into the diagonally opposite cell. Assumes that the path
doesn't cross the antimeridian or pass over a pole.
"""
function greatcircle_waypoints(point1::Tuple, point2::Tuple, grid_resolution::Float64)
    waypoints = @NamedTuple{cell::Tuple{Float64, Float64}, len::Float64, mean_bearing::Float64, bearing_error::Float64}[]
    angular_dist = calc_angular_distance(point1, point2)
    angular_dist == 0 && return waypoints

    # Parametrize the path as p(t) = cosd(t)*a + sind(t)*b for t in [0, angular_dist], with a, b orthonormal.
    a, p2 = lonlat2xyz(point1), lonlat2xyz(point2)
    b = normalize3(p2 .- dot3(a, p2) .* a)
    pathpoint(t) = xyz2lonlat(cosd(t) .* a .+ sind(t) .* b)

    crossings = gridline_crossings(point1, point2, a, b, angular_dist, grid_resolution)
    ts = sort!([0.0; crossings; angular_dist])

    # Each segment between consecutive crossings lies within a single cell, so look up the cell of its midpoint.
    min_segment = 1e-9      # [degrees] skip zero length segments, e.g. when the path hits a grid corner exactly
    cells, t_entry, t_exit = Tuple{Float64, Float64}[], Float64[], Float64[]
    for i in 1:length(ts)-1
        t1, t2 = ts[i], ts[i+1]
        t2 - t1 < min_segment && continue
        lon, lat = pathpoint((t1 + t2) / 2)
        cell = (round_res(lon, grid_resolution), round_res(lat, grid_resolution))
        if !isempty(cells) && cell == cells[end]
            t_exit[end] = t2    # still in the same cell (e.g. after skipping a zero length segment)
        else
            push!(cells, cell); push!(t_entry, t1); push!(t_exit, t2)
        end
    end

    for (cell, t1, t2) in zip(cells, t_entry, t_exit)
        entry_point, exit_point = pathpoint(t1), pathpoint(t2)

        len = greatcircledistance(entry_point, exit_point)
        bearings = greatcirclebearings(entry_point, exit_point)

        mean_bearing = mean(bearings)
        bearing_error = maximum(abs.(bearings .- mean_bearing))

        push!(waypoints, (; cell, len, mean_bearing, bearing_error))
    end

    if !all_cells_adjacent(waypoints, grid_resolution)
        @warn "Consecutive cells along line path are not adjacent (this shouldn't happen)." point1 point2
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

"Convert (lon, lat) in degrees to a unit vector (x, y, z)."
lonlat2xyz((lon, lat)) = (cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))

"Convert a vector (x, y, z) to (lon, lat) in degrees."
xyz2lonlat((x, y, z)) = (atand(y, x), atand(z, hypot(x, y)))

dot3(u, v) = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
normalize3(u) = u ./ sqrt(dot3(u, u))

"""
Find the parameters t (degrees, 0 < t < tmax) where the great circle arc p(t) = cosd(t)*a + sind(t)*b
crosses a grid meridian or parallel (at multiples of grid_resolution).
"""
function gridline_crossings(point1, point2, a, b, tmax, grid_resolution)
    (lon1, lat1), (lon2, lat2) = point1, point2
    res = grid_resolution
    ts = Float64[]

    # Longitude changes monotonically along the arc, so each grid meridian strictly between the endpoints
    # is crossed exactly once. The meridian at longitude λ lies in the plane with normal (-sind(λ), cosd(λ), 0).
    lonmin, lonmax = minmax(lon1, lon2)
    for k in (floor(Int, lonmin / res) + 1):(ceil(Int, lonmax / res) - 1)
        λ = k * res
        append!(ts, arc_crossings(a[2]*cosd(λ) - a[1]*sind(λ), b[2]*cosd(λ) - b[1]*sind(λ), 0.0, tmax))
    end

    # The parallel at latitude φ is where z = sind(φ). Latitude isn't monotonic along the arc (it can bulge
    # poleward of both endpoints), so first find its latitude range. Along the arc z(t) = R*cosd(t - δ).
    R, δ = hypot(a[3], b[3]), atand(b[3], a[3])
    latmin, latmax = minmax(lat1, lat2)
    mod(δ, 360) < tmax && (latmax = asind(min(R, 1.0)))          # northernmost point of great circle is on the arc
    mod(δ + 180, 360) < tmax && (latmin = -asind(min(R, 1.0)))   # southernmost point of great circle is on the arc
    for k in ceil(Int, latmin / res):floor(Int, latmax / res)
        append!(ts, arc_crossings(a[3], b[3], sind(k * res), tmax))
    end

    return ts
end

"Solve α*cosd(t) + β*sind(t) = γ for t (degrees) in the open interval (0, tmax)."
function arc_crossings(α, β, γ, tmax)
    R = hypot(α, β)
    abs(γ) >= R && return Float64[]     # no crossing (or just touching the grid line)
    δ, Δ = atand(β, α), acosd(γ / R)
    return filter(t -> 0 < t < tmax, unique(mod.([δ - Δ, δ + Δ], 360)))
end

# Round a value to the nearest multiple of the given resolution, offset by resolution/2, treating -0.0 as 0.0.
function round_res(value, resolution)
    res2 = resolution / 2
    rounded = round((value - res2) / resolution) * resolution + res2
    return rounded == -0.0 ? 0.0 : rounded
end

# Cells sharing an edge or a corner (a path through a grid corner goes directly to the diagonal cell).
adjacentcells(cell1, cell2, res) = cell1 != cell2 && all(abs.(cell1 .- cell2) .< 1.5 * res)
all_cells_adjacent(waypoints, res) = all(adjacentcells(waypoints[i].cell, waypoints[i+1].cell, res) for i = 1:length(waypoints)-1)

function download_era5_DLR(year)
    date1, date2 = "$year-01-01", "$year-12-31"

    # Split into two requests, one for instantaneous variables (wind/temp)
    # and one for accumulated (solar) (otherwise Copernicus returns a zip file)
    windvars = ["100m_u_component_of_wind", "100m_v_component_of_wind", "10m_u_component_of_wind", "10m_v_component_of_wind", "2m_temperature"]
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
    u10 = permutedims(nc_wt["u10"][:,:,:], [3,1,2])     # instantaneous, eastward [m/s]
    v10 = permutedims(nc_wt["v10"][:,:,:], [3,1,2])     # instantaneous, northward [m/s]
    u100 = permutedims(nc_wt["u100"][:,:,:], [3,1,2])   # instantaneous, eastward [m/s]
    v100 = permutedims(nc_wt["v100"][:,:,:], [3,1,2])   # instantaneous, northward [m/s]
    t2m = permutedims(nc_wt["t2m"][:,:,:], [3,1,2])     # instantaneous [K]
    ssrd = permutedims(nc_s["ssrd"][:,:,:], [3,1,2])    # accumulated [J/m2/period]
    fdir = permutedims(nc_s["fdir"][:,:,:], [3,1,2])    # accumulated [J/m2/period]

    lons = bbox[1] + res/2 : res : bbox[2] - res/2
    lats = bbox[4] - res/2 : -res : bbox[3] + res/2     # lats in descending order
    return (; geo, u10, v10, u100, v100, t2m, ssrd, fdir, lons, lats, res, year)
end

"Interpolate wind between 10 m and 100 m: speed by power law (hourly shear), direction along the shortest arc."
function interpolate_wind(u10, v10, u100, v100, height)
    w = log(height / 10) / log(10)
    speed = hypot(u10, v10)^(1 - w) * hypot(u100, v100)^w
    dir10, dir100 = atand(v10, u10), atand(v100, u100)
    dir = dir10 + w * (mod(dir100 - dir10 + 180, 360) - 180)
    return speed * cosd(dir), speed * sind(dir)
end

function get_cell_weather!(cell_weather, cell, line_azimuth, line_height, weatherdata)
    lon, lat = cell
    (; geo, u10, v10, u100, v100, t2m, ssrd, fdir, year) = weatherdata   # u is eastward, v is northward
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
    for i in time
        wind_u[i], wind_v[i] = interpolate_wind(u10[i, index], v10[i, index], u100[i, index], v100[i, index], line_height)
    end
    wind_speed .= sqrt.(wind_u.^2 + wind_v.^2)           # [m/s]
    wind_angle .= mod.(atand.(wind_u, wind_v), 360)      # direction wind blows toward, clockwise from North

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
        insolation[i] = GNI * sinθ      # effective insolation on the conductor [W/m2] (IEEE 738: qs = α * Q * sinθ * D)
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
    xlsx_buslines = in_datafolder("DLR", "Bus_and_line_data_EHUB400_future_data_v1_12.xlsx")
    lines = XLSX.readtable(xlsx_buslines, "lines"; infer_eltypes=true) |> DataFrame
    buses = XLSX.readtable(xlsx_buslines, "buses"; infer_eltypes=true) |> DataFrame
    rename!(lines, ["line_id", "start_node", "end_node", "resistance", "reactance", "voltage", "transformer",
                    "c_rating", "length", "geometry", "conductor_count", "conductor_type", "circuits", "b_f", "b_t", "series_compensation", "country"])
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
