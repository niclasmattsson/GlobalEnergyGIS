using DataFrames, XLSX

"Calculate static and dynamic line ratings according to the IEEE 738 standard for overhead conductors."
function calculate_line_ratings(lines::DataFrame, weatherdata::NamedTuple)
    # units of lines columns: diameter [mm], dc_resistance [ohm/km], c_rating [A], voltage [kV], reactance [p.u.]
    # weatherdata fields: (; u100, v100, t2m, ssrd, fdir, lons, lats, res)
    line_params = (
        temp_line = 50.0,       # Max conductor surface temperature [°C]
        elevation = 0.0,        # Elevation above sea level [m]
        emissivity = 0.9,
        absorptivity = 0.5
    )
    (; temp_line, elevation, emissivity, absorptivity) = line_params

    nhours, nlines = 8760, nrow(lines)
    mean_bearings = zeros(nlines)
    min_line_ampacity = zeros(nhours)
    segment_ampacities = zeros(nhours)
    ampacity = zeros(nhours, nlines)                # max current carrying capacity [A]
    thermal_capacity = zeros(nhours, nlines)        # dynamic thermal capacity [MW]

    cell_weather = (;   # initialize cell weather vectors to be filled in get_cell_weather! (to avoid repeated allocations)
        temp_air=zeros(nhours), wind_speed=zeros(nhours), wind_angle=zeros(nhours), insolation=zeros(nhours),
        wind_u=zeros(nhours), wind_v=zeros(nhours), SSRD=zeros(nhours), FDIR=zeros(nhours)
    )
    mean_line_weather = (;  # to accumulate mean weather along the line
        temp_air=zeros(nhours, nlines), wind_speed=zeros(nhours, nlines),
        wind_angle=zeros(nhours, nlines), insolation=zeros(nhours, nlines)
    )
    dimensioning_line_weather = deepcopy(mean_line_weather)     # to store weather of the dimensioning segment for each hour

    updateprogress = Progress(nlines, 1)
    for (i, line) in enumerate(eachrow(lines))
        next!(updateprogress)
        line.transformer && continue   # skip transformers

        linestart, lineend = getlinecoords(line)
        linesegments = greatcircle_waypoints(linestart, lineend, weatherdata.res)

        min_line_ampacity .= Inf
        segment_ampacities .= 0.0

        for segment in linesegments
            (; cell, mean_bearing) = segment    # NM: don't we need len anywhere???
            get_cell_weather!(cell_weather, cell, mean_bearing, line.diameter, weatherdata)
            mean_bearings[i] += mean_bearing

            Threads.@threads for hour in 1:nhours
                weather = hourly_weather(hour, cell_weather)
                segment_ampacities[hour] = calculate_ampacity(line, mean_bearing, weather, line_params)  # [A]
            end

            new_minimum = segment_ampacities .< min_line_ampacity   # hours where this segment is dimensioning
            min_line_ampacity[new_minimum] .= segment_ampacities[new_minimum]
            for k in keys(dimensioning_line_weather)
                dimensioning_line_weather[k][new_minimum, i] .= cell_weather[k][new_minimum]
            end

            for k in keys(mean_line_weather)
                mean_line_weather[k][:, i] .+= cell_weather[k]      # accumulate weather along the line
            end
        end

        mean_bearings[i] /= length(linesegments)
        for k in keys(mean_line_weather)
            mean_line_weather[k][:, i] ./= length(linesegments)     # average hourly weather along the line
        end

        ampacity[:, i] .= line.circuits * min_line_ampacity                             # [A]
        thermal_capacity[:, i] .= thermal_capacity_limit(line.voltage, ampacity[:, i])  # [MW]
    end

    calculate_static_line_ratings!(lines; max_power_angle=30)       # Static line rating capacities (MW)
    lines.mean_bearing .= mean_bearings
    CSV.write(in_datafolder("DLR", "line_data.csv"), lines)

    thermal_ratio = calculate_thermal_ratio(thermal_capacity, lines)
    thermal_ratio_noangle = thermal_capacity ./ lines.SLR_thermal'  # ratio of IEEE dynamic to static thermal max
    
    calculate_static_line_ratings!(lines; max_power_angle=20)
    thermal_ratio_20 = calculate_thermal_ratio(thermal_capacity, lines)
    calculate_static_line_ratings!(lines; max_power_angle=40)
    thermal_ratio_40 = calculate_thermal_ratio(thermal_capacity, lines)
    calculate_static_line_ratings!(lines; max_power_angle=50)
    thermal_ratio_50 = calculate_thermal_ratio(thermal_capacity, lines)

    line_ids = lines.line_id
    write_csv("ampacity", ampacity, line_ids)
    write_csv("thermal_ratio", thermal_ratio, line_ids)
    write_csv("thermal_ratio_noangle", thermal_ratio_noangle, line_ids)
    write_csv("thermal_ratio_20", thermal_ratio_20, line_ids)
    write_csv("thermal_ratio_40", thermal_ratio_40, line_ids)
    write_csv("thermal_ratio_50", thermal_ratio_50, line_ids)
    write_csv("mean_temp_air", mean_line_weather.temp_air, line_ids)
    write_csv("mean_wind_speed", mean_line_weather.wind_speed, line_ids)
    write_csv("mean_wind_angle", mean_line_weather.wind_angle, line_ids)
    write_csv("mean_insolation", mean_line_weather.insolation, line_ids)
    write_csv("dimensioning_temp_air", dimensioning_line_weather.temp_air, line_ids)
    write_csv("dimensioning_wind_speed", dimensioning_line_weather.wind_speed, line_ids)
    write_csv("dimensioning_wind_angle", dimensioning_line_weather.wind_angle, line_ids)
    write_csv("dimensioning_insolation", dimensioning_line_weather.insolation, line_ids)

    return ampacity
end

"Calculate the total current carrying capacity of a line."
function calculate_ampacity(line, mean_bearing, cell_hourly_weather, line_params)
    (; line_id, conductor_count, conductor_type, diameter, dc_resistance) = line
    (; temp_air, wind_speed, wind_angle, insolation) = cell_hourly_weather
    (; temp_line, elevation, emissivity, absorptivity) = line_params

    temp_film = (temp_line + temp_air) / 2          # i.e. the thermal boundary layer around the conductor
    wind_attack_angle = abs(mod(mean_bearing - wind_angle + 90, 180) - 90)

    k_f = air_conductivity(temp_film)           # [W/m·K]
    rho_f = air_density(temp_film, elevation)   # [kg/m^3]
    mu_f = air_viscosity(temp_film)             # [Pa·s]
    k_angle = wind_cooling_factor(wind_attack_angle)

    diam_m = diameter * 1e-3                                # Outer diameter [m]
    N_re = reynolds_number(diam_m, rho_f, wind_speed, mu_f)

    # Heat terms, all three in [W/m] (i.e per unit length of the line)
    qc = convective_cooling(k_angle, N_re, k_f, temp_line, temp_air, rho_f, diam_m)
    qr = radiative_cooling(temp_line, temp_air, diam_m, emissivity)
    qs = solar_heat_absorption(absorptivity, insolation, diam_m)

    resistance = temp_adjusted_resistance(dc_resistance, temp_line)        # electrical resistance at temp_line [ohm/m]

    I_per = max_current(qc, qr, qs, resistance)
    I_total = I_per * conductor_count                   # total current for all conductors [A]

    return I_total
end

# -------------------------------------------------------
# Fluid properties vs film temperature (°C)
# -------------------------------------------------------

"Thermal conductivity of air [W/m·K]"
function air_conductivity(temp_film)
    return 2.424e-2 + 7.477e-5 * temp_film - 4.407e-9 * (temp_film^2)
end

"Air density [kg/m^3] with elevation correction"
function air_density(temp_film, elevation)
    numerator = 1.293 - 1.525e-4 * elevation + 6.379e-9 * (elevation^2)
    denominator = 1 + 0.00367 * temp_film
    return numerator / denominator
end

"Dynamic viscosity of air [Pa·s]"
function air_viscosity(temp_film)
    return (1.458e-6 * ((temp_film + 273.0)^1.5)) / (temp_film + 383.4)
end

# -------------------------------------------------------
# Heat transfer pieces (IEEE-style)
# -------------------------------------------------------

"Reynolds number (assumes diam_m in [m])"
function reynolds_number(diam_m, rho_f, wind_speed, mu_f)
    return (diam_m * rho_f * wind_speed) / mu_f
end

"Wind angle correction factor (k_angle), IEEE 738 empirical model of cooling effectiveness."
function wind_cooling_factor(phi)
    return 1.194 - cosd(phi) + 0.194 * cosd(2 * phi) + 0.368 * sind(2 * phi)
end

"Convective cooling per unit length [W/m]. Three correlations; pick the maximum."
function convective_cooling(k_angle, N_re, k_f, temp_line, temp_air, rho_f, diam_m)
    deltaT = temp_line - temp_air

    # Note: k_f already in W/mK; diam_m implicit in correlations per user's math.
    qc1 = k_angle * (1.01 + 1.35 * (N_re^0.52)) * k_f * deltaT
    qc2 = k_angle * 0.754 * (N_re^0.6) * k_f * deltaT
    qc3 = 3.645 * (rho_f^0.5) * (diam_m^0.75) * (abs(deltaT)^1.25) * (deltaT >= 0 ? 1 : -1)

    return max(qc1, qc2, qc3)
end

"Radiative cooling per unit length [W/m]."
function radiative_cooling(temp_line, temp_air, diam_m, emissivity)
    eps = max(0.0, min(1.0, float(emissivity)))
    term_s = ((temp_line + 273.0) / 100.0)^4
    term_a = ((temp_air + 273.0) / 100.0)^4
    return 17.8 * diam_m * eps * (term_s - term_a)
end

"Solar heat gain per unit length [W/m]."
function solar_heat_absorption(absorptivity, insolation, diam_m)
    return absorptivity * insolation * diam_m
end

"DC resistance at temperature T_c [ohm/m]. Uses fixed temp coeff 0.0039 per the provided IEEE math block."
function temp_adjusted_resistance(line_dc_resistance_ohm_per_km, T_c)
    return line_dc_resistance_ohm_per_km * 0.001 * (1 + 0.0039 * (T_c - 20.0))
end

"Maximum current per conductor [A] from heat balance."
function max_current(qc, qr, qs, R_ohm_per_m)
    numerator = qc + qr - qs
    if R_ohm_per_m <= 0 || numerator < 0
        return 1000.0
        # error("Invalid parameters for calculate_imax")
    end
    return sqrt(numerator / R_ohm_per_m)
end

function calculate_thermal_ratio(thermal_capacity, lines::DataFrame)
    max_capacity = min.(thermal_capacity, lines.SLR_angle')         # [MW]
    thermal_ratio = max_capacity ./ lines.SLR_max'                  # ratio of IEEE dynamic max to static max
    return thermal_ratio
end

"Calculate static line rating capacities (all in MW)."
function calculate_static_line_ratings!(lines::DataFrame; max_power_angle)
    lines.SLR_thermal .= thermal_capacity_limit.(lines.voltage, lines.c_rating)
    lines.SLR_angle .= angle_capacity_limit.(lines.reactance, max_power_angle)
    lines.SLR_max .= min.(lines.SLR_thermal, lines.SLR_angle)
    return nothing
end

"Calculate thermal line capacities."
function thermal_capacity_limit(voltage, current)
    return sqrt(3) * voltage * 1000 .* current / 1e6  # [MW]  (voltages [kV], currents [A])
end

"Calculate power angle limits."
function angle_capacity_limit(reactance, max_power_angle)
    return (1 ./ reactance) * max_power_angle / 57.29 * 1000   # [MW]
end
