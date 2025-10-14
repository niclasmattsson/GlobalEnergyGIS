using DataFrames, XLSX

"Calculate static and dynamic line ratings according to the IEEE 738 standard for overhead conductors."
function calculate_line_ratings(lines::DataFrame, weatherdata::NamedTuple)
    # units of lines columns: diameter [mm], dc_resistance [ohm/km], c_rating [A], voltage [kV], reactance [p.u.]
    # weatherdata fields: (; u100, v100, t2m, ssrd, fdir, lons, lats, res)
    line_params = (
        temp_line = 50.0,       # Max conductor surface temperature [°C]
        max_power_angle = 30,
        elevation = 0.0,        # Elevation above sea level [m]
        emissivity = 0.9,
        absorptivity = 0.5
    )
    (; temp_line, max_power_angle, elevation, emissivity, absorptivity) = line_params

    nhours, max_segments, nlines = 8760, 50, nrow(lines)
    maxcurrent_segment = zeros(nhours, max_segments)
    ampacity = zeros(nhours, nlines)            # max current carrying capacity [A]
    thermal_capacity = zeros(nhours, nlines)    # dynamic thermal capacity [MW]

    cell_weather = (;   # initialize cell weather vectors to be filled in get_cell_weather! (to avoid repeated allocations)
        temp_air=zeros(nhours), wind_speed=zeros(nhours), wind_angle=zeros(nhours), insolation=zeros(nhours),
        wind_u=zeros(nhours), wind_v=zeros(nhours), SSRD=zeros(nhours), FDIR=zeros(nhours)
    )

    updateprogress = Progress(nlines, 1)
    for (i, line) in enumerate(eachrow(lines))
        next!(updateprogress)
        line.transformer && continue   # skip transformers

        linestart, lineend = getlinecoords(line)
        linesegments = greatcircle_waypoints(linestart, lineend, weatherdata.res)

        maxcurrent_segment .= 0.0
        for (seg, segment) in enumerate(linesegments)
            (; cell, len, mean_bearing) = segment   # NM: don't we need len anywhere???
            get_cell_weather!(cell_weather, cell, mean_bearing, line.diameter, weatherdata)
            Threads.@threads for hour in 1:nhours
                weather = hourly_weather(hour, cell_weather)
                maxcurrent_segment[hour, seg] = calculate_ampacity(line, weather, line_params)      # [A]
            end
        end

        ampacity[:, i] .= minimum((@view maxcurrent_segment[:, 1:length(linesegments)]), dims=2)    # [A]
        thermal_capacity[:, i] .= thermal_capacity_limit(line.voltage, ampacity[:, i])              # [MW]
    end

    # Static line rating capacities (MW)
    lines.SLR_thermal .= thermal_capacity_limit.(lines.voltage, lines.c_rating)
    lines.SLR_angle .= angle_capacity_limit.(lines.reactance, max_power_angle)
    lines.SLR_max .= min.(lines.SLR_thermal, lines.SLR_angle)

    CSV.write(in_datafolder("downloads", "Processed_line_data.csv"), lines)

    max_capacity = min.(thermal_capacity, lines.SLR_angle')     # [MW]
    thermal_ratio = max_capacity ./ lines.SLR_max'              # ratio of IEEE dynamic max to static max

    df_ampacity = DataFrame(ampacity, string.(lines.line_id))
    df_thermal_ratio = DataFrame(thermal_ratio, string.(lines.line_id))
    CSV.write(in_datafolder("downloads", "ampacity_test.csv"), df_ampacity)
    CSV.write(in_datafolder("downloads", "thermal_ratio_test.csv"), df_thermal_ratio)
    return ampacity

    # # If Voltage == 0 or NaN → use angle limit only
    # # mask_V0 = ismissing.(df.voltage) .|| df.voltage .== 0
    # # df[mask_V0, :max_capacity] .= angle_capacity_limits(df[mask_V0, :reactance], max_power_angle)
    # # df[mask_V0, :SLR_max] .= angle_capacity_limits(df[mask_V0, :reactance], max_power_angle)
end

"Calculate the total current carrying capacity of a line."
function calculate_ampacity(line, cell_weather, line_params)
    (; line_id, conductor_count, conductor_type, diameter, dc_resistance) = line
    (; temp_air, wind_speed, wind_angle, insolation) = cell_weather
    (; temp_line, elevation, emissivity, absorptivity) = line_params

    temp_film = (temp_line + temp_air) / 2          # i.e. the thermal boundary layer around the conductor

    k_f = air_conductivity(temp_film)           # [W/m·K]
    rho_f = air_density(temp_film, elevation)   # [kg/m^3]
    mu_f = air_viscosity(temp_film)             # [Pa·s]
    k_angle = wind_cooling_factor(wind_angle)

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

"Calculate thermal line capacities."
function thermal_capacity_limit(voltage, current)
    return sqrt(3) * voltage * 1000 .* current / 1e6  # [MW]  (voltages [kV], currents [A])
end

"Calculate power angle limits."
function angle_capacity_limit(reactance, max_power_angle)
    return (1 ./ reactance) * max_power_angle / 57.29 * 1000   # [MW]
end
