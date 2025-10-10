# Unfinished. Needs renaming of identifiers for consistency with other code, complete thermal ratio calculation, etc.
# check inputs to calculate_thermal_ratio (cf python code), add more lines data to dataframe

"Calculates the dynamic thermal capacity ratio for a single time step."
function calculate_thermal_ratio(
        Ta, solar_radiation, wind_speed, wind_angle,    # Weather parameters
        conductor_params,                               # Conductor properties
        line_params                                     # Line parameters
    )

    (; material, nominal_area, diameter, outer_strand_diameter, num_layers, dc_resistance, resistance_coeff) = conductor_params
    (; voltage_kV, num_conductors) = line_params

    defaults = (
        emissivity = 0.5,
        absorptivity = 0.5,
        Tc = 80.0,            # Max conductor temperature [°C]
        elevation = 0.0,      # Elevation [m]
    )
    (; Tc, elevation, emissivity, absorptivity) = defaults

    t_film = (Tc + Ta) / 2.0

    Pc = convection_power(Tc, Ta, t_film, elevation, wind_speed, wind_angle, diameter, outer_strand_diameter)
    Pr = radiation_power(diameter, Tc, Ta, emissivity)
    Ps = solar_gain(diameter, absorptivity, solar_radiation)

    R_dc = resistance_dc(dc_resistance, resistance_coeff, Tc)
    
    numerator = Pc + Pr - Ps
    if numerator <= 0 || R_dc == 0
        return 0.0      # No thermal capacity if heat loss is negative or resistance is zero
    end

    current_dc = sqrt(numerator / R_dc)
    Ik = current_dc / nominal_area
    strom_corr = ac_correct_current(current_dc, Ik, nominal_area, num_layers, material)
    strom_total = strom_corr * num_conductors   # Total current for the line (e.g., duplex/triplex)

    thermal_capacity_MVA = sqrt(3) * voltage_kV * strom_total / 1000.0  # Thermal capacity [MVA]
    
    return thermal_capacity_MVA
end

"Calculates the convective heat loss (power) from the conductor."
function convection_power(
        Tc,                         # Conductor temperature [°C]
        Ta,                         # Ambient temperature [°C]
        t_film,                     # Film temperature [°C], typically (Tc + Ta) / 2
        elevation,                  # Elevation [m]
        wind_speed,                 # Wind speed [m/s]
        wind_angle_deg,             # Wind attack angle [°]
        conductor_diameter,         # Conductor diameter [mm]
        outer_strand_diameter       # Diameter of outer strands [mm]
    )

    ypsilon = 0.0000132 + 0.000000095 * t_film  # Dynamic viscosity of air [Pa·s]
    lambda_f = 0.0242 + 0.000072 * t_film       # Thermal conductivity of air [W/(m·K)]
    rho_r = exp(-0.000116 * elevation)          # Air density correction factor for elevation

    D_m = conductor_diameter * 0.001         # Conductor diameter [m]
    Re = (rho_r * wind_speed * D_m) / ypsilon   # Reynolds number

    # Stranded conductor roughness factor
    R_f = outer_strand_diameter / (2 * (conductor_diameter - 2 * outer_strand_diameter))

    B1, n = if Re <= 2650
        0.641, 0.471
    elseif R_f > 0.05
        0.048, 0.8
    else
        0.178, 0.633
    end

    Nu_90 = B1 * Re^n     # Nusselt number for wind angle of 90°

    A1, B2, M1 = if wind_angle_deg > 24
        0.42, 0.58, 0.9
    else
        0.42, 0.68, 1.08
    end

    Nu = if wind_speed < 0.5
        # Natural convection case (low wind)
        Gr = (D_m^3 * (Tc - Ta) * 9.807) / ((t_film + 273) * ypsilon^2)     # Grashof number
        Pr_2 = 0.715 - 0.00025 * t_film                                     # Prandtl number
        
        A2, m2 = (Gr * Pr_2 < 1e4) ? (0.85, 0.188) : (0.48, 0.25)

        Nu2 = A2 * (Gr * Pr_2)^m2
        Nu_45 = Nu_90 * (0.42 + 0.58 * (sind(45))^0.9)
        Nucorr = 0.55 * Nu_90
        max(Nu2, Nu_45, Nucorr)
    else
        # Forced convection case
        Nu_90 * (A1 + B2 * sind(wind_angle_deg)^M1)
    end

    Pc = pi * lambda_f * (Tc - Ta) * Nu         # Convective heat loss [W/m]
    return Pc
end

"Calculates the radiative heat loss from the conductor."
function radiation_power(
        conductor_diameter,     # Conductor diameter [mm]
        Tc,                     # Conductor temperature [°C]
        Ta,                     # Ambient temperature [°C]
        emissivity              # Emissivity of the conductor surface
    )
    D_m = conductor_diameter * 0.001    # Conductor diameter [m]
    σ = 5.67e-8                         # Stefan–Boltzmann constant [W·m⁻²·K⁻⁴]
    return pi * D_m * emissivity * σ * ((Tc + 273.15)^4 - (Ta + 273.15)^4)  # Radiative heat loss [W/m]
end

"Calculates the heat gain from solar radiation."
function solar_gain(
        conductor_diameter,     # Conductor diameter [mm]
        absorptivity,           # Absorptivity of the conductor surface
        solar_radiation         # Global solar radiation [W/m²]
    )
    return absorptivity * conductor_diameter * solar_radiation * 0.001  # Solar heat gain [W/m]
end

"Calculates the DC resistance of the conductor at a given temperature."
function resistance_dc(
        dc_resistance_at_20C,   # DC resistance at 20°C [Ω/km]
        resistance_temp_coeff,  # Resistance temperature coefficient [1/K]
        Tc                      # Conductor temperature [°C]
    )
    return dc_resistance_at_20C * 0.001 * (1 + resistance_temp_coeff * (Tc - 20.0)) # DC resistance [Ω/m]
end

"Corrects the DC current capacity to AC, accounting for skin and proximity effects."
function ac_correct_current(
        current_dc,            # DC current [A]
        Ik,                    # Current density [A/mm²]
        nominal_area,          # Conductor nominal area [mm²]
        num_layers,            # Number of layers in the conductor
        material_str           # Conductor material ("ACSR", "FE", etc.)
    )
    
    mat = uppercase(strip(material_str))
    if (mat == "ACSR") || (mat == "FE")
        if num_layers == 3
            return current_dc / sqrt(1.0123 + 0.00002319 * current_dc)
        elseif nominal_area >= 175
            return current_dc / sqrt(1.0045 + 0.00000009 * current_dc)
        elseif Ik > 3.908
            return current_dc / sqrt(1.1)
        elseif Ik > 2.486
            denom = 1 + 0.02 * (2.978 - 22.02*Ik + 24.87*Ik^2 - 11.64*Ik^3 + 2.973*Ik^4 - 0.4135*Ik^5 + 0.02445*Ik^6)
            return current_dc / sqrt(denom)
        elseif Ik > 0.742
            denom = 1 + 0.02 * (25.62 - 133.9*Ik + 288.8*Ik^2 - 334.5*Ik^3 + 226.5*Ik^4 - 89.73*Ik^5 + 19.31*Ik^6 - 1.744*Ik^7)
            return current_dc / sqrt(denom)
        else
            return current_dc
        end
    else # For other materials like AAAC
        R_ac_factor = 1.0123
        return current_dc / sqrt(R_ac_factor)
    end
end
