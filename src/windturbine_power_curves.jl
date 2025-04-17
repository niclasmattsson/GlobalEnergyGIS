function speed2capacityfactor(windspeed, powercurve)
    if 2*windspeed >= length(powercurve) || windspeed < 0
        return 0.0
    end
    fw = floor(Int, 2 * windspeed)  # *2 because the power curves are given in 0.5 m/s steps
    frac = windspeed - fw
    return (1-frac).*powercurve[fw+1] + frac.*powercurve[fw+2]
end

# new power curves are given in 0.5 m/s steps, the old one has 1 m/s steps
function interpolate_old_powercurve(pc_old)
    pc = zeros(length(pc_old)*2 - 1)
    pc[1:2:end] .= pc_old
    pc[2:2:end-1] .= (pc_old[1:end-1] + pc_old[2:end]) ./ 2
    return pc
end

# Turbine curves for wind farms, including electrical and wake losses.
const powercurves = Dict(
    # Original GlobalEnergyGIS turbine curve, 0 - 29 m/s in 1 m/s steps
    "Vestas V112-3" => interpolate_old_powercurve([
        0.0, 0.0014, 0.0071, 0.0229, 0.0545, 0.1067, 0.1831, 0.2850, 0.4085, 0.5434,
        0.6744, 0.7847, 0.8614, 0.9048, 0.9266, 0.9353, 0.9373, 0.9375, 0.9375, 0.9375,
        0.9375, 0.9375, 0.9375, 0.9311, 0.8683, 0.6416, 0.2948, 0.0688, 0.0063, 0.0
    ]),

    # A floating offshore turbine for the Mareld project, 0 - 40 m/s in 0.5 m/s steps
    # round.(generateSyntheticPowerCurve(; specificPower=0, cutin=3, cutout=31, v_rated=11.1,
    #           elec_losses=0.06, wake_losses=0.115, wakemodel=:shift, v_sample=0:0.5:40), digits=4)
    "Vestas V136-15" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0006, 0.0039, 0.0133, 0.0308, 0.057, 0.092, 0.1363, 0.1908, 0.2562,
        0.3329, 0.4202, 0.5153, 0.6129, 0.7054, 0.7853, 0.8473, 0.8901, 0.9161, 0.9299, 0.9363, 0.9388, 0.9397, 0.9399,
        0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
        0.94, 0.94, 0.94, 0.94, 0.94, 0.9397, 0.9389, 0.9362, 0.9285, 0.9103, 0.874, 0.8118, 0.72, 0.6021, 0.47, 0.3398,
        0.2264, 0.1386, 0.0778, 0.0401, 0.0189, 0.0082, 0.0033, 0.0012, 0.0004, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ]
)
