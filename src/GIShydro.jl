# using CSV, Dates, HDF5, Statistics, Parameters

hydrooptions() = Dict(
    :gisregion => "Europe8",                    # "Europe8", "Eurasia38", "Scand3"

    :costclasses_min => [ 0,  50, 100],         # US $/MWh
    :costclasses_max => [50, 100, 999],

    :storageclasses_min => [   0, 1e-6,  12],   # weeks (discharge time)
    :storageclasses_max => [1e-6,   12, 9e9]
)

mutable struct HydroOptions
    gisregion               ::String
    costclasses_min         ::Vector{Float64}
    costclasses_max         ::Vector{Float64}
    storageclasses_min      ::Vector{Float64}
    storageclasses_max      ::Vector{Float64}
end

HydroOptions() = HydroOptions("",[],[],[],[])

function HydroOptions(d::Dict{Symbol,Any})
    options = HydroOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function getregion(lon, lat, regions, lonrange=1:36000, latrange=1:18000)
    res = 0.01
    res2 = res/2
    lons = (-180+res2:res:180-res2)[lonrange]         # longitude values (pixel center)
    lats = (90-res2:-res:-90+res2)[latrange]          # latitude values (pixel center)
    (lon < lons[1]-res2 || lon > lons[end]+res2) && return 0
    (lat > lats[1]+res2 || lat < lats[end]-res2) && return 0
    flon = findfirst(lons .> lon-res2)
    flat = findfirst(lats .<= lat+res2)
    # println(lons[flon],", ",lats[flat])
    # println(flon, " ", flat)
    if regions[flon, flat] > 0
        return regions[flon, flat]
    else
        return maximum(regions[max(1,flon-2):min(length(lons),flon+2), max(1,flat-2):min(length(lats),flat+2)])
    end
end


function GIShydro(; optionlist...)
    options = HydroOptions(merge(hydrooptions(), optionlist))

    potential, existing, hydroplants, WECcapacity, WECpotential = readhydrodatabases()

    dayspermonth, nordicprofile, eleccost, capacity, reservoirenergy, dischargetime, monthlyinflow =
        potentialhydrovars(potential)

    countries, countrynames, regions, regionlist, lonrange, latrange = read_countries_and_regions(options)

    existingcapac, existinginflowcf =
        existinginflow(existing, countries, countrynames, regions, regionlist, lonrange, latrange, hydroplants, WECcapacity,
                        dayspermonth)

    potentialcapac, potentialinflowcf, potentialmeancost, potentialmeandischargetime =
        potentialinflow(options, potential, capacity, regions, regionlist, lonrange, latrange, monthlyinflow, eleccost,
                        dischargetime, dayspermonth)

    println("\nSaving...")

    @unpack gisregion = options
    datafolder = getconfig("datafolder")
    outputfolder = joinpath(datafolder, "output")
    mkpath(outputfolder)
    filename = joinpath(outputfolder, "GISdata_hydro_$gisregion.mat")
    matopen(filename, "w") do file
        write(file, "existingcapac", existingcapac)
        write(file, "existinginflowcf", existinginflowcf)
        write(file, "potentialcapac", potentialcapac)
        write(file, "potentialinflowcf", potentialinflowcf)
        write(file, "potentialmeancost", potentialmeancost)
        write(file, "potentialmeandischargetime", potentialmeandischargetime)
    end
    nothing
end


function readhydrodatabases()
    println("\nReading hydro databases...")
    datafolder = getconfig("datafolder")

    # lat,lon,COE,Production_GWh,Lake_surface_m2,Lake_volume_m3,
    #   Qm1,Qm2,Qm3,Qm4,Qm5,Qm6,Qm7,Qm8,Qm9,Qm10,Qm11,Qm12,Qm13,ContID,BasinID,SysID,CapCost
    # original headers:
    # lat,lon,COE ($/kWh),Production (GWh),Lake surface (m2),Lake volume (m3),
    #   Qm1,Qm2,Qm3,Qm4,Qm5,Qm6,Qm7,Qm8,Qm9,Qm10,Qm11,Qm12,Qm13,ContID,BasinID,
    #   SysID (1=DiversionalCanalPower/2=RiverPower), CapCost ($/kW)
    potential = CSV.read(joinpath(datafolder, "Hydro database (Gernaat) - potential.csv"))

    # GrandID,lat,lon,Production_kWh,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13
    # original headers:
    # GrandID,lat,lon,Production_GWh,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13
    existing = CSV.read(joinpath(datafolder, "Hydro database (Gernaat) - existing (GRanD).csv"))

    # country,country_long,name,gppd_idnr,capacity_mw,latitude,longitude,fuel1,fuel2,fuel3,fuel4,commissioning_year,owner,source,url,geolocation_source,year_of_capacity_data,generation_gwh_2013,generation_gwh_2014,generation_gwh_2015,generation_gwh_2016,estimated_generation_gwh
    elecplants = CSV.read(joinpath(datafolder, "WRI - Global Power Plant Database v1.10", "global_power_plant_database.csv"), copycols=true)

    # clean up wrong coordinates
    wrong = findall((elecplants.fuel1 .== "Hydro") .& (
            (elecplants.latitude.>90) .| (elecplants.latitude.<-90) .| (elecplants.longitude.>180) .| (elecplants.longitude.<-180)
        ))
    # elecplants.url[wrong]     # 8 results, check urls online
    # 1: DDMMSS in lat and lon
    # 2-8: lat & lon flipped
    w = wrong[1]
    elecplants.latitude[w], elecplants.longitude[w] = dms2deg(elecplants.latitude[w]), dms2deg(elecplants.longitude[w])
    w = wrong[2:8]
    elecplants.latitude[w], elecplants.longitude[w] = elecplants.longitude[w], elecplants.latitude[w] 
    wrong = findall((elecplants.fuel1 .== "Hydro") .& (
            (elecplants.latitude.>90) .| (elecplants.latitude.<-90) .| (elecplants.longitude.>180) .| (elecplants.longitude.<-180)
        ))
    if !isempty(wrong)
        error("Cleanup of power plant database coordinates failed.")
    end

    # capacity_mw,latitude,longitude
    hydroplants = filter(r -> !ismissing(r.fuel1) && r.fuel1 == "Hydro", elecplants)[:,5:7]

    # Country, Capacity_MW, Pumped_MW, Other_MW, Generation_GWh, Generation_BP_GWh   % [0 = no data]
    # original headers:
    # Country, Total Hydropower Capacity (MW) in 2015, Pumped Storage Capacity (MW) in 2015, Excluding Pumped Storage (MW) in 2015, Estimated Net Hydropower Generation (GWh) in 2015, Consumption of Hydroelectricity (GWh) in 2015 (BP 2016)
    WECcapacity = CSV.read(joinpath(datafolder, "WEC hydro capacity 2015.csv"))

    # Country, Undeveloped_GWh, Potential_GWh, Utilisation
    # original headers:
    # Country, Undeveloped (GWh/year), Total Potential (GWh/year), Current Utilisation (%)
    WECpotential = CSV.read(joinpath(datafolder, "WEC hydro potentials.csv"))

    return potential, existing, hydroplants, WECcapacity, WECpotential
end


function potentialhydrovars(potential)
    dayspermonth = [Dates.daysinmonth(Date("2018-$m")) for m = 1:12]

    # No data on Nordic sites (hydro data limited to 60 degrees north), so take Nordic monthly inflow from a European model instead.
    nordicprofile = [2224, 1746, 2032, 4089, 15163, 23469, 16126, 9740, 9026, 9061, 5438, 3544]

    production = potential[:, :Production_GWh]              # GWh
    profile = Matrix(potential[:, 7:18])                    # m3/s  (columns Qm1-Qm12)
    eleccost = potential[:, :COE]                           # $/kWh
    capcost = potential[:, :CapCost]                        # $/kW
    type = potential[:, :SysID]
    reservoirsize = potential[:, :Lake_volume_m3] / 1e6     # Mm3

    cfriver = mean(profile, dims=2) ./ maximum(profile, dims=2)
    crf = CRF(0.1, 40)                                      # Gernaat assumes 10% discount rate and 40 year lifetime
    cf = capcost * crf/8760 ./ eleccost
    effic = 0.9*(type .== 1) .+ 0.7*(type .== 2)
    capacity = production ./ cf * 1000/8760                 # MW
    sortedflow = sort(profile, dims=2)                      # m3/s
    designflow = sortedflow[:, 9]                           # m3/s
    waterdensity = 997                                      # kg/m3
    grav = 9.81                                             # m/s2

    # My estimate of fallheight using eq 2-3 in the Gernaat paper.
    # Maybe ask for a new version of the database with this field included.
    fallheight = production*1e9/8760 ./ cf ./ effic ./ (designflow*waterdensity*grav)     # m

    # remove outliers
    f = findall(fallheight .> 2000)
    fallheight[f] .= 2000
    production[f] = fallheight[f] .* cf[f] .* effic[f] .* designflow[f] * waterdensity*grav*8760/1e9
    capacity[f] = production[f] ./ cf[f] * 1000/8760
    eleccost[f] = eleccost[f] .* potential[f, :Production_GWh] ./ production[f]
    capcost[f] = capcost[f] .* potential[f, :Production_GWh] ./ production[f]

    # J = kg m2/s2
    # m * m3/s * kg/m3 * m/s2 = J/s = W
    energyprofile = fallheight * dayspermonth' * 24/1e9 .* profile * waterdensity*grav    # GWh/month
    waterenergy = sum(energyprofile, dims=2)    # GWh/year

    # m * Mm3 * kg/m3 * m/s2 = MJ
    reservoirenergy = fallheight .* reservoirsize * waterdensity*grav / 3600 / 1000     # GWh
    dischargetime = reservoirenergy ./ capacity * 1000     # h
    # Matlab plots:  figure; hist(log10(reservoirenergy(reservoirenergy>0)),1000)
    # hist(log10(dischargetime(reservoirenergy>0)),1000)

    # Try calculating a monthly cf based on water inflow and turbine size.
    # Maybe not always accurate for large reservoirs since they can store between months,
    # but should work for run-of-river. Annual production and CF should still be fine.
    # (But why is my production estimate sometimes much higher than Gernaat's?)
    monthlyinflow = energyprofile .* effic
    monthlyprod = min.(monthlyinflow, capacity * dayspermonth' * 24/1000)   # GWh/month

    return dayspermonth, nordicprofile, eleccost, capacity, reservoirenergy, dischargetime, monthlyinflow
end


function read_countries_and_regions(options)
    println("Reading Global GADM country-level database...")
    countries, _, countrynames, _, _ = loadregions("Global_GADM0")
    numcountries = length(countrynames)

    @unpack gisregion = options
    println("Reading regions for $gisregion...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    # bbox = [-90 -180; 90 180]

    return countries, countrynames, regions, regionlist, lonrange, latrange
end


function existinginflow(existing, countries, countrynames, regions, regionlist, lonrange, latrange, hydroplants, WECcapacity,
                        dayspermonth)
    println("\nCalculate monthly inflow for existing hydro...");

    existingprofile = Matrix(existing[:, 5:16])             # m3/s  (columns m1-m12)

    numcountries = length(countrynames)
    capacwri = zeros(numcountries,1)        # GW
    annualcf = zeros(numcountries,1)
    scalefactor_capacwri = zeros(numcountries,1)

    for i = 1:size(hydroplants,1)
        country = getregion(hydroplants[i,:longitude], hydroplants[i,:latitude], countries)
        i == 4867 && continue   # weird data point in the Indian Ocean
        (country == 0 || country == NOREGION) && error("Can't identify country: $(hydroplants[i,:])")
        capacwri[country] += hydroplants[i,:capacity_mw] / 1e3    # GW
    end

    # calculate scale factor for national hydro capacity (WRI capac/WEC capac)
    # calculate annual cf for each country from WEC data
    for i = 1:size(WECcapacity,1)-1     # last row is "World"
        cty = findfirst(countrynames .== Symbol(WECcapacity[i,:Country]))
        cty == nothing && error("Can't find country name: $(Symbol(WECcapacity[i,:Country]))")
        scalefactor_capacwri[cty] = capacwri[cty] / WECcapacity[i,:Capacity_MW] * 1000
        annualcf[cty] = WECcapacity[i,:Generation_GWh] / WECcapacity[i,:Capacity_MW] * 1000/8760
    end

    numreg = length(regionlist)
    existingcapac = zeros(numreg,1)         # GW
    existinginflow = zeros(numreg,12)       # GWh/month

    # calculate annual generation at each site using capacity and scale factor & CF from above
    # then get monthly inflow profile from nearest site in Gernaat/Grand
    # finally calculate monthly inflow for existing sites
    for i = 1:size(hydroplants,1)
        cty = getregion(hydroplants[i,:longitude], hydroplants[i,:latitude], countries)
        reg = getregion(hydroplants[i,:longitude], hydroplants[i,:latitude], regions, lonrange, latrange) 
        (reg == 0 || reg == NOREGION || cty == 0 || cty == NOREGION) && continue
        scalefactor = scalefactor_capacwri[cty] > 0 ? scalefactor_capacwri[cty] : 1.0 
        existingcapacity = hydroplants[i,:capacity_mw] / 1e3 / scalefactor      # GW
        annualgeneration = existingcapacity * annualcf[cty] * 8760              # GWh/year
        existingcapac[reg] += existingcapacity
        
        # approx distance in km to all existing plants in Gernaat/Grand
        dist = 111*sqrt.((existing[:,:lon] .- hydroplants[i,:longitude]).^2 + (existing[:,:lat] .- hydroplants[i,:latitude]).^2)
        d, ndx = findmin(dist)
        if regionlist[reg] == "NOR"
            inflowprofile = nordicprofile
        else
            inflowprofile = existingprofile[ndx,:]
        end
        monthlyprofile = inflowprofile .* dayspermonth
        monthlyprofile = monthlyprofile / sum(monthlyprofile)
        
        monthlygeneration = annualgeneration * monthlyprofile       # GWh/month
        existinginflow[reg,:] += monthlygeneration
    end

    existinginflowcf = existinginflow ./ existingcapac ./ dayspermonth' / 24

    return existingcapac, existinginflowcf
end


function potentialinflow(options, potential, capacity, regions, regionlist, lonrange, latrange, monthlyinflow, eleccost,
            dischargetime, dayspermonth)
    println("Calculate monthly inflow for potential hydro...")

    @unpack costclasses_min, costclasses_max, storageclasses_min, storageclasses_max, gisregion = options

    ncostclasses = length(costclasses_min)
    nstorageclasses = length(storageclasses_min)

    numreg = length(regionlist)
    potentialcapac = zeros(numreg,ncostclasses,nstorageclasses)        # GW
    potentialinflow = zeros(numreg,ncostclasses,nstorageclasses,12)    # GWh/month

    potentialmeancost = zeros(numreg,ncostclasses,nstorageclasses)             # $/kWh
    potentialmeandischargetime = zeros(numreg,ncostclasses,nstorageclasses)    # hours
    nobservations = zeros(numreg,ncostclasses,nstorageclasses)

    for i = 1:size(potential,1)
        reg = getregion(potential[i,:lon], potential[i,:lat], regions, lonrange, latrange)
        (reg == 0 || reg == NOREGION) && continue

        weeks = dischargetime[i] / 168
        storageclass = findfirst((weeks .>= storageclasses_min) .& (weeks .< storageclasses_max))
        cost = eleccost[i]*1000     # US $/MWh
        costclass = findfirst((cost .>= costclasses_min) .& (cost .< costclasses_max))
        
        if regionlist[reg] == "NOR"
            inflowprofile = nordicprofile/sum(nordicprofile) * sum(monthlyinflow[i,:])
        else
            inflowprofile = monthlyinflow[i,:]
        end
        
        potentialcapac[reg,costclass,storageclass] += capacity[i]/1000
        potentialinflow[reg,costclass,storageclass,:] += inflowprofile

        potentialmeancost[reg,costclass,storageclass] += eleccost[i]
        potentialmeandischargetime[reg,costclass,storageclass] += min(8760, dischargetime[i])
        nobservations[reg,costclass,storageclass] += 1
    end

    potentialinflowcf = potentialinflow ./ potentialcapac ./ reshape(dayspermonth, (1,1,1,12)) / 24
    potentialmeancost = potentialmeancost ./ nobservations
    potentialmeandischargetime = potentialmeandischargetime ./ nobservations

    return potentialcapac, potentialinflowcf, potentialmeancost, potentialmeandischargetime
end

# function mapurl(latlon1,latlon2)
#     url = "https://www.google.com/maps/dir/?api=1&origin=$(latlon1[1]),$(latlon1[2])" *
#             "&destination=$(latlon2[1]),$(latlon2[2])&travelmode=walking"
# end
