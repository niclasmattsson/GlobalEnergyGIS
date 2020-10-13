using DataFrames

export buildtrainingdata, loadtrainingdata, loaddemanddata, savetrainingdata

# Can call this with ssp5 or ssp32 (global constants)
# Since not all GADM level 0 regions are assigned to SSP regions, some regions are missed
# (mostly small island nations and the like). This captures 99.85% of global population.
function make_sspregionlookup(ssp)
    println("Lookup SSP R5.2 region names for each GADM level 0 region...")
    _, _, regionlist, _, _ = loadregions("Global_GADM0")
    sspregions = fill("", length(regionlist))
    for (sspregionname, gadmregionnames) in ssp
        for gname in gadmregionnames
            index = findfirst(isequal(Symbol(gname)), regionlist)
            if index == nothing
                error("Region $gname is not a GADM level 0 region.")
            end
            sspregions[index] = sspregionname
        end
    end
    return sspregions # Vector{String}(length numcountries, indexed by gadm country code)
end

function getnationalpopulation(scenarioyear)
    filename = in_datafolder("nationalpopulation_$scenarioyear.jld")
    natpop = isfile(filename) ? JLD.load(filename, "natpop") : savenationalpopulation(scenarioyear)
    return natpop
end

# load population dataset and sum globally by country
# Vector{Float64}(length numcountries, indexed by gadm country code)
function savenationalpopulation(scenarioyear)
    println("Calculating population in all GADM level 0 regions...")
    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")   # persons per pixel
    regions, _, regionlist, lonrange, latrange = loadregions("Global_GADM0")
    natpop = zeros(length(regionlist))
    for j in latrange
        for i in lonrange
            reg = regions[i,j]
            if reg > 0
                natpop[reg] += pop[i,j]
            end
        end
    end
    JLD.save(in_datafolder("nationalpopulation_$scenarioyear.jld"), "natpop", natpop, compress=true)
    return natpop
end

# Get current national electricity demand from IEA online statistics (actually "domestic supply").
# https://www.iea.org/statistics/?country=MOROCCO&year=2016&category=Electricity&indicator=ElecGenByFuel&mode=table&dataTable=ELECTRICITYANDHEAT
function ieademand()
    println("Get current national electricity demand from IEA statistics...")
    iea = CSV.read(in_datafolder("ieademand_2016.csv"))      # GWh/year
    _, _, regionlist, _, _ = loadregions("Global_GADM0")
    nationaldemand = zeros(length(regionlist))
    for row in eachrow(iea)
        country, demand = row
        country == uppercase(country) && continue   # skip IEA aggregated regions in uppercase, e.g. "ASIA"
        f = findfirst(regionlist .== Symbol(country))
        f == nothing && error("IEA country $country missing in GADM regions.")
        nationaldemand[f] = demand/1000
    end
    return nationaldemand        # TWh/year 
end

# Lookup electricity demand from a SSP database dataframe. Return a single dataframe row.
ssplookup(ssp, model, scen, region) =
    ssp[(ssp[:,:MODEL] .== model) .& (ssp[:,:SCENARIO] .== scen) .& (ssp[:,:REGION] .== region), :][1,:]

# Calculate regional multiplers for demand from 2016 to target year
function calcdemandmultipliers(sspscenario::String, year::Int)
    println("Calculate demand multipliers...")
    # First read electricity demand from the SSP database into a dataframe.
    ssp = CSV.read(in_datafolder("SSP v2 Final Energy - Electricity.csv"))

    # SSP scenarios also include radiative forcing targets, e.g. SSP2-34
    # (only for IAM energy system results, not underlying population & GDP scenario)
    # We'll hardcode the 3.4 W/m2 scenario variant for now and make it a configurable option later. 
    scen = uppercase(sspscenario)   # e.g. "SSP2-34"

    # We'll take the average result from IMAGE and MESSAGE models.
    demandmult2050 = Dict{String, Float64}()
    for sspreg in keys(ssp5)
        reg = "R5.2$sspreg"
        image = ssplookup(ssp, "IMAGE", scen, reg)
        message = ssplookup(ssp, "MESSAGE-GLOBIOM", scen, reg)
        # Take the year 2020 to represent 2016, then adjust for this difference.
        image_mult = image[Symbol(year)] / image[Symbol(2020)]
        message_mult = message[Symbol(year)] / message[Symbol(2020)]
        demandmult2050[sspreg] = ((image_mult + message_mult)/2) ^ ((year-2016)/(year-2020))
    end
    return demandmult2050
end

function makeregionaldemanddata(gisregion, sspscenario::String, year::Int)
    res = 0.01          # resolution of auxiliary datasets [degrees per pixel]

    scenarioyear = "$(sspscenario[1:4])_$year"
    datasetinfo = Dict(:gisregion=>gisregion, :scenarioyear=>scenarioyear, :res=>res)
    regions, _, regionlist, _, pop, _, _, _, lonrange, latrange = read_datasets(datasetinfo)     # pop unit: people/grid cell

    gdp = JLD.load(in_datafolder("gdp_$(scenarioyear).jld"))["gdp"][lonrange,latrange]       # unit: USD(2010)/grid cell, PPP

    gadm0regions, _, gadm0regionlist, _, _ = loadregions("Global_GADM0")
    countrycodes = gadm0regions[lonrange,latrange]

    lats = (90-res/2:-res:-90+res/2)[latrange]
    cellarea = rastercellarea.(lats, res)
    popdens = pop ./ cellarea'

    natpop = getnationalpopulation(scenarioyear)            # people/GADMregion (vector, length numGADMcountries)
    demandpercapita = ieademand()./natpop * 1e6             # MWh/year/capita (vector, length numGADMcountries)
    ssp5region = make_sspregionlookup(ssp5)                 # SSP 5-region names (vector, length numGADMcountries)
    demandmult = calcdemandmultipliers(sspscenario, year)   # Dict: SSP region name => multiplier 

    numreg = length(regionlist)
    regionaldemand = zeros(numreg)
    regionalpop = zeros(numreg)
    regionalgdp = zeros(numreg)

    println("Calculating annual electricity demand for model regions...")
    nlons, nlats = size(regions)
    updateprogress = Progress(nlats, 1)
    for j = 1:nlats
        for i = 1:nlons
            reg = regions[i,j]
            countrycode = countrycodes[i,j]
            if reg > 0 && reg != NOREGION && countrycode > 0
                sspreg = ssp5region[countrycode]
                if isempty(sspreg)
                    # error("Oops, no SSP region assigned to region $(gadm0regionlist[countrycode]).")
                    continue
                end
                regionaldemand[reg] += demandpercapita[countrycode] * pop[i,j]/1e6 * demandmult[sspreg]     # TWh/year
                regionalpop[reg] += pop[i,j]        # unit: people
                regionalgdp[reg] += gdp[i,j]        # unit: USD(2010) 
            end
        end
        next!(updateprogress)
    end
    regionaldemandpercapita = regionaldemand ./ regionalpop * 1e6   # MWh/year/capita
    regionalgdppercapita = regionalgdp ./ regionalpop               # USD(2010)/capita

    return regionlist, regionaldemandpercapita, regionalgdppercapita
end

function buildtrainingdata(; gisregion="Europe8", sspscenario="ssp2-34", sspyear=2050, era_year=2018, numcenters=3, mindist=3.3)
    println("\nBuilding training data for $gisregion...")
    regionlist, demandpercapita, gdppercapita = 
                    makeregionaldemanddata(gisregion, sspscenario, sspyear)
    scenarioyear = "$(sspscenario[1:4])_$sspyear"
    hours, temp_popcenters = GIStemp(gisregion, scenarioyear, era_year, numcenters, mindist)
    offsets, zone_maxpop, population = regional_timezone_offsets_Jan1(gisregion=gisregion, scenarioyear=scenarioyear, era_year=era_year)

    numreg, numhours = length(regionlist), length(hours)
    firsttime = ZonedDateTime.(hours[1], zone_maxpop)
    zonedtime = hcat(collect.([firsttime[i]:Hour(1):firsttime[i]+Hour(8759) for i = 1:numreg])...)[:]

    println("\nShifting hourly temperatures from UTC to local time...")
    temperature_top3_mean = dropdims(mean(temp_popcenters, dims=3), dims=3)
    quantiles = hcat([quantile(temp_popcenters[:,c,1], [0.05, 0.5, 0.95]) for c = 1:numreg]...)
    shifted_temperature_top3_mean = similar(temperature_top3_mean)
    shifted_temp1 = similar(temperature_top3_mean)
    for r = 1:numreg
        shifted_temperature_top3_mean[:,r] = circshift(temperature_top3_mean[:,r], round(Int, offsets[r]))
        shifted_temp1[:,r] = circshift(temp_popcenters[:,r,1], round(Int, offsets[r]))
    end

    # dataframe with hourly data
    df_time = DataFrame(
        localtime = DateTime.(zonedtime, Local),
        country = repeat(string.(regionlist), inner=numhours),
        temp_top3 = shifted_temperature_top3_mean[:],
        temp1 = shifted_temp1[:],
        localhour = hour.(zonedtime),
        month = month.(zonedtime),
        weekend01 = Int.(dayofweek.(zonedtime) .>= 6)
    )

    # sort by average monthly temperature (in popcenter 1), store rank in ranked_month
    df_monthlytemp = by(df_time, [:country, :month], temp_monthly = :temp1 => mean) |>
            d -> sort!(d, [:country, :temp_monthly]) |>
            d -> insertcols!(d, 4, ranked_month=repeat(1:12, outer=numreg))

    # dataframe with regional data
    df_reg = DataFrame(country=string.(regionlist), demandpercapita=demandpercapita, gdppercapita=gdppercapita,
                        temp1_qlow=quantiles[1,:], temp1_mean=quantiles[2,:], temp1_qhigh=quantiles[3,:])

    # join everything together
    df = join(df_time, df_monthlytemp, on=[:country, :month]) |>
                    d -> join(d, df_reg, on=:country)
    return df, offsets, population
end

function loadtrainingdata()
    filename = in_datafolder("syntheticdemand_trainingdata.csv")
    if !isfile(filename)
        savetrainingdata()
    end
    df_train = CSV.read(filename)
    offsets = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets.csv"))[:,1]
    return df_train, offsets
end

function savetrainingdata(; numcenters=3, mindist=3.3)
    create_scenario_datasets("SSP2", 2020)
    println("\nCreating training dataset for synthetic demand...")
    println("(This requires ERA5 temperature data for the year 2015 and scenario datasets for SSP2 2020.)")
    df_train, offsets, _ = buildtrainingdata(gisregion="SyntheticDemandRegions",
                sspscenario="SSP2-34", sspyear=2020, era_year=2015,
                numcenters=numcenters, mindist=mindist)
    CSV.write(in_datafolder("syntheticdemand_timezoneoffsets.csv"), DataFrame(offsets=offsets))
    CSV.write(in_datafolder("syntheticdemand_trainingdata.csv"), df_train)
end

loaddemanddata() = CSV.read(in_datafolder("syntheticdemand_demanddata.csv"))

function regional_timezone_offsets_Jan1(; gisregion="Europe8", scenarioyear="ssp2_2050", era_year=2018)
    println("\nCalculating population-weighted regional time zone offsets...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    tzindices, tznames = loadtimezones(lonrange, latrange)
    popscale = 1e5
    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange] ./ popscale   # scale down for better precision
    numreg = length(regionlist)
    numhours = 24*daysinyear(era_year)
    offsets = zeros(numreg)
    population = zeros(numreg)
    zone_maxpop = fill(tz"Europe/London", numreg)
    updateprogress = Progress(numreg, 1)
    for r = 1:numreg
        reg = (regions .== r)
        regindices = tzindices[reg]
        regpop = pop[reg]
        zoneindices = unique(regindices)
        weightedoffset = 0.0
        zones = TimeZone[]
        pops = Float64[]
        for idx in zoneindices
            tzname = tznames[idx]
            zone = tzname[1:3] == "Etc" ? TimeZone(tzname, TimeZones.Class(:LEGACY)) : TimeZone(tzname)
            push!(zones, zone)
            firsthour = ZonedDateTime(DateTime(era_year,1,1,0), zone)
            hours = firsthour : Hour(1) : firsthour + Hour(numhours-1)
            offset = tzoffset(hours[1])
            localpop = sum(regpop[regindices.==idx])
            push!(pops, localpop)
            weightedoffset += localpop * offset
        end
        _, i = findmax(pops)    # find time zone with the most population
        zone_maxpop[r] = zones[i]
        population[r] = sum(pops) * popscale    # scale up again
        offsets[r] = weightedoffset / sum(pops)
        next!(updateprogress)
    end

    return offsets, zone_maxpop, population
end

function saveVilledemand()
    vv = CSV.read(in_datafolder("syntheticdemand", "data", "df_model_features.csv"))
    hours = DateTime(2015, 1, 1, 0) : Hour(1) : DateTime(2015, 12, 31, 23)
    demand1D = hcat(DataFrame(hours=repeat(hours, outer=44)), select(vv, [:country, :demand_total_mwh]))
    demand2D = unstack(demand1D, :country, :demand_total_mwh)
    rename!(demand2D, Dict("Bosnia and Herz." => "Bosnia and Herzegovina",
                                "Czech Rep." => "Czech Republic",
                                "Korea" => "South Korea"))
    select!(demand2D, [:hours; sort(names(demand2D)[2:end])])

    dem2 = stack(demand2D, variable_name=:country, value_name=:demand_MW)
    select!(dem2, [3,1,2])
    meandemand = by(dem2, :country, meandemand = :demand_MW => mean)
    demand = join(dem2, meandemand, on=:country)

    insertcols!(demand, 5, normdemand=demand[:,:demand_MW]./demand[:,:meandemand])
    CSV.write(in_datafolder("syntheticdemand_demanddata.csv"), demand)
end

tzoffset(tz::FixedTimeZone) = tz.offset.std.value / 3600
tzoffset(tz::VariableTimeZone) = tzoffset(tz.transitions[end].zone)
tzoffset(dt::ZonedDateTime) = TimeZones.value(dt.zone.offset) / 3600

function timezonetest()
    # zones = TimeZone.(["Europe/London", "Europe/Stockholm", "Europe/Helsinki", "Europe/Minsk"])
    zones = TimeZone.(["Australia/Perth", "Australia/Adelaide", "Australia/Darwin", "Australia/Brisbane", "Australia/Sydney"])
    # all hours in 2015 for each zone in local time
    hours = hcat([ZonedDateTime(DateTime(2015,1,1,0), z) : Hour(1) : ZonedDateTime(DateTime(2015,12,31,23), z) for z in zones]...)
    # shift hours by the time zone offset to ensure simultaneous rows
    # take the offset of ZonedDateTime in hour 1 (in case hour 1 is DST), not the UTC time zone offset
    shifted = hcat([circshift(hours[:,i], round(Int, -tzoffset(hours[1,i]))) for i = 1:length(zones)]...)
    # convert to UTC and ensure all columns identical (except last 24 hours because of shifted hours from previous year)
    utctime = DateTime.(shifted, UTC)
    # @assert all([all(utctime[1:end-24,1] .== utctime[1:end-24,i]) for i=2:length(zones)])  # doesn't work for non-integer time offsets
    # show the local time in each zone around the shift to summer time (daylight savings) in the EU
    localtime = Time.(DateTime.(shifted, Local))
    # display(localtime[2085:2095,:])   # EU
    display(localtime[6615:6620,:])     # Australia
    return hours, shifted, utctime, localtime
end
