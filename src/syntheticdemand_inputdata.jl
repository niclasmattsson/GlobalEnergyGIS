using DataFrames

export buildtrainingdata, gettrainingdata, savetrainingdata

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
function calcdemandmultipliers(scenarioyear)
    println("Calculate demand multipliers...")
    # First read electricity demand from the SSP database into a dataframe.
    ssp = CSV.read(in_datafolder("SSP v2 Final Energy - Electricity.csv"))

    # SSP scenarios also include radiative forcing targets, e.g. SSP2-34
    # (only for IAM energy system results, not underlying population & GDP scenario)
    # We'll hardcode the 3.4 W/m2 scenario variant for now and make it a configurable option later. 
    scen = uppercase("$(scenarioyear[1:4])-34")
    year = parse(Int, scenarioyear[6:9])

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

function makeregionaldemanddata(gisregion, scenarioyear)
    res = 0.01          # resolution of auxiliary datasets [degrees per pixel]

    datasetinfo = Dict(:gisregion=>gisregion, :scenarioyear=>scenarioyear, :res=>res)
    regions, _, regionlist, _, pop, _, _, _, lonrange, latrange = read_datasets(datasetinfo)     # pop unit: people/grid cell

    gdp = JLD.load(in_datafolder("gdp_$(scenarioyear).jld"))["gdp"][lonrange,latrange]       # unit: USD(2010)/grid cell, PPP

    gadm0regions, _, gadm0regionlist, _, _ = loadregions("Global_GADM0")
    countrycodes = gadm0regions[lonrange,latrange]

    lats = (90-res/2:-res:-90+res/2)[latrange]
    cellarea = rastercellarea.(lats, res)
    popdens = pop ./ cellarea'

    natpop = getnationalpopulation(scenarioyear)        # people/GADMregion (vector, length numGADMcountries)
    demandpercapita = ieademand()./natpop * 1e6         # MWh/year/capita (vector, length numGADMcountries)
    ssp5region = make_sspregionlookup(ssp5)             # SSP 5-region names (vector, length numGADMcountries)
    demandmult = calcdemandmultipliers(scenarioyear)    # Dict: SSP region name => multiplier 

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

function buildtrainingdata(; gisregion="Europe8", scenarioyear="ssp2_2050", era_year=2018, numcenters=3, mindist=3.3)
    println("\nBuilding training data for $gisregion...")
    regionlist, demandpercapita, gdppercapita = makeregionaldemanddata(gisregion, scenarioyear)
    hours, temp_popcenters = GIStemp(gisregion, scenarioyear, era_year, numcenters, mindist)

    numreg, numhours = length(regionlist), length(hours)
    temperature_top3_mean = dropdims(mean(temp_popcenters, dims=3), dims=3)
    quantiles = hcat([quantile(temp_popcenters[:,c,1], [0.05, 0.5, 0.95]) for c = 1:numreg]...)

    # dataframe with hourly data
    df_time = DataFrame(
        time = repeat(hours, outer=numreg),
        country = repeat(string.(regionlist), inner=numhours),
        temp_top3 = temperature_top3_mean[:],
        temp1 = temp_popcenters[:,:,1][:],
        hour = repeat(hour.(hours), outer=numreg),
        month = repeat(month.(hours), outer=numreg),
        weekend01 = repeat(Int.(dayofweek.(hours).>=6), outer=numreg)
    )

    # sort by average monthly temperature (in popcenter 1), store rank in ranked_month
    df_monthlytemp = by(df_time, [:country, :month], temp_monthly = :temp1 => mean) |>
            d -> sort!(d, [:country, :temp_monthly]) |>
            d -> insertcols!(d, 4, ranked_month=repeat(13:24, outer=numreg))

    # dataframe with regional data
    df_reg = DataFrame(country=string.(regionlist), demandpercapita=demandpercapita, gdppercapita=gdppercapita,
                        temp1_qlow=quantiles[1,:], temp1_mean=quantiles[2,:], temp1_qhigh=quantiles[3,:])

    # join everything together
    df = join(df_time, df_monthlytemp, on=[:country, :month]) |>
                d -> join(d, df_reg, on=:country)

    # get rid of columns not used for the training
    # select!(df, Not([:temp1, :month]))
    println("\n\nDon't forget to create SSP2_2020 dataset.\n\n")

    return df
end

function gettrainingdata()
    # filename = in_datafolder("syntheticdemand_trainingdata.jld")
    # df_train = isfile(filename) ? JLD.load(filename, "df_train") : savetrainingdata()
    filename = in_datafolder("syntheticdemand_trainingdata.csv")
    df_train = isfile(filename) ? CSV.read(filename) : savetrainingdata()
    return df_train
end

function savetrainingdata(; numcenters=3, mindist=3.3)
    println("\nCreating training dataset for synthetic demand...")
    println("(This requires ERA5 temperature data for the year 2015.)")
    df_train = buildtrainingdata(gisregion="SyntheticDemandRegions", scenarioyear="SSP2_2020", era_year=2015, numcenters=numcenters, mindist=mindist)
    # JLD.save(in_datafolder("syntheticdemand_trainingdata.jld"), "df_train", df_train, compress=true)
    CSV.write(in_datafolder("syntheticdemand_trainingdata.csv"), df_train) 
end

gettrainingdemand() = CSV.read(in_datafolder("syntheticdemand_demanddata.csv"))
# gettrainingdemand() = JLD.load(in_datafolder("syntheticdemand_demanddata.jld"), "demand")

function saveVilledemand()
    vv = CSV.read("C:/GISdata/syntheticdemand/data/df_model_features.csv")
    hours = DateTime(2015, 1, 1, 0) : Hour(1) : DateTime(2015, 12, 31, 23)
    dem = hcat(DataFrame(time=repeat(hours, outer=44)), select(vv, [:country, :demand_total_mwh]))
    un_dem = unstack(dem, :country, :demand_total_mwh)
    rename!(un_dem, Dict("Bosnia and Herz." => "Bosnia and Herzegovina",
                                "Czech Rep." => "Czech Republic",
                                "Korea" => "South Korea"))
    select!(un_dem, [:time; sort(names(un_dem)[2:end])])
    dem2 = stack(un_dem, variable_name=:country, value_name=:demand_MW)
    select!(dem2, [3,1,2])
    meandemand = by(dem2, :country, meandemand = :demand_MW => mean)
    demand = join(dem2, meandemand, on=:country)
    insertcols!(demand, 5, normdemand=demand[:,:demand_MW]./demand[:,:meandemand])
    # CSV.write("tempdemand.csv", demand)     # next lines to get around subtype weirdness in original dataframes
    # demand = CSV.read("tempdemand.csv")
    # rm("tempdemand.csv")
    # JLD.save(in_datafolder("syntheticdemand_demanddata.jld"), "demand", demand, compress=true)
    CSV.write(in_datafolder("syntheticdemand_demanddata.csv"), demand) 
end
