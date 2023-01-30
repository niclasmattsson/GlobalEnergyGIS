export predictheatdemand, trainheatmodel, crossvalidateheat

function getheatdemand()
    df = CSV.read(in_datafolder("when2heat.csv"), DataFrame, delim=';', decimal=',')
    allcountries = unique(first.(names(df)[3:end], 2))
    countries = setdiff(allcountries, ["CH", "NO"])
    profilenames = ["heat_profile_$(a)_$b" for a in ["space", "water"] for b in ["COM", "MFH", "SFH"]]
    for cc in countries
        profiles = ["$(cc)_$pp" for pp in profilenames]
        df[!, "$(cc)_profile"] .= sum(Matrix(df[!, profiles]), dims=2)
    end
    timeUTC = DateTime.(first.(df.utc_timestamp, 19))
    demandcols = ["$(cc)_heat_demand_total" for cc in countries]
    profilecols = ["$(cc)_profile" for cc in countries]
    years = 2008 .<= year.(timeUTC) .<= 2014
    ndf = df[years, sort([demandcols; profilecols])]
    rename!(ndf, Dict("$(cc)_heat_demand_total" => "$(cc)_demand" for cc in countries))
    insertcols!(ndf, 1, :timeUTC => timeUTC[years]) 
    disallowmissing!(ndf)
    return ndf
end

function loadheatdemanddata(column="demand")    # or "profile"
    ndf = getheatdemand()
    allcountries = unique(first.(names(ndf)[2:end], 2))
    columns = allcountries .* "_$column"
    select!(ndf, ["timeUTC"; columns])
    rename!(ndf, ["timeUTC"; allcountries])
    return stack(ndf, allcountries, variable_name="country", value_name="demand")
end

function predictheatdemand(; variables=defaultvariables, gisregion="Europe8",
            sspscenario="ssp2-34", sspyear=2050, era_year=2018, numcenters=3, mindist=3.3,
            nrounds=100, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"], more_xgoptions...)
    # df, offsets, pop = buildheattrainingdata(; gisregion=gisregion, sspscenario=sspscenario,
    #         sspyear=sspyear, era_year=era_year, numcenters=numcenters, mindist=mindist)
    # regionlist = unique(df[:, :country])
    # numhours = 24*daysinyear(era_year)
    # demandpercapita = df[1:numhours:end, :demandpercapita]      # MWh/year/capita
    # select!(df, variables)
    
    # traindata = Matrix(df)
    # model = trainmodel(; nrounds=nrounds, max_depth=max_depth, eta=eta, subsample=subsample, metrics=metrics, more_xgoptions...)
    # normdemand = XGBoost.predict(model, traindata)              # mean(normdemand) == 1
    # numreg = length(regionlist)
    # demand = reshape(normdemand, (numhours, numreg)) .* (demandpercapita/8760)' .* pop'    # MW
    # println("\nConverting synthetic demand to UTC...")
    # for r = 1:numreg
    #     demand[:,r] = circshift(demand[:,r], round(Int, -offsets[r]))
    # end
    # println("\nSaving...")
    # JLD.save(in_datafolder("output",
    #         "SyntheticDemand_$(gisregion)_$sspscenario-$(sspyear)_$era_year.jld"),
    #         "demand", demand, compress=true)
    # nothing
end

function trainheatmodel(; variables=defaultvariables, nrounds=100, xgoptions...)
    df_train, offsets = loadheattrainingdata()
    println("\nTraining model...")
    select!(df_train, variables)
    traindata = Matrix(df_train)
    demand = loadheatdemanddata()[!, :demand]

    model = xgboost(traindata, nrounds; label=normdemand, xgoptions...)
end

# crossvalidateheat("FR", variables=[:localhour, :weekend01, :temp_monthly, :temp_top3, :temp_topN, :ranked_month, :month, :temp2, :temp1, :temp3, :temp1_mean, :temp1_qlow, :temp1_qhigh], nrounds=150, max_depth=8, eta=0.05, subsample=0.75, metrics=["mae"])
function crossvalidateheat(country; variables=defaultvariables, demandtype="demand", nrounds=100, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"], more_xgoptions...)
    df_train, offsets = loadheattrainingdata()
    sort!(df_train, [:country, :localtime])
    hours = df_train.timeUTC
    numyears = year(df_train.localtime[end]) - year(df_train.localtime[1]) + 1
    keeprows = .!(month.(hours) .== 2 .&& day.(hours) .== 29)    # don't use leap days here (to get same hours each year)
    if !isempty(country)
        keeprows = keeprows .&& (df_train.country .== country)
    end

    countries = unique(df_train.country)
    countrylookup = Dict(c => i for (i, c) in enumerate(countries))
    df_train.country = get.(Ref(countrylookup), df_train.country, 0)    # convert country name to number for training

    select!(df_train, variables)
    traindata = Matrix(df_train[keeprows, :])

    df_demand = loadheatdemanddata(demandtype)  # "demand" or "profile"
    sort!(df_demand, [:country, :timeUTC])
    demand = df_demand.demand[keeprows]
    demand .= circshift(demand, -round(Int, offsets[countrylookup[country]]))    # shift from UTC to local time
    # for r = 1:length(countries)
    #     rows = (keeprows .&& df_demand.country .== countries[r])
    #     demand[rows] .= circshift(demand[rows], round(Int, offsets[r]))    # shift from UTC to local time
    # end

    params = Any["max_depth"=>round(Int, max_depth), "eta"=>eta, "subsample"=>subsample, "metrics"=>metrics, more_xgoptions...]

    models = nfold_cv_return(traindata, nrounds, numyears; label=demand, metrics=metrics, param=params)   # "rmse" or "mae"
    display(importance(models[1], string.(variables)))
    println("\n\nOnly prints feature importance for 2008, collect all years for all countries.\n\n")

    return models
end

function loadheattrainingdata()
    filename = in_datafolder("syntheticdemand_heattrainingdata_new.csv")
    if !isfile(filename)
        saveheattrainingdata()
    end
    df_train = CSV.read(filename, DataFrame)
    offsets = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets_heatregions_new.csv"), DataFrame)[:,1]
    return df_train, offsets
end

function saveheattrainingdata(; numcenters=5, mindist=3.3)
    # create_scenario_datasets("SSP2", 2020)
    println("\nCreating training dataset for synthetic demand...")
    println("(This requires ERA5 temperature data for years 2008-2014 and scenario datasets for SSP2 2020.)")
    allyears = 2008:2014
    println("\nYear $(allyears[1]):")
    df_train, offsets, _ = buildheattrainingdata(gisregion="HeatDemandRegions",
                sspscenario="SSP2-34", sspyear=2020, era_year=allyears[1],
                numcenters=numcenters, mindist=mindist)
    for year in allyears[2:end]
        println("\nYear $year:")
        df, offset, _ = buildheattrainingdata(gisregion="HeatDemandRegions",
                sspscenario="SSP2-34", sspyear=2020, era_year=year,
                numcenters=numcenters, mindist=mindist)
        df_train = vcat(df_train, df)
        offsets .+= offset
    end
    offsets ./= length(allyears)
    CSV.write(in_datafolder("syntheticdemand_timezoneoffsets_heatregions_new.csv"), DataFrame(offsets=offsets))
    CSV.write(in_datafolder("syntheticdemand_heattrainingdata_new.csv"), df_train)
end

function createheatregions()
    ndf = getheatdemand()
    allcountries = unique(first.(names(ndf)[2:end], 2))
    exceptions = Dict("GB" => "UK", "GR" => "EL")
    nutscountries = [get(exceptions, c, c) for c in allcountries]
    regiondefinitionarray = [allcountries NUTS.(nutscountries)]
    saveregions("HeatDemandRegions", regiondefinitionarray, bbox=[34 -9; 72 32])     # exclude colonies and remote islands
end

function buildheattrainingdata(; gisregion="Europe8", sspscenario="ssp2-34", sspyear=2050, era_year=2018, numcenters=3, mindist=3.3)
    println("\nBuilding training data for $gisregion...")
    regionlist, demandpercapita, gdppercapita = 
                    makeregionaldemanddata(gisregion, sspscenario, sspyear)
    scenarioyear = "$(sspscenario[1:4])_$sspyear"
    hours, temp_popcenters = GIStemp(gisregion, scenarioyear, era_year, numcenters, mindist)
    offsets, zone_maxpop, population = regional_timezone_offsets_Jan1(gisregion=gisregion, scenarioyear=scenarioyear, era_year=era_year)

    numreg, numhours = length(regionlist), length(hours)
    firsttime = ZonedDateTime.(hours[1], zone_maxpop)
    zonedtime = hcat(collect.([firsttime[i]:Hour(1):firsttime[i]+Hour(numhours-1) for i = 1:numreg])...)[:]

    println("\nShifting hourly temperatures from UTC to local time...")
    temperature_topN_mean = dropdims(mean(temp_popcenters, dims=3), dims=3)
    temperature_top3_mean = dropdims(mean(temp_popcenters[:,:,1:3], dims=3), dims=3)
    quantiles = hcat([quantile(temp_popcenters[:,c,1], [0.05, 0.5, 0.95]) for c = 1:numreg]...)
    shifted_temperature_topN_mean = similar(temperature_top3_mean)
    shifted_temperature_top3_mean = similar(temperature_top3_mean)
    shifted_temp1 = similar(temperature_top3_mean)
    shifted_temp2 = similar(temperature_top3_mean)
    shifted_temp3 = similar(temperature_top3_mean)
    for r = 1:numreg
        shifted_temperature_topN_mean[:,r] = circshift(temperature_topN_mean[:,r], round(Int, offsets[r]))
        shifted_temperature_top3_mean[:,r] = circshift(temperature_top3_mean[:,r], round(Int, offsets[r]))
        shifted_temp1[:,r] = circshift(temp_popcenters[:,r,1], round(Int, offsets[r]))
        shifted_temp2[:,r] = circshift(temp_popcenters[:,r,2], round(Int, offsets[r]))
        shifted_temp3[:,r] = circshift(temp_popcenters[:,r,3], round(Int, offsets[r]))
    end

    # dataframe with hourly data
    df_time = DataFrame(
        timeUTC = DateTime.(zonedtime, UTC),
        localtime = DateTime.(zonedtime, Local),
        country = repeat(string.(regionlist), inner=numhours),
        temp_topN = shifted_temperature_topN_mean[:],
        temp_top3 = shifted_temperature_top3_mean[:],
        temp1 = shifted_temp1[:],
        temp2 = shifted_temp2[:],
        temp3 = shifted_temp3[:],
        localhour = hour.(zonedtime),
        month = month.(zonedtime),
        weekend01 = Int.(dayofweek.(zonedtime) .>= 6)
    )

    # sort by average monthly temperature (in popcenter 1), store rank in ranked_month
    df_monthlytemp = combine(groupby(df_time, [:country, :month]), :temp1 => mean => :temp_monthly) |>
            d -> sort!(d, [:country, :temp_monthly]) |>
            d -> insertcols!(d, 4, :ranked_month => repeat(1:12, outer=numreg))

    # dataframe with regional data
    df_reg = DataFrame(country=string.(regionlist), demandpercapita=demandpercapita, gdppercapita=gdppercapita,
                        temp1_qlow=quantiles[1,:], temp1_mean=quantiles[2,:], temp1_qhigh=quantiles[3,:])

    # join everything together
    df = innerjoin(df_time, df_monthlytemp, on=[:country, :month]) |>
                    d -> innerjoin(d, df_reg, on=:country)
    return df, offsets, population
end
