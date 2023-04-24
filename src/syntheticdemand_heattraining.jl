export predictheatdemand, crossvalidateheat, saveheattrainingdata, axiskeys, mean, meandrop

const defaultheatvariables = [:localhour, :temp_monthly, :temp_topN, :temp1, :month, :wind_topN]

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

function load_heatdemanddata(column="demand")    # or "profile"
    ndf = getheatdemand()
    allcountries = unique(first.(names(ndf)[2:end], 2))
    columns = allcountries .* "_$column"
    select!(ndf, ["timeUTC"; columns])
    rename!(ndf, ["timeUTC"; allcountries])
    return stack(ndf, allcountries, variable_name="country", value_name="demand")
end

function predictheatdemand(;years=2015:2019, decimals=3, shift_to_universaltime=true, variables=defaultheatvariables, demandtype="demand", iter=Int[], numcenters=5,
            nrounds=300, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"], more_xgoptions...)
    df_countrydata = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets_heatregions.csv"), DataFrame)
    countries = string.(df_countrydata.country)
    offsets = df_countrydata.offsets
    ncountries = length(countries)
    println("\nPredicting heat demand for $ncountries countries for years $years...\n")
    iter = round.(Int, iter .* 1.1)     # perform 10% more iterations than in CV training to account for extra (non-CV) year
    models = trainheatmodel(countries; demandtype, numcenters, iter,
                             nrounds, max_depth, eta, subsample, metrics, more_xgoptions...)
    #copy the model for sweden and use it for norway, add norway to countries from now on (also increase ncountries)
    # prepare_heatfuturedata is prepped to return data for norway if included in countries
    # but the csv files must be created with saveheattrainingdata(allyears=my_years, gisregion="Norway") after making a Norway region
    traindata, localtime = prepare_heatfuturedata(years, countries, variables, demandtype; numcenters)
    demand_predicted = zeros(length(localtime), ncountries)
    for c = 1:ncountries
        # No circshift needed because when2heat is trained and predicted in localtime
        demand_predicted[:,c] .= XGBoost.predict(models[c], traindata[:,:,c])
        shift_to_universaltime && (demand_predicted[:,c] = circshift(demand_predicted[:,c], round(Int, -offsets[c])))
    end
    df_results = DataFrame(demand_predicted, countries)
    # convert MWs to GWs and round to reduce file space
    for col in names(df_results)
        df_results[:, col] = round.(df_results[:, col] ./ 1000, digits=decimals)
    end
    insertcols!(df_results, 1, :localtime => localtime)
    filename = in_datafolder("output", "SyntheticHeatDemand_$(years[1])-$(years[end]).csv")
    println("\nSaving to $filename...")
    CSV.write(filename, df_results)
    nothing
end

function trainheatmodel(countries; variables=defaultheatvariables, iter=Int[], nrounds=300, demandtype="demand", numcenters=5, xgoptions...)
    ncountries = length(countries)
    if isempty(iter)
        iter = fill(nrounds, ncountries)
    end
    traindata, demand = prepare_heattrainingdata(countries, variables, demandtype; numcenters)
    println("\nTraining model for $ncountries countries...")

    models = [xgboost(traindata[:,:,c], iter[c]; label=demand[:,c], xgoptions...) for c = 1:ncountries]
end

function crossvalidateheat(country_or_all::String; variables=defaultheatvariables, demandtype="demand", numcenters=5,
            nrounds=300, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"], exit_on_loss_increase=true, more_xgoptions...)
    if country_or_all == "all"
        df_countrydata = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets_heatregions.csv"), DataFrame)
        countries = string.(df_countrydata.country)
    else
        countries = [country_or_all]
    end
    crossvalidateheat(countries; variables, demandtype, numcenters, nrounds, max_depth, eta, subsample,
                        metrics, exit_on_loss_increase, more_xgoptions...)
end

# allmodels, gain, err = crossvalidateheat("all", variables=[:localhour, :temp_monthly, :temp_topN, :temp1, :month], max_depth=6, eta=0.05);
function crossvalidateheat(countries::Vector{<:AbstractString}; variables=defaultheatvariables, demandtype="demand", numcenters=5,
            nrounds=300, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"], exit_on_loss_increase=true, more_xgoptions...)
    years = 2008:2014
    nvars, nyears, ncountries = length(variables), length(years), length(countries)
    params = Any["max_depth"=>round(Int, max_depth), "eta"=>eta, "subsample"=>subsample, "metrics"=>metrics, more_xgoptions...]
    string_variables = string.(variables)
    iv = Dict(v => i for (i, v) in enumerate(string_variables))

    traindata, demand = prepare_heattrainingdata(countries, variables, demandtype; numcenters, skipleapdays=true)
    
    allmodels = Matrix{XGBoost.Booster}(undef, nyears, ncountries)
    loss = zeros(nyears, ncountries)
    gain = zeros(nvars, nyears, ncountries)
    iterations = zeros(Int, ncountries)
    for (c, country) in enumerate(countries)
        models, minloss, iters = nfold_cv_return(traindata[:,:,c], nrounds, nyears;
                    label=demand[:,c], metrics=metrics, param=params, exit_on_loss_increase)
        allmodels[:,c] .= models
        loss[:,c] .= minloss
        iterations[c] = iters
        println("\nCrossvalidation results for country $country:")
        for (y, model) in enumerate(models)
            featuredata = importance(model, string_variables)
            for f in featuredata
                gain[iv[f.fname], y, c] += f.gain
            end
        end
        gaininfo = [variables mean(gain[:, :, c], dims=2)]
        display(sortslices(gaininfo, dims=1, by=x->x[2], rev=true))
        println()
    end

    ka_models = KeyedArray(allmodels; years, countries)
    ka_loss = KeyedArray(loss; years, countries)
    ka_gain = KeyedArray(gain; variables, years, countries)
    println("\nMean gain over all countries and years:")
    gainsummary = meandrop(ka_gain, dims=(2,3))
    display(sort(gainsummary, rev=true))
    println("\nMean test error over all countries and years: ", mean(ka_loss), "\n")

    return ka_models, ka_gain, ka_loss, iterations
end

function prepare_heatfuturedata(trainingyears, countries, variables, demandtype="demand"; numcenters=5, skipleapdays=false)
    norway_in_countries = "NO" in countries
    df_train, df_countrydata = load_heattrainingdata(; numcenters, include_norway=norway_in_countries)
    sort!(df_train, [:country, :localtime])
    trainingrows = in.(year.(df_train.localtime), Ref(trainingyears))
    df_train = df_train[trainingrows, :]

    nhours, ncountries = size(df_train,1)÷length(df_countrydata.country), length(countries)
    traindata = zeros(nhours, length(variables), ncountries)

    for (c, country) in enumerate(countries)
        keeprows = (df_train.country .== country)
        # df_train.country = get.(Ref(countrylookup), df_train.country, 0)    # convert country name to number for training
        traindata[:,:,c] .= Matrix(df_train[keeprows, variables])
    end

    return traindata, df_train.localtime[df_train.country .== "AT"]
end

function prepare_heattrainingdata(countries, variables, demandtype="demand"; numcenters=5, skipleapdays=false)
    df_train, df_countrydata = load_heattrainingdata(; numcenters)
    sort!(df_train, [:country, :localtime])
    trainingyears = 2008:2014
    trainingrows = in.(year.(df_train.localtime), Ref(trainingyears))
    df_train = df_train[trainingrows, :]   # should put df_train in sync with df_demand below

    df_demand = load_heatdemanddata(demandtype)  # "demand" or "profile"
    sort!(df_demand, [:country, :timeUTC])

    countrylookup = Dict(c => i for (i, c) in enumerate(df_countrydata.country))

    hours = df_train.timeUTC
    leapdays = (month.(hours) .== 2 .&& day.(hours) .== 29)     # skip Feb 29 (to get same hours each year)
    allrows = skipleapdays ? .!leapdays : BitVector(ones(Bool, length(hours)))
    
    nhours, ncountries = sum(allrows)÷length(df_countrydata.country), length(countries)
    traindata = zeros(nhours, length(variables), ncountries)
    demand = zeros(nhours, ncountries)

    for (c, country) in enumerate(countries)
        keeprows = allrows .&& (df_train.country .== country)
        # df_train.country = get.(Ref(countrylookup), df_train.country, 0)    # convert country name to number for training
        traindata[:,:,c] .= Matrix(df_train[keeprows, variables])

        ic = countrylookup[country]
        demand[:,c] .= circshift(df_demand.demand[keeprows], -round(Int, df_countrydata.offsets[ic]))    # shift from UTC to local time
    end

    return traindata, demand
end

function predict_heat_from_cv(models, country, variables, demandtype="demand", numcenters=5)
    trainingyears = 2008:2014
    numyears = length(trainingyears)
    traindata, demand = prepare_heattrainingdata([country], variables, demandtype; numcenters, skipleapdays=true)
    demand_predicted = similar(demand)
    numhours = 8760
    for y = 1:numyears
        rows = numhours*(y-1) + 1 : numhours*y
        demand_predicted[rows] = XGBoost.predict(models[y], traindata[rows,:])          # mean(normdemand) == 1
    end
    return reshape(demand, (numhours,numyears)), reshape(demand_predicted, (numhours,numyears))
end

function load_heattrainingdata(; numcenters=5, include_norway=false)
    filename = in_datafolder("syntheticdemand_heattrainingdata_$(numcenters)centers.csv")
    if !isfile(filename)
        saveheattrainingdata(; numcenters)
    end
    df_train = CSV.read(filename, DataFrame)
    df_countrydata = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets_heatregions.csv"), DataFrame)
    if include_norway
        df_train_norway = CSV.read(in_datafolder("syntheticdemand_heattrainingdata_Norway.csv"), DataFrame)
        df_train = vcat(df_train, df_train_norway)
        df_countrydata_norway = CSV.read(in_datafolder("syntheticdemand_timezoneoffsets_Norway.csv"), DataFrame)
        df_countrydata = vcat(df_countrydata, df_countrydata_norway)
    end
    return df_train, df_countrydata
end

function saveheattrainingdata(; allyears=2008:2019, gisregion="HeatDemandRegions", numcenters=5, mindist=3.3)
    # create_scenario_datasets("SSP2", 2020)
    println("\nCreating training dataset for synthetic demand...")
    println("(This requires ERA5 temperature data for years 2008-2014 and scenario datasets for SSP2 2020.)")
    println("\nYear $(allyears[1]):")
    df_train, offsets, _ = buildheattrainingdata(gisregion=gisregion,
                sspscenario="SSP2-34", sspyear=2020, era_year=allyears[1],
                numcenters=numcenters, mindist=mindist)
    for year in allyears[2:end]
        println("\nYear $year:")
        df, offset, _ = buildheattrainingdata(gisregion=gisregion,
                sspscenario="SSP2-34", sspyear=2020, era_year=year,
                numcenters=numcenters, mindist=mindist)
        df_train = vcat(df_train, df)
        offsets .+= offset
    end
    offsets ./= length(allyears)
    df_countrydata = DataFrame(country=unique(df_train.country), offsets=offsets)
    filename_suffix = gisregion == "HeatDemandRegions" ? "" : "_$(gisregion)"
    CSV.write(in_datafolder("syntheticdemand_timezoneoffsets_heatregions$filename_suffix.csv"), df_countrydata)
    CSV.write(in_datafolder("syntheticdemand_heattrainingdata_$(numcenters)centers$filename_suffix.csv"), df_train)
end

function createheatregions()
    ndf = getheatdemand()
    allcountries = unique(first.(names(ndf)[2:end], 2))
    exceptions = Dict("GB" => "UK", "GR" => "EL")
    nutscountries = [get(exceptions, c, c) for c in allcountries]
    regiondefinitionarray = [allcountries NUTS.(nutscountries)]
    saveregions("HeatDemandRegions", regiondefinitionarray, bbox=[34 -9; 72 32])     # exclude colonies and remote islands
end

function buildheattrainingdata(; gisregion="Europe8", sspscenario="ssp2-34", sspyear=2050, era_year=2018, numcenters=5, mindist=3.3)
    println("\nBuilding training data for $gisregion...")
    regionlist, demandpercapita, gdppercapita = 
                    makeregionaldemanddata(gisregion, sspscenario, sspyear)
    scenarioyear = "$(sspscenario[1:4])_$sspyear"
    hours, temp_popcenters = GIStemp(gisregion, scenarioyear, era_year, numcenters, mindist)
    _, wind_popcenters = GISwind_popcenters(gisregion, scenarioyear, era_year, numcenters, mindist)
    offsets, zone_maxpop, population = regional_timezone_offsets_Jan1(gisregion=gisregion, scenarioyear=scenarioyear, era_year=era_year)

    numreg, numhours = length(regionlist), length(hours)
    firsttime = ZonedDateTime.(hours[1], zone_maxpop)
    zonedtime = hcat(collect.([firsttime[i]:Hour(1):firsttime[i]+Hour(numhours-1) for i = 1:numreg])...)[:]

    println("\nShifting hourly temperatures and wind speeds from UTC to local time...")
    temperature_topN_mean = dropdims(mean(temp_popcenters, dims=3), dims=3)
    temperature_top3_mean = dropdims(mean(temp_popcenters[:,:,1:3], dims=3), dims=3)
    quantiles = hcat([quantile(temp_popcenters[:,c,1], [0.05, 0.5, 0.95]) for c = 1:numreg]...)
    shifted_temperature_topN_mean = similar(temperature_top3_mean)
    shifted_temperature_top3_mean = similar(temperature_top3_mean)
    shifted_temp1 = similar(temperature_top3_mean)
    shifted_temp2 = similar(temperature_top3_mean)
    shifted_temp3 = similar(temperature_top3_mean)

    wind_topN_mean = dropdims(mean(wind_popcenters, dims=3), dims=3)
    wind_top3_mean = dropdims(mean(wind_popcenters[:,:,1:3], dims=3), dims=3)
    shifted_wind_topN_mean = similar(wind_top3_mean)
    shifted_wind_top3_mean = similar(wind_top3_mean)
    shifted_wind1 = similar(wind_top3_mean)
    shifted_wind2 = similar(wind_top3_mean)
    shifted_wind3 = similar(wind_top3_mean)
    for r = 1:numreg
        shifted_temperature_topN_mean[:,r] = circshift(temperature_topN_mean[:,r], round(Int, offsets[r]))
        shifted_temperature_top3_mean[:,r] = circshift(temperature_top3_mean[:,r], round(Int, offsets[r]))
        shifted_temp1[:,r] = circshift(temp_popcenters[:,r,1], round(Int, offsets[r]))
        shifted_temp2[:,r] = circshift(temp_popcenters[:,r,2], round(Int, offsets[r]))
        shifted_temp3[:,r] = circshift(temp_popcenters[:,r,3], round(Int, offsets[r]))

        shifted_wind_topN_mean[:,r] = circshift(wind_topN_mean[:,r], round(Int, offsets[r]))
        shifted_wind_top3_mean[:,r] = circshift(wind_top3_mean[:,r], round(Int, offsets[r]))
        shifted_wind1[:,r] = circshift(wind_popcenters[:,r,1], round(Int, offsets[r]))
        shifted_wind2[:,r] = circshift(wind_popcenters[:,r,2], round(Int, offsets[r]))
        shifted_wind3[:,r] = circshift(wind_popcenters[:,r,3], round(Int, offsets[r]))
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
        wind_topN = shifted_wind_topN_mean[:],
        wind_top3 = shifted_wind_top3_mean[:],
        wind1 = shifted_wind1[:],
        wind2 = shifted_wind2[:],
        wind3 = shifted_wind3[:],
        localhour = hour.(zonedtime),
        month = month.(zonedtime),
        weekend01 = Int.(dayofweek.(zonedtime) .>= 6)
    )

    # sort by average monthly temperature (in popcenter 1), store rank in ranked_month
    df_monthlytemp = combine(groupby(df_time, [:country, :month]), :temp1 => mean => :temp_monthly) |>
            d -> sort!(d, [:country, :temp_monthly]) |>
            d -> insertcols!(d, 4, :ranked_month => repeat(1:12, outer=numreg))
    df_monthlywind = combine(groupby(df_time, [:country, :month]), :wind1 => mean => :wind_monthly) |>
            d -> sort!(d, [:country, :wind_monthly]) |>
            d -> insertcols!(d, 4, :ranked_month_wind => repeat(1:12, outer=numreg))

    # dataframe with regional data
    df_reg = DataFrame(country=string.(regionlist), demandpercapita=demandpercapita, gdppercapita=gdppercapita,
                        temp1_qlow=quantiles[1,:], temp1_mean=quantiles[2,:], temp1_qhigh=quantiles[3,:])

    # join everything together
    df = innerjoin(df_time, df_monthlytemp, df_monthlywind, on=[:country, :month]) |>
                    d -> innerjoin(d, df_reg, on=:country)
    return df, offsets, population
end

# Need mean hourly wind speeds in N largest population centers per region for synthetic heat demand.
# Basically just a copy/paste of GIStemp but using wind instead of temperature.
function GISwind_popcenters(gisregion::String, scenarioyear::String, era_year::Int, numcenters::Int, mindist::Float64)
    println("\nReading wind data for $gisregion...")

    res = 0.01          # resolution of auxiliary datasets [degrees per pixel]
    erares = 0.28125    # resolution of ERA5 datasets [degrees per pixel]
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange]
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    @time meanwind, wind = h5open(in_datafolder("era5wind$era_year.h5"), "r") do file
        if length(eralonranges) == 1
            file["meanwind"][eralonranges[1], eralatrange],
                file["wind"][:, eralonranges[1], eralatrange]
        else
            [file["meanwind"][eralonranges[1], eralatrange]; file["meanwind"][eralonranges[2], eralatrange]],
                [file["wind"][:, eralonranges[1], eralatrange] file["wind"][:, eralonranges[2], eralatrange]]
        end
    end

    erapop = rescale_population_to_ERA5_res(era_year, regions, offshoreregions, regionlist, lonrange, latrange, pop, wind)
    windmask = dropdims(sum(abs.(wind[1000:1100,:,:]), dims=1), dims=1) .> 0    # detect cells with any wind data (coasts?)

    numreg = length(regionlist)
    popcenters = [findpopcenters(erapop[i,:,:].*windmask, numcenters, mindist) for i = 1:numreg]
    # if any region lacks enough population centers, just repeat the ones found (only happens for tiny regions, e.g. Malta)
    for (i,p) in enumerate(popcenters)
        length(p) == 0 && error("No population centers found for region $(regionlist[i]).")
        if length(p) < numcenters
            popcenters[i] = repeat(p, outer=numcenters)[1:numcenters]
        end
    end

    hours = DateTime(era_year, 1, 1, 0) : Hour(1) : DateTime(era_year, 12, 31, 23)
    numhours = length(hours)
    numhours != size(wind,1) && error("Inconsistent number of hours.")

    wind_popcenters = [wind[h, popcenters[r][i]...] for h = 1:numhours, r = 1:numreg, i = 1:numcenters]

    return hours, wind_popcenters
end
