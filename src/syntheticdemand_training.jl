using XGBoost, Printf

export predict, trainmodel, crossvalidate

# List of variables in the training dataset that can be use for the regression:
#
# Calendar variables:
#   :hour               hour of day
#   :month              month of year
#   :weekend01          weekend indicator
#
# Hourly temperatures:
#   :temp1              temperature in the largest population center of each region
#   :temp_top3          average temperature in the three largest population centers of each region
# 
# Monthly temperatures (season indicators):
#   :temp_monthly       average monthly temperature in the largest population center of each region
#   :ranked_month       rank of the average monthly temperature of each month (1-12)
#
# Annual temperature levels and variability:
#   :temp1_mean         average annual temperature in the largest population center of each region
#   :temp1_qlow         low annual temperature - 5% quantile of hourly temperatures
#   :temp1_qhigh        high annual temperature - 95% quantile of hourly temperatures
#
# Economic indicators:
#   :demandpercapita    level of annual average electricity demand [MWh/year/capita] in each region
#   :gdppercapita       level of annual average GDP per capita [USD(2010)/capita] in each region
# 
# Other variables (don't use these):
#    :time
#    :country

# training tests:
# ranked month  13-24  (test 1-12)
# try adding: month, temp1, temp_monthly, demandpercapita

const defaultvariables = [:hour, :weekend01, :temp_top3, :ranked_month, :temp1_mean, :temp1_qlow, :temp1_qhigh, :gdppercapita]

function predict(; variables=defaultvariables, gisregion="Europe8", scenarioyear="ssp2_2050", era_year=2018, numcenters=3, mindist=3.3,
                    nrounds=100, xgoptions...)
    df = buildtrainingdata(; gisregion=gisregion, scenarioyear=scenarioyear, era_year=era_year, numcenters=numcenters, mindist=mindist)
    regionlist = unique(df[:, :country])
    demandpercapita = df[1:8760:end, :demandpercapita]
    select!(df, variables)
    traindata = Matrix(df)
    model = trainmodel(; nrounds=nrounds, xgoptions...)
    normdemand = XGBoost.predict(model, traindata)
    numreg = length(regionlist)
    numhours = length(normdemand) ÷ numreg
    demand = reshape(normdemand, (numhours, numreg)) .* demandpercapita'
    println("\nSaving synthetic demand...")
    JLD.save(in_datafolder("output", "SyntheticDemand_$(gisregion)_$era_year.jld"), "demand", demand, compress=true)
    return demand
end

function trainmodel(; variables=defaultvariables, nrounds=100, xgoptions...)
    println("\nTraining model...")
    df_train = gettrainingdata()
    select!(df_train, variables)
    traindata = Matrix(df_train)
    normdemand = gettrainingdemand()[:, :normdemand]

    model = xgboost(traindata, nrounds; label=normdemand, xgoptions...)
end

function crossvalidate(; variables=defaultvariables, nrounds=100, metrics=["mae"], xgoptions...)
    df_train = gettrainingdata()
    regionlist = unique(df_train[:, :country])
    select!(df_train, variables)
    traindata = Matrix(df_train)
    normdemand = gettrainingdemand()[:, :normdemand]

    numreg = length(regionlist)
    params = Any[xgoptions...]

    models = nfold_cv_return(traindata, nrounds, numreg, label=normdemand, metrics=metrics, param=params)   # "rmse" or "mae"
    # pp = XGBoost.predict(models[1], traindata)
    # pp2 = XGBoost.predict(models[end], traindata)
    display(importance(models[1], string.(variables)))
    # return showstats(y, pp, pp2)
end



# Same as XGBoost.nfold_cv() but deterministic: split the training data into nfold equal parts using sequential indices.
# Also returns the models used in the last iteration.
function nfold_cv_return(data, num_boost_round::Integer = 10, nfold::Integer = 3; label = Union{},
                  param=[], metrics=[], obj = Union{}, feval = Union{}, fpreproc = Union{},
                  show_stdv = true, seed::Integer = 0, kwargs...)
    dtrain = XGBoost.makeDMatrix(data, label)
    results = String[]
    cvfolds = mknfold_deterministic(dtrain, nfold, param, metrics, fpreproc=fpreproc, kwargs = kwargs)
    for i in 1:num_boost_round
        for f in cvfolds
            XGBoost.update(f.bst, 1, f.dtrain, obj = obj)
        end
        res = XGBoost.aggcv([XGBoost.eval_set(f.bst, f.watchlist, i, feval = feval) for f in cvfolds],
                    show_stdv = show_stdv)
        push!(results, res)
        @printf(stderr, "%s\n", res)
    end
    models = [f.bst for f in cvfolds]
    return models
end

# Same as XGBoost.mknfold() but deterministic: split the training data into nfold equal parts using sequential indices
function mknfold_deterministic(dall::DMatrix, nfold::Integer, param, evals=[]; fpreproc = Union{},
                 kwargs = [])
    idx = collect(1:XGBoost.XGDMatrixNumRow(dall.handle))
    kstep = size(idx)[1] / nfold
    idset = [idx[round(Int64, (i-1) * kstep) + 1 : min(size(idx)[1],round(Int64, i * kstep))] for i in 1:nfold]
    ret = XGBoost.CVPack[]
    for k in 1:nfold
        selected = Int[]
        for i in 1:nfold
            if k != i
                selected = vcat(selected, idset[i])
            end
        end
        dtrain = XGBoost.slice(dall, selected)
        dtest = XGBoost.slice(dall, idset[k])
        if typeof(fpreproc) == Function
            dtrain, dtest, tparam = XGBoost.fpreproc(dtrain, dtest, deepcopy(param))
        else
            tparam = param
        end
        plst = vcat([itm for itm in param], [("eval_metric", itm) for itm in evals])
        plst = vcat(plst, [(string(itm[1]), string(itm[2])) for itm in kwargs])
        push!(ret, XGBoost.CVPack(dtrain, dtest, plst))
    end
    return ret
end