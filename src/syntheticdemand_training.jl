using XGBoost, Printf

export predictdemand, trainmodel, crossvalidate, defaultvariables

const defaultvariables = [:localhour, :weekend01, :temp_monthly, :ranked_month, :temp_top3,
                            :temp1_mean, :temp1_qlow, :temp1_qhigh, :demandpercapita]

function predictdemand(; variables=defaultvariables, gisregion="Europe8", scenarioyear="ssp2_2050", era_year=2018, numcenters=3, mindist=3.3,
                    nrounds=100, max_depth=8, eta=0.05, subsample=0.75, metrics=["mae"], more_xgoptions...)
    df, offsets = buildtrainingdata(; gisregion=gisregion, scenarioyear=scenarioyear, era_year=era_year, numcenters=numcenters, mindist=mindist)
    regionlist = unique(df[:, :country])
    numhours = 24*daysinyear(era_year)
    demandpercapita = df[1:numhours:end, :demandpercapita]
    select!(df, variables)
    traindata = Matrix(df)
    model = trainmodel(; nrounds=nrounds, max_depth=max_depth, eta=eta, subsample=subsample, metrics=metrics, more_xgoptions...)
    normdemand = XGBoost.predict(model, traindata)
    numreg = length(regionlist)
    demand = reshape(normdemand, (numhours, numreg)) .* demandpercapita'
    println("\nConverting synthetic demand to UTC...")
    for r = 1:numreg
        demand[:,r] = circshift(demand[:,r], round(Int, -offsets[r]))
    end
    println("\nSaving...")
    JLD.save(in_datafolder("output", "SyntheticDemand_$(gisregion)_$era_year.jld"), "demand", demand, compress=true)
    nothing
end

function trainmodel(; variables=defaultvariables, nrounds=100, xgoptions...)
    df_train, offsets = loadtrainingdata()
    println("\nTraining model...")
    select!(df_train, variables)
    traindata = Matrix(df_train)
    normdemand = loaddemanddata()[:, :normdemand]

    model = xgboost(traindata, nrounds; label=normdemand, xgoptions...)
end

function crossvalidate(; variables=defaultvariables, nrounds=100, max_depth=8, eta=0.05, subsample=0.75, metrics=["mae"], more_xgoptions...)
    df_train, offsets = loadtrainingdata()
    regionlist = unique(df_train[:, :country])
    select!(df_train, variables)
    traindata = Matrix(df_train)
    normdemand = loaddemanddata()[:, :normdemand]

    numreg = length(regionlist)
    params = Any["max_depth"=>round(Int, max_depth), "eta"=>eta, "subsample"=>subsample, "metrics"=>metrics, more_xgoptions...]

    models = nfold_cv_return(traindata, nrounds, numreg; label=normdemand, metrics=metrics, param=params)   # "rmse" or "mae"
    display(importance(models[1], string.(variables)))

    return models
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