# run this file using "include" after adding Plots to the default Julia environment
using GlobalEnergyGIS, Plots, Plots.PlotMeasures

function plotlines(models, startday; variables=defaultvariables, era_year=2018)
    normdemand, normdemand_predicted, regionlist = getcvdata(models; variables=defaultvariables, era_year=era_year)
    regionlist[5] = "Bosnia & Herz."
    plot(normdemand, linecolor = "#0000FFFF", size=(1800,1350), lw=4, layout = (4,11), label="", xticks=[0,4000,8000],
    		title=reshape(regionlist, (1,44)), ylims=(0.5,1.7), titlefont=(13,:Calibri), tickfont=(11,:Calibri), bottom_margin = 10px)
    display(plot!(normdemand_predicted, linecolor = "#FF00007F", lw=4, label=""))
    # alternative colors "#00FC" and "#D008"
    plot(normdemand[startday*24 .+ (1:24*7), :], linecolor = "#0000FFFF", size=(1800,1350), layout = (4,11), label="",
    		title=reshape(regionlist, (1,44)), ylims=(0.5,1.7), titlefont=(13,:Calibri), tickfont=(11,:Calibri), bottom_margin = 10px)
    plot!(normdemand_predicted[startday*24 .+ (1:24*7), :], linecolor = "#FF00007F", label="")
end

plotly()

models = crossvalidate(variables=defaultvariables, nrounds=1000, max_depth=7, eta=0.005, subsample=0.05, metrics=["mae"])
plotlines(models, 123)
