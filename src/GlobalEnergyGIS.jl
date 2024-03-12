module GlobalEnergyGIS

export GISwind, GISsolar, GIShydro, GIStemp, GISturbines, makedistances, annualwindindex

using MAT, HDF5, ProgressMeter, Random, Interpolations, BenchmarkTools, Images,
    Statistics, DelimitedFiles, Dates, NCDatasets, JLD, Parameters, ImageSegmentation,
    StatsBase, CSV, Distances, Printf, TimeZones, DataFrames

include("GeoArray.jl")
include("rasterize_shapefiles.jl")
include("make_auxiliary_datasets.jl")
include("era5download.jl")
include("helperfunctions.jl")
include("make_regions.jl")
include("regiondefinitions.jl")
include("solarposition.jl")
include("makewindera5.jl")
include("makesolarera5.jl")
include("maketempera5.jl")
include("makedistances.jl")
include("GISwind.jl")
include("GISsolar.jl")
include("GIShydro.jl")
include("GIStemp.jl")
include("GISturbines.jl")
include("downloaddatasets.jl")
include("sspregiondefinitions.jl")
include("mapping.jl")
include("map_test_plots.jl")
include("syntheticdemand_inputdata.jl")
include("syntheticdemand_training.jl")
include("readclimatedata.jl")

function GISsequence(regionname, regionmat)
    saveregions(regionname, regionmat)
    makedistances(regionname)
    createmaps(regionname)
    GISsolar(gisregion=regionname, plotmasks=true)
    GISwind(gisregion=regionname, plotmasks=true)
    GIShydro(gisregion=regionname)
    predictdemand(gisregion=regionname, sspscenario="ssp2-26", sspyear=2050, era_year=2018)
end

end
