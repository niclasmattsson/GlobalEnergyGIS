module GlobalEnergyGIS

export GISwind, GISsolar, GIShydro, makedistances

using MAT, HDF5, ProgressMeter, Random, Interpolations, BenchmarkTools, Images, Statistics, DelimitedFiles, Dates,
        NCDatasets, JLD, Parameters, ImageSegmentation, StatsBase, CSV, Distances, Printf, Plots

include("rasterize_shapefiles.jl")
include("regiondefinitions.jl")
include("helperfunctions.jl")
include("solarposition.jl")
include("makewindera5.jl")
include("makesolarera5.jl")
include("makedistances.jl")
include("GISwind.jl")
include("GISsolar.jl")
include("GIShydro.jl")
include("era5download.jl")
include("testio.jl")

end
