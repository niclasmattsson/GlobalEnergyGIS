module GlobalEnergyGIS

export GISwind, GISsolar

using MAT, HDF5, ProgressMeter, Interpolations, BenchmarkTools, Images, Statistics, DelimitedFiles, Dates, NCDatasets, JLD

include("rasterize_shapefiles.jl")
include("regiondefinitions.jl")
include("helperfunctions.jl")
include("solarposition.jl")
include("makewindera5.jl")
include("makesolarera5.jl")
include("GISwind.jl")
include("GISsolar.jl")
include("era5download.jl")
include("testio.jl")

end
