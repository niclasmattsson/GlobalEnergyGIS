module GlobalEnergyGIS

export GISwind, GISsolar

using MAT, HDF5, ProgressMeter, Interpolations, BenchmarkTools, Images, Statistics

# include("rasterize_shapefiles.jl")
include("helperfunctions.jl")
include("GISwind.jl")
include("GISsolar.jl")

end
