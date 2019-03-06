module GlobalEnergyGIS

export GISwind

using MAT, HDF5, ProgressMeter, Interpolations, BenchmarkTools, Images, Statistics

# include("rasterize_shapefiles.jl")
include("helperfunctions.jl")
include("GISwind.jl")

end
