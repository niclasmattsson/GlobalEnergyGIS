export GeoArray, lonlat_index, lookup, getlon, getlat, lonindex, latindex, crop

const EPS = 1e-8

# Assumes the latitudes of the internal array are in reverse order.
struct GeoArray{T,N} <: AbstractArray{T,N}
    arr::AbstractArray{T,N}
    res::Float64
    lonlim::Tuple{Float64,Float64}
    latlim::Tuple{Float64,Float64}

    # Check if lon-lat limits are valid before constructing.
    function GeoArray(arr, res, lonlim, latlim)
        calcsize = (lonlim[2] - lonlim[1]) / res, (latlim[2] - latlim[1]) / res
        if all(size(arr)[1:2] .== calcsize)      # maybe ≈
            return new{eltype(arr), ndims(arr)}(arr, res, lonlim, latlim)
        else
            error("lon-lat limits incompatible with GeoArray size")
        end
    end
end

GeoArray(arr, res, extent::Vector{<:Real}) = GeoArray(arr, res, (extent[1],extent[3]), (extent[2],extent[4]))

function Base.show(io::IO, ::MIME"text/plain", a::GeoArray)
    summary(io, a)
    println(io, ": res = $(a.res), lonlim = $(a.lonlim), latlim = $(a.latlim)")
    print("Array: ")
    show(io, "text/plain", a.arr)
end

Base.size(a::GeoArray) = size(a.arr)
Base.@propagate_inbounds Base.getindex(a::GeoArray, i::Int) = getindex(a.arr, i)
Base.IndexStyle(::Type{<:GeoArray{T,N}}) where T where N = IndexStyle(Array{T,N})

getlon(a::GeoArray, i::Int) = a.lonlim[1] + a.res*(i - 0.5)  # i - 1 + 0.5
getlat(a::GeoArray, i::Int) = a.latlim[2] - a.res*(i - 0.5)

lonindex(a::GeoArray, lon) = floor(Int, (lon - a.lonlim[1])/a.res) + 1
latindex(a::GeoArray, lat) = floor(Int, (a.latlim[2] - lat)/a.res) + 1

lonlat_index(a::GeoArray, lon, lat, t; restrict=false) = CartesianIndex(lonlat_index(a, lon, lat; restrict), t)
function lonlat_index(a::GeoArray, lon, lat; restrict=false)
    ndx = CartesianIndex(lonindex(a, lon), latindex(a, lat))
    return restrict ? clamp(ndx, zero(ndx), CartesianIndex(size(a)[1:2])) : ndx
end

lookup(a::GeoArray, args...; restrict=false) = a[lonlat_index(a, args...; restrict)]

lon_indices_within(a::GeoArray, lonlow, lonhigh) = 
        max(1, lonindex(a, lonlow)):min(size(a,1), lonindex(a, lonhigh) - 1)
lat_indices_within(a::GeoArray, latlow, lathigh) =
        max(1, latindex(a, lathigh)):min(size(a,2), latindex(a, latlow) - 1)

function crop(a::GeoArray, targetextent::Vector{<:Real})
    ilons = lonindex(a, targetextent[1]) : lonindex(a, targetextent[3])
    ilats = latindex(a, targetextent[4]) : latindex(a, targetextent[2])
    newextent = [getlon(a, ilons[1]) - a.res/2, getlat(a, ilats[end]) - a.res/2,
                 getlon(a, ilons[end]) + a.res/2, getlat(a, ilats[1]) + a.res/2]
    return GeoArray(a[ilons, ilats], a.res, newextent)
end

# Also write: resize_categorical(), rescale()
