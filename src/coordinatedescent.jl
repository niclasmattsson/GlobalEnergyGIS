function choosedirection(f::Function, index, deltas)
    best, i = findmin(f(index + d) for d in deltas)
    dir = deltas[i]
    return dir, best, index + dir
end

function descend(f::Function, dir, best, index)
    while true
        newindex = index + dir
        val = f(newindex)
        if val < best
            best = val
            index = newindex
        else
            return best, index
        end
    end
end

coordinate_descent(A::AbstractArray, index) =
    coordinate_descent(ndx -> get(A, ndx, typemax(eltype(A))), index)

# Custom version of findmin(A) for a matrix A which is sampled from a convex function.
# Expanded to support a function f to evaluate (to avoid precomputing all elements of A).
# https://discourse.julialang.org/t/is-there-a-findmin-a-that-uses-local-search-instead-of-visiting-all-elements-of-matrix-a/74045/4
function coordinate_descent(f::Function, index)
    ci = CartesianIndex(index)
    same = zero(ci)
    deltas = -oneunit(ci):oneunit(ci)
    while true
        dir, best, ci = choosedirection(f, ci, deltas)
        dir == same && return best, ci
        best, ci = descend(f, dir, best, ci)
    end
end

function testdescent()
    nc = Dataset(simname("mohc", 85, "v", 100, 2050, false))
    xylon = nc["lon"][:,:]
    xylat = nc["lat"][:,:]
    lons = -10.5:0.1:32
    lats = 34.5:0.1:71.5
    indices = zeros(CartesianIndex{2}, length(lons), length(lats))
    distances = zeros(length(lons), length(lats))
    lon, lat = 15, 60
    # d = @time @. (xylon - lon)^2 + (xylat - lat)^2
    # dmin = @time findmin(d)
    # f = ndx -> (xylon[ndx] - lon)^2 + (xylat[ndx] - lat)^2
    # f = ndx -> (get(xylon, ndx, Inf) - lon)^2 + (get(xylat, ndx, Inf) - lat)^2
    f = ndx -> checkbounds(Bool, xylon, ndx) ? (xylon[ndx] - lon)^2 + (xylat[ndx] - lat)^2 : Inf
    # dmin = @time coordinate_descent(f, CartesianIndex(3,3))
    dmin = @btime coordinate_descent($f, $CartesianIndex(250,310))
end

# function choosedirection(A::AbstractArray, index, deltas)
#     best, i = findmin(get(A, index + d, typemax(eltype(A))) for d in deltas)
#     dir = deltas[i]
#     return dir, best, index + dir
# end

# function descend(A::AbstractArray, dir, best, index)
#     while true
#         newindex = index + dir
#         val = get(A, newindex, typemax(eltype(A)))
#         if val < best
#             best = val
#             index = newindex
#         else
#             return best, index
#         end
#     end
# end

# function coordinate_descent(A::AbstractArray, index)
#     ci = CartesianIndex(index)
#     same = zero(ci)
#     deltas = -oneunit(ci):oneunit(ci)
#     while true
#         dir, best, ci = choosedirection(A, ci, deltas)
#         dir == same && return best, ci
#         best, ci = descend(A, dir, best, ci)
#     end
# end
