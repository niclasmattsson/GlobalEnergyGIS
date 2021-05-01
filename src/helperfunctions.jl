# bbox is a "bounding box" representing outer limits of latitude & longitude: [bottomlat leftlon; toplat rightlon]
# rasterdensity is number of pixels per degree, so rasterdensity = 1/res with res = raster resolution measured in degrees
# TO DO: make these consistent with lon, lat convention in rest of code, use resolution instead of rasterdensity
function bbox2ranges(bbox, rasterdensity)
    latindexes, lonindexes = bbox2indexes(bbox, rasterdensity)
    latindex = latindexes[1] : latindexes[2]
    lonindex = lonindexes[1] : lonindexes[2]
    return latindex, lonindex
end

function bbox2indexes(bbox, rasterdensity)
    latindexes = round.(Int, reverse(rasterdensity*(90 .- bbox[:,1]), dims=1) + [1 0]')
    lonindexes = round.(Int, rasterdensity*(bbox[:,2] .- (-180)) + [1 0]')
    return latindexes, lonindexes
end

function roundbbox(bbox, rasterdensity)
    newbbox = zeros(2,2)
    newbbox[1,:] = floor.(bbox[1,:]*rasterdensity)/rasterdensity
    newbbox[2,:] = ceil.(bbox[2,:]*rasterdensity)/rasterdensity
    return newbbox
end

filterrange(range, lowhigh) = range[(range .>= lowhigh[1]) .& (range .<= lowhigh[2])]

function getlats(bbox, rasterdensity, halfshift)
	bottomlat, toplat = bbox[:,1]
	res = 1/rasterdensity
	shift = halfshift ? res/2 : 0.0
	return bottomlat + shift : res : toplat - shift
	# return toplat - shift : -res : bottomlat + shift
end

function getlons(bbox, rasterdensity, halfshift)
	leftlon, rightlon = bbox[:,2]
	res = 1/rasterdensity
	shift = halfshift ? res/2 : 0.0
	return leftlon + shift : res : rightlon - shift
end

# Boolean function, true for a leap year (for Gregorian calendar, assumes year >= 1583).
# The Gregorian calendar has 97 leap years every 400 years: 
#     Every year divisible by 4 is a leap year. 
#     However, every year divisible by 100 is not a leap year. 
#     However, every year divisible by 400 is a leap year after all. 
# So, 1700, 1800, 1900, 2100, and 2200 are not leap years, but 1600, 2000, and 2400 are leap years.
leapyear(year) = divisible(year,4) && (!divisible(year,100) || divisible(year,400))
divisible(x,y) = x % y == 0


function meshgrid2d(x,y)
	xx = repeat(x', length(y), 1)
	yy = repeat(y, 1, length(x))
	return xx, yy
end

# Borrowed from Matlab, see fspecial('disk', rad)
function diskfilterkernel(rad)
	crad  = ceil(Int, rad-0.5)
	x,y = meshgrid2d(-crad:crad,-crad:crad)
	maxxy = max.(abs.(x),abs.(y))
	minxy = min.(abs.(x),abs.(y))
	m1 = (rad^2 .< (maxxy .+ 0.5).^2 + (minxy .- 0.5).^2) .* (minxy .- 0.5) +
	  (rad^2 .>= (maxxy .+ 0.5).^2 + (minxy .- 0.5).^2) .* sqrt.(max.(0, rad^2 .- (maxxy .+ 0.5).^2))
	m2 = (rad^2 .> (maxxy .- 0.5).^2 + (minxy .+ 0.5).^2) .* (minxy .+ 0.5) + 
	  (rad^2 .<= (maxxy .- 0.5).^2 + (minxy .+ 0.5).^2) .* sqrt.(max.(rad^2 .- (maxxy .- 0.5).^2))
	sgrid = (rad^2 * (0.5*(asin.(m2/rad) .- asin.(m1/rad)) + 0.25*(sin.(2*asin.(m2/rad)) - sin.(2*asin.(m1/rad)))) +
	     		- (maxxy .- 0.5).*(m2-m1) + (m1 .- minxy .+ 0.5)) .*
	      ((((rad^2 .< (maxxy .+ 0.5).^2 + (minxy .+ 0.5).^2) .&
	      	(rad^2 .> (maxxy .- 0.5).^2 + (minxy .- 0.5).^2)) .|
	      	((minxy .== 0) .& (maxxy .- 0.5 .< rad) .& (maxxy .+ 0.5 .>= rad))))
	sgrid = sgrid .+ ((maxxy .+ 0.5).^2 + (minxy .+ 0.5).^2 .< rad^2)
	sgrid[crad+1,crad+1] = min.(pi*rad^2, pi/2)
	if ((crad>0) && (rad > crad-0.5) && (rad^2 < (crad-0.5)^2+0.25)) 
		m1a  = sqrt(rad^2 - (crad - 0.5).^2)
		m1n = m1a/rad
		sg0 = 2*(rad^2*(0.5*asin(m1n) + 0.25*sin(2*asin(m1n)))-m1a*(crad-0.5))
		sgrid[2*crad+1,crad+1] = sg0
		sgrid[crad+1,2*crad+1] = sg0
		sgrid[crad+1,1]        = sg0
		sgrid[1,crad+1]        = sg0
		sgrid[2*crad,crad+1]   = sgrid[2*crad,crad+1] - sg0
		sgrid[crad+1,2*crad]   = sgrid[crad+1,2*crad] - sg0
		sgrid[crad+1,2]        = sgrid[crad+1,2]      - sg0
		sgrid[2,crad+1]        = sgrid[2,crad+1]      - sg0
	end
	sgrid[crad+1,crad+1] = min(sgrid[crad+1,crad+1], 1)
	h = centered(sgrid/sum(sgrid))
end


meandrop(x; dims=dims) = dropdims(mean(x, dims=dims), dims=dims)
sumdrop(x; dims=dims) = dropdims(sum(x, dims=dims), dims=dims)

shifthalfcell(data, kernel=centered([.25 .25; .25 .25])) = imfilter(data, kernel)

function shiftallcells!(windCF)
    yearlength = size(windCF, 1)
    updateprogress = Progress(yearlength, 1)
    kernel = centered([.25 .25; .25 .25])
    for i=1:yearlength
        windCF[i,:,:] = shifthalfcell(windCF[i,:,:], kernel)
        next!(updateprogress)
    end
end

function max!(x::AbstractArray, val)
    @inbounds for i in eachindex(x)
        x[i] = max(x[i], val)
    end
	x
end

selfmap!(f,x) = map!(f,x,x)

# Apply a function fn to the data matrix by chunks to avoid memory issues
function gridsplit(data::AbstractArray, fn::Function, resulttype::DataType; nmax=9000, overlap=250)
    rows, cols = size(data)
    nparts = length(1:nmax:rows) * length(1:nmax:cols)
    nparts > 1 && println("\nSplitting data matrix into $nparts parts to avoid memory issues...")
    result = zeros(resulttype, size(data))
    part = 0
    for r = 1:nmax:rows, c = 1:nmax:cols
        part += 1
        nparts > 1 && println("\nPart $part/$nparts:")
        rowrange_in = max(r-overlap, 1):min(r-1+nmax+overlap, rows)
        colrange_in = max(c-overlap, 1):min(c-1+nmax+overlap, cols)
        rowrange_out = r:min(r-1+nmax, rows)
        colrange_out = c:min(c-1+nmax, cols)
        rowdelta = rowrange_out[1] - rowrange_in[1] + 1
        coldelta = colrange_out[1] - colrange_in[1] + 1
        rowrange_delta = rowdelta:(rowdelta + length(rowrange_out) - 1)
        colrange_delta = coldelta:(coldelta + length(colrange_out) - 1)
        result[rowrange_out,colrange_out] = fn(data[rowrange_in,colrange_in])[rowrange_delta,colrange_delta]
    end
    return result
end

# function gridsurface(x, y; color=color, nmax=9000, args...)
#     rows, cols = size(color)
#     rowchunks, colchunks = ceil(Int, rows/nmax), ceil(Int, cols/nmax)
#     rowsplit, colsplit = ceil(Int, rows/rowchunks), ceil(Int, cols/colchunks)
#     nparts = rowchunks * colchunks
#     nparts > 1 && println("\nSplitting surface plat into $nparts parts to avoid Makie GPU limits...")
#     part = 0
#     scene = Scene()
#     for r = 1:rowsplit:rows, c = 1:colsplit:cols
#         part += 1
#         nparts > 1 && println("\nPart $part/$nparts:")
#         rowrange = r:min(r-1+rowsplit, rows)
#         colrange = c:min(c-1+colsplit, cols)
#         surface!(x[rowrange,colrange], y[rowrange,colrange]; color=color[rowrange,colrange], args...)
#     end
#     return scene
# end

row2lon(row::Int, res) = (row - 0.5) * res - 180
col2lat(col::Int, res) = 90 - (col - 0.5) * res
lon2row(lon::Real, res) = floor(Int, mod(180+lon, 360) / res) + 1
lat2col(lat::Real, res) = floor(Int, min(180-res/2, 90-lat) / res) + 1

function rowcol2lonlat(rowcol::Tuple{Int,Int}, res)
    row, col = rowcol
    @assert 1 <= row <= 360/res
    @assert 1 <= col <= 180/res
    return row2lon(row,res), col2lat(col,res)
end

function lonlat2rowcol(lonlat::Tuple{<:Real,<:Real}, res)
    lon, lat = lonlat
    @assert -90 <= lat <= 90
    return lon2row(lon,res), lat2col(lat, res)
end

# convert indexes of datasets with 0.01 degree resolution to indexes of ERA5 resolution (0.28125 degrees)
function eraranges(lonrange, latrange, res, erares)
    topleft = rowcol2lonlat((lonrange[1],latrange[1]), res)
    bottomright = rowcol2lonlat((lonrange[end],latrange[end]), res)
    row1, col1 = lonlat2rowcol(topleft, erares)
    row2, col2 = lonlat2rowcol(bottomright, erares)
    nlons = round(Int, 360/erares)
    eralonranges = row2 >= row1 ? [row1:row2] : [row1:nlons, 1:row2] 
    eralatrange = col1:col2
    return eralonranges, eralatrange
end

function eraindexlookup(lons, lats, eralonranges, eralatrange)
    eralonrange = length(eralonranges) == 1 ? eralonranges[1] : [eralonranges[1]; eralonranges[2]]
    lonindexlookup = zeros(Int,1280)
    latindexlookup = zeros(Int,640)
    for (i,r) in enumerate(eralonrange)
        lonindexlookup[r] = i
    end
    for (i,r) in enumerate(eralatrange)
        latindexlookup[r] = i
    end
    eralonindexes = lonindexlookup[lon2row.(lons,0.28125)]
    eralatindexes = latindexlookup[lat2col.(lats,0.28125)]
    return eralonindexes, eralatindexes
end

# Automatically detect the "bounding box" of nonzero data in a matrix.
# Returns a tuple of indexes of the box containing data, (lon,lat) or (row,col).
function getbboxranges(regions::AbstractMatrix, padding::Int=0)
    data = (regions .> 0) .& (regions .!= NOREGION)
    lonrange = dataindexes_lon(vec(any(data, dims=2)), padding)     # all longitudes with region data
    latrange = dataindexes_lat(vec(any(data, dims=1)), padding)     # all latitudes with region data
    return lonrange, latrange
end

# Given a vector of Bool indicating nonzero data along latitudes, return a range
# of indexes indicating the area where there is data.
# Optionally add some padding elements to the data area (to ensure that there is
# some offshore territory around the regions).
function dataindexes_lat(latdata::AbstractVector, padding::Int=0)
    len = length(latdata)
    first, last = extrema(findall(latdata))         # first & last indexes of the data region
    return max(1, first-padding):min(len, last+padding)   # lat indexes of region elements
end

# Given a vector of Bool indicating nonzero data along longitudes, return a vector
# of indexes indicating the area where there is data. Do this by eliminating the
# longest contiguous sequence of zero data, considering wraparound (at lon=180).
# Optionally add some padding elements to the data area (to ensure that there is
# some offshore territory around the regions).
function dataindexes_lon(londata::AbstractVector, padding::Int=0)
    len = length(londata)
    seq = longest_circular_sequence(londata, false) # first & last indexes of the longest contiguous sequence of NONregion elements  
    first, last = seq[2]+1, seq[1]-1                # the rest are region elements (including "holes" in the sequence)
    last = last >= first ? last : last + len        # make sure last > first so the loop below works
    return [mod1(i,len) for i = first-padding:last+padding]   # lon indexes of region elements
end

# Return the start and end indexes of the longest contiguous sequence of v such that v[i] == x,
# allowing the sequence to wrap around the end. 
function longest_circular_sequence(v::AbstractVector, x)
    longest = Int[]
    current = Int[]
    haswrapped = false
    i = 1
    len = length(v)
    sequencelength(seq) = isempty(seq) ? 0 : mod(seq[2]-seq[1], len) + 1
    while true
        if v[i] == x
            if isempty(current)
                current = [i,i]
            else
                current[2] = i
            end
        else
            longest = sequencelength(current) > sequencelength(longest) ? current : longest
            haswrapped && return longest
            current = Int[]
        end
        if !haswrapped && i==len
            haswrapped = true
            i = 1
            sequencelength(current) == len && return current
        else
            i += 1
        end
    end
end

function coordmap(len::Int, croprange)
    coords = zeros(Int, len)
    indexes = (1:len)[croprange]
    for (i,n) = enumerate(indexes)
        coords[n] = i
    end
    return coords
end

# area of a grid cell in km2 for a given latitude (in degrees) and raster resolution
rastercellarea(lat, res) = cosd(lat) * (2*6371*π/(360/res))^2

function eralonlat(options, lonrange, latrange)
    @unpack res, erares = options

    lons = (-180+res/2:res:180-res/2)[lonrange]         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    eralonrange = length(eralonranges) == 1 ? eralonranges[1] : [eralonranges[1]; eralonranges[2]]
    eralons = (-180+erares/2:erares:180-erares/2)[eralonrange]     # longitude values (pixel center)
    eralats = (90-erares/2:-erares:-90+erares/2)[eralatrange]      # latitude values (pixel center)
    
    lonmap = coordmap(round(Int, 360/res), lonrange)      # map full longitude indexes to cropped indexes
    latmap = coordmap(round(Int, 180/res), latrange)      # ditto latitude

    return eralons, eralats, lonmap, latmap, cellarea
end

function uncrop(croppedarray, lonrange, latrange, res)
    nlon, nlat = round(Int, 360/res), round(Int, 180/res)
    full = zeros(eltype(croppedarray), nlon, nlat)
    full[lonrange, latrange] = croppedarray
    return full
end

function drawmap(mapdata)
    skip = ceil(Int, maximum(size(mapdata))/3000)
    mirrormap = reverse(mapdata[1:skip:end,1:skip:end], dims=2)
    # display(heatmap(mirrormap, size=(1200, 900), c=[cgrad(:viridis)[x] for x in 0.0:0.2:1.0]))
    # display(heatmap(mirrormap, size=(1200, 900), c=cgrad(:viridis, [0, 0.5, 1])))
    heatmap(mirrormap)
    # display(countmap(mirrormap[:]))
end

# I can't believe Makie can't do this more conveniently.
function plotlines(x,y)
    scene = Scene()
    for i = 1:size(y,2)
        lines!(scene, x, y[:,i], color=ColorSchemes.Set1_9[i])
    end
    scene
end

function drawregionmap(regionname)
    plotly()
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(regionname)
    skip = ceil(Int, maximum(size(mapdata))/3000)
    reg = reverse(regions[1:skip:end,1:skip:end]', dims=1)
    display(heatmap(reg .+ (reg.>0).*20, size=(1200, 900)))
    reg = reverse(offshoreregions[1:4:end,1:4:end]', dims=1)
    display(heatmap(reg .+ (reg.>0).*20, size=(1200, 900)))
end

function resize_categorical(regions, regionlist, lonrange, latrange; skipNOREGION=false)
    res = 0.01          # resolution of auxiliary datasets [degrees per pixel]
    erares = 0.28125    # resolution of ERA5 datasets [degrees per pixel]

    eralons, eralats, lonmap, latmap, cellarea = eralonlat(Dict(:res=>res, :erares=>erares), lonrange, latrange)
    numreg = length(regionlist)
    nlons, nlats = length(eralons), length(eralats)

    RegType = eltype(regions)
    smallregions = zeros(RegType, nlons, nlats)
    count0 = Dict(i => 0 for i = 0:numreg)
    count0[NOREGION] = 0

    for (j, eralat) in enumerate(eralats)
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for (i, eralon) in enumerate(eralons)
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]
            count = copy(count0)
            foundnonzero = false
            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                skipNOREGION && (reg == 0 || reg == NOREGION) && continue
                count[reg] += 1
                foundnonzero = true
            end
            smallregions[i,j] = foundnonzero ? findmax(count)[2] : 0
        end
    end

    return smallregions
end

function matlab2elin(; gisregion="Europe8", year=2018)
    filenamesuffix = ""
    _, _, regionlist, _, _ = loadregions(gisregion)
    
    # CF_pvrooftop, capacity_pvrooftop
    region = string.(regionlist)
    tech = ["PVplantA", "PVplantB", "PVroof", "CSPA", "CSPB", "WindonshoreA", "WindonshoreB", "Windoffshore"]
    classname = ["PVPA", "PVPB", "PVR", "CSPA", "CSPB", "WONA", "WONB", "WOFF"]
    capvar = ["capacity_pvplantA", "capacity_pvplantB", "capacity_pvrooftop", "capacity_cspplantA", "capacity_cspplantB",
                "capacity_onshoreA", "capacity_onshoreB", "capacity_offshore"]
    cfvar = ["CFtime_pvplantA", "CFtime_pvplantB", "CFtime_pvrooftop", "CFtime_cspplantA", "CFtime_cspplantB",
                "CFtime_windonshoreA", "CFtime_windonshoreB", "CFtime_windoffshore"]

    winddata = matread(in_datafolder("output", "GISdata_wind$(year)_$gisregion$filenamesuffix.mat"))
    solardata = matread(in_datafolder("output", "GISdata_solar$(year)_$gisregion$filenamesuffix.mat"))
    data = merge(winddata, solardata)

    # read number of classes from wind & solar GIS output 
    nwindclasses = [size(winddata[varname], 2) for varname in capvar[6:8]] 
    nsolarclasses = [size(solardata[varname], 2) for varname in capvar[1:5]] 
    nclasses = [nsolarclasses; nwindclasses]

    for t = 1:length(tech)
        open(in_datafolder("output", "capacity_$(tech[t]).inc"), "w") do f
            for (r,reg) in enumerate(region)
                for c = 1:nclasses[t]
                    val = data[capvar[t]][r,c]
                    !isnan(val) && val > 0 && @printf(f, "%-3s . %s%-2d %12.6f\n", reg, classname[t], c, val)
                end
            end
        end
        open(in_datafolder("output", "cf_$(tech[t]).inc"), "w") do f
            for (r,reg) in enumerate(region)
                for c = 1:nclasses[t]
                    for h = 1:8760
                        val = data[cfvar[t]][h,r,c]
                        !isnan(val) && val > 0 && @printf(f, "%-3s . %s%-2d . h%04d %10.6f\n", reg, classname[t], c, h, val)
                    end
                end
            end
        end
    end
end

function dms2deg(dms)
    sg = sign(dms)
    dms = abs(dms)
    d = sg*floor(dms/10000)
    dms = dms-10000*d
    m = floor(dms/100)
    return d + m/60 + (dms-100*m)/3600
end

# Capital recovery factor for investments
CRF(r, T) = r / (1 - 1/(1+r)^T)

in_datafolder(names...) = joinpath(getconfig("datafolder"), names...)
