function bbox2ranges(bbox, rasterdensity)
    latindexes, lonindexes = bbox2indexes(bbox, rasterdensity)
    latindex = latindexes[1] : latindexes[2]
    lonindex = lonindexes[1] : lonindexes[2]
    return latindex, lonindex
end

function bbox2indexes(bbox, rasterdensity)
    latindexes = Int.(reverse(rasterdensity*(90 .- bbox[:,1]), dims=1) + [1 0]')
    lonindexes = Int.(rasterdensity*(bbox[:,2] .- (-180)) + [1 0]')
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