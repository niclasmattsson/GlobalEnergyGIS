export solarposition, zenith, zenith_azimuth, sines_and_cosines

# Calculate the position of the sun as (zenith,azimuth) in degrees for a location
# on earth given as (latitude,longitude) in degrees. The DateTime should be given
# as Universal Time, e.g. UTC. 
function solarposition_wiki(datetime::DateTime, latitude, longitude)
	# Position of the Sun as seen from Earth.
	# Ported from Matlab, original helpstring below.

	#   This is the most basic algorithm. It is documented in Seinfeld &
	#   Pandis, Duffie & Beckman and Wikipedia.
	#
	# [ANGLES,PROJECTION] = SOLARPOSITION(DATE,TIME,LATITUDE,LONGITUDE,TIME_ZONE)
	# returns ZENITH & AZIMUTH for all DATE & TIME pairs at LATITUDE, LONGITUDE.
	# ANGLES = [ZENITH,AZIMUTH] and PROJECTION = [PHI_X, PHI_Y]
	# PHI_X is projection on x-z plane & PHI_Y is projection on y-z plane.
	# DATETIME can be string, vector [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS],
	#   cellstring or matrix N x [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS] for N
	#   times.
	# LATITUDE [degrees] and LONGITUDE [degrees] are the coordinates of the site.
	# TIME_ZONE [hours] of the site.
	# ROTATION [degrees] clockwise rotation of system relative to north.
	# DST [logical] flag for daylight savings time, typ. from March to November
	#   in the northern hemisphere.
	#
	# References:
	# http://en.wikipedia.org/wiki/Solar_azimuth_angle
	# http://en.wikipedia.org/wiki/Solar_elevation_angle
	#
	# Mark A. Mikofski
	# Copyright (c) 2013



	# Equation of Time
	# "Alternative calculation" on wikipedia:
	# https://en.wikipedia.org/wiki/Equation_of_time#Calculating_the_equation_of_time
	# Supposedly accurate to 6 seconds (0.1 minutes)

	year = Dates.year(datetime)
	# date in days starting at zero on 1 January (the subtraction produces a difference in milliseconds)
	d_ms = datetime - DateTime(year,1,1)
	d = d_ms.value/1000/3600/24

	# Earth's mean angular orbital velocity [degrees/day]
	w = 360/365.24
	# angle the earth moves on its orbit at its average speed from the December solstice
	# 10 is the approximate number of days from the December solstice to January 1
	a = w * (d+10)
	# The angle the Earth moves from the solstice to date D, including a first-order correction
	# for the Earth's orbital eccentricity, 0.0167. The number 2 is the number of days from 1 January
	# to the date of the Earth's perihelion.
	b = a + (360/pi*0.0167)*sind(w*(d-2))
	# C is the difference between the angles moved at mean speed, and at the corrected speed projected
	# onto the equatorial plane, and divided by 180 to get the difference in "half turns". The value 23.44°
	# is the obliquity (tilt) of the Earth's axis. The subtraction gives the conventional sign to the equation of time. 
	c = (a - atand(tand(b)/cosd(23.44))) / 180
	# For any given value of x, arctan x (sometimes written as tan−1 x) has multiple values, differing from each other
	# by integer numbers of half turns. This may cause C to be wrong by an integer number of half turns. The excess
	# half turns are removed in the next step of the calculation to give the equation of time:
	ET = 720*(c - round(c))

	# approximate solar time [hours]
	solarTime = 24*(d - floor(d)) + longitude*24/360 + ET/60
	t_h = 15*(solarTime - 12) # [degrees] hour angle

	# declination [degrees]
	# accurate to 0.2 degrees, see https://en.wikipedia.org/wiki/Position_of_the_Sun#Calculations
	delta = -asind(sind(23.44)*cosd(b))
	zenith = acosd(sind(latitude)*sind(delta) + cosd(latitude)*cosd(delta)*cosd(t_h)) # [degrees] zenith

	# azimuth [0, 180], absolute value measured from due south, so east = west = 90,
	# south = 0, north = 180
	cos_phi = (cosd(zenith)*sind(latitude) - sind(delta)) ./ (sind(zenith)*cosd(latitude)) # cosine(azimuth)
	phi_south = acosd(clamp(cos_phi,-1,1))
	#phi_south = acosd(cos_phi);

	# azimuth [0, 360], measured clockwise from due north,
	# so east = 90, south = 180, and west = 270 degrees
	azimuth = 180 + sign(t_h)*phi_south # Shift domain to 0-360 deg

	#angles = [theta, phi]; # [degrees] zenith, azimuth
	return zenith, azimuth, solarTime, delta
end

function solarposition_wiki_faster(datetime::DateTime, latitude, longitude)
	year = Dates.year(datetime)
	d_ms = datetime - DateTime(year,1,1)
	d = d_ms.value/1000/3600/24
	w = 360/365.24
	a = w * (d+10)
	b = a + (360/pi*0.0167)*sind(w*(d-2))
	c = (a - atand(tand(b)/cosd(23.44))) / 180
	ET = 720*(c - round(c))
	solarTime = 24*(d - floor(d)) + longitude*24/360 + ET/60
	t_h = 15*(solarTime - 12) # [degrees] hour angle
	delta = -asind(sind(23.44)*cosd(b))
	slat, clat = sind(latitude), cosd(latitude)
	sdel, cdel = sind(delta), cosd(delta)
	czen = slat*sdel + clat*cdel*cosd(t_h)
	zenith = acosd(czen) # [degrees] zenith

	# azimuth [0, 180], absolute value measured from due south, so east = west = 90,
	# south = 0, north = 180
	cos_phi = (czen*slat - sdel) / (sind(zenith)*clat) # cosine(azimuth)
	phi_south = acosd(clamp(cos_phi,-1,1))
	#phi_south = acosd(cos_phi);

	# azimuth [0, 360], measured clockwise from due north,
	# so east = 90, south = 180, and west = 270 degrees
	azimuth = 180 + sign(t_h)*phi_south # Shift domain to 0-360 deg

	#angles = [theta, phi]; # [degrees] zenith, azimuth
	return zenith, azimuth, solarTime, delta
end



# Splits the DateTime into (year, month, day) (as integers) and the hour (as a floating point)
function splitdatetime(datetime::DateTime)
	# date = floor(datetime, Dates.Day)
	y, m, d = Dates.yearmonthday(datetime)
	return y, m, d, (datetime - DateTime(y,m,d)).value/1000/3600
end

# time correction term, difference between Terrestrial Time and Universal Time, see Grena (2012) section 2.1.
# Conversion to TT is necessary since UT follows the rotation of the earth, which is not uniform.
# t [days] days from January 1, 2060 (in UT)
# returns Δτ [seconds]
differenceTTUT(t) = 96.4 + 0.00158*t

# Calculate time variable needed for Grena (2012) algorithms, see section 3.1.
# The time t is the input DateTime expressed as days from January 1, 2060 (in UT).
# returns t [days]
calctime(datetime::DateTime) = (datetime-DateTime(2060,1,1)).value/1000/3600/24

# Same result as above, but this is the function in Grena's paper. The Julia version is faster.
function calctime_grena(datetime::DateTime)
	y, m, d, h = splitdatetime(datetime)
	if m <= 2
		m += 12
		y -= 1
	end
	return trunc(365.25*(y-2000)) + trunc(30.6001*(m+1)) - trunc(0.01*y) + d + 0.0416667*h - 21958
end

# Full version of algorithm 1 of Grena (2012), see section 3.2.
# Calculate the position of the sun as (zenith,azimuth) in radians for a location
# on earth given as (latitude,longitude) in degrees. The DateTime should be given
# as Universal Time, e.g. UTC. 
function solarposition_grena1(datetime::DateTime, latitude, longitude)
	ω = 0.017202786 # [1/day]
	t = calctime(datetime)
	t_e = t + 1.1574e-5 * differenceTTUT(t)
	s1, c1 = sin(ω*t_e), cos(ω*t_e)
	s2, c2 = 2*s1*c1, (c1+s1)*(c1-s1)
	α = mod(-1.38880 + 1.72027920e-2*t_e + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2 + 1.525e-2*c2, 2π)
	δ = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2
	H = mod(1.75283 + 6.3003881*t + deg2rad(longitude) - α + π, 2π) - π
	sϕ, cϕ = sind(latitude), cosd(latitude)
	sδ, cδ = sin(δ), cos(δ)
	sH, cH = sin(H), cos(H)
	se0 = sϕ*sδ + cϕ*cδ*cH
	ep = asin(se0) - 4.26e-5*sqrt(1-se0^2)  #cos(e0)
	azimuth = atan(sH, cH*sϕ - sδ*cϕ/cδ)	#azimuth = 0 towards south and positive direction towards west
	zenith = π/2 - ep
	return zenith, azimuth				# radians
	# return zenith, azimuth, δ, H		# radians
	# return rad2deg(zenith), rad2deg(mod(azimuth+π,2π)), rad2deg(H)/15+12, rad2deg(δ) # degrees, same azimuth origin as solarposition0
end

# First half of algorithm 1 of Grena (2012), see section 3.2.
# Calculate the position of the sun as (declination, hour angle) in radians for a location
# on earth given by its longitude in degrees. The DateTime should be given
# as Universal Time, e.g. UTC. 
function solarposition(datetime::DateTime, longitude)
	ω = 0.017202786 # [1/day]
	t = calctime(datetime)
	t_e = t + 1.1574e-5 * differenceTTUT(t)
	s1, c1 = sin(ω*t_e), cos(ω*t_e)
	s2, c2 = 2*s1*c1, (c1+s1)*(c1-s1)
	α = mod(-1.38880 + 1.72027920e-2*t_e + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2 + 1.525e-2*c2, 2π)
	δ = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2
	H = mod(1.75283 + 6.3003881*t + deg2rad(longitude) - α + π, 2π) - π
	return δ, H
end

# Full version of algorithm 1 of Grena (2012), see section 3.2.
# Calculate the position of the sun as (zenith, azimuth) in degrees for a location
# on earth given as (latitude,longitude) in degrees. The DateTime should be given
# as Universal Time, e.g. UTC. 
function solarposition(datetime::DateTime, latitude, longitude)
	δ, H = solarposition(datetime, longitude)
	zenith, azimuth = zenith_azimuth(latitude, sines_and_cosines(δ, H)...)
	return rad2deg(zenith), rad2deg(azimuth)
end

sines_and_cosines(δ, H) = sin(δ), cos(δ), sin(H), cos(H)

function zenith_azimuth(latitude, sδ, cδ, sH, cH)
	sϕ, cϕ = sind(latitude), cosd(latitude)
	se0 = sϕ*sδ + cϕ*cδ*cH
	azimuth = atan(sH, cH*sϕ - sδ*cϕ/cδ)	# azimuth = 0 towards south and positive direction towards west
	zenith = π/2 - asin(se0) + 4.26e-5*sqrt(1-se0^2)  #cos(e0)
	return zenith, azimuth
end

function zenith(latitude, sδ, cδ, sH, cH)
	sϕ, cϕ = sind(latitude), cosd(latitude)
	se0 = sϕ*sδ + cϕ*cδ*cH
	return π/2 - asin(se0) + 4.26e-5*sqrt(1-se0^2)  # sqrt(1-se0^2) = cos(e0), but slightly faster
end

solarconstant(dt::DateTime) = 1361 * (1 + 0.034*cos(2*π*(dt - DateTime(year(dt))).value/1000/3600/24/365.25))

# dt = DateTime(2019,4,15,13,46)
# δ, H = solarposition(dt, 12)	# return rad2deg(zenith), rad2deg(mod(azimuth+π,2π)), rad2deg(H)/15+12, rad2deg(δ)

# pos = sines_and_cosines(δ, H)
# zen = zenith(58, pos...)




# ERA5 radiations are those accumulated during one hour onto a square meter of a flat horizontal plane.
# This holds for both SSRD and FDIR (see Hogan, "Radiation Quantities in the ECMWF model and MARS").
# ERA5 vars:  SSRD = Surface Solar Radiation Downwards, FDIR = DIRect solar radiation at the surface.
#
# We will use the following abbreviations and identifiers:
#
# GHI = Global Horizontal Irradiance (here "global" means direct + diffuse)
# DHI = Diffuse Horizontal Irradiance
# DNI = Direct Normal Irradiance
#
# AOI = angle of incidence of sun on panel
# β = panel tilt angle
# z = solar zenith angle = 90 - solar elevation angle
# ρ = ground albedo ≈ 0.2	(we can try estimating albedo from land cover data later)
# Rb = cosAOI / cosz    (projection coefficient of horizontal beam radiation on tilted surface)
#
# We have:
# FDIR = DNI * cosz
# GHI = SSRD = FDIR + DHI = DNI*cosz + DHI
# cosAOI = cosz*cosβ + sinz*sinβ*cos(AZsolar-AZpv)
#
# For solar PV, we need to estimate Global Tilted Irradiance (GTI) from GHI.
# There are three components to the radiation that hits a tilted PV panel:
# 1. direct beam radiation from the sun
# 2. diffuse radiation from the sky
# 3. diffuse reflected radiation from the ground
#
# Components 1 and 3 are easy (see any textbook or paper on solar radiation, e.g. Lave 2015):
# 1.  FDIR * Rb = FDIR * cosAOI / cosz = DNI * cosAOI
# 3.  GHI * ρ * (1 - cosβ)/2
#
# For an evaluation of approaches on estimation of diffuse sky radiation, see Loutzenhiser 2007.
# For now we will use the simplest approach and assume the sky is isotropic (i.e. diffuse radiation
# is uniformly distributed over the sky dome). Then we have:
#
# 2. DHI * (1 + cosβ)/2
#
# But the sky around the sun is significantly brighter ("circumsolar diffuse radiation"). Later,
# if we also download the TISR variable from ERA5 (top of atmosphere incident solar radiation),
# we can use the Hay and Davies model from 1980 to consider circumsolar diffuse radiation:
#
# 2. DHI*AI*Rb + DHI*(1-AI)*(1+cosβ)/2, with AI = FDIR / TISR = DNI / solarconstant
#
# AI is the "anisotropy index", i.e. beam radiation transmittance through atmosphere.
# (TISR, like FDIR, is also defined as radiation hitting a horizontal plane)
#
# Relations from textbook by Petros Axaopoulos (chapter 4)
# http://www.labri.fr/perso/billaud/Helios2/docs/page-all-tiles.php
#
# Similarly for CSP, the insolation per m2 on a 2-axis solar tower collector is approximately FDIR/cos(zenith).
# For a 1-axis parabolic trough, the axis is usually (but not always) oriented north-south, and the rotation
# of the trough follows the sun east-west.

# North-south oriented axis: The insolation per m2 of collector is FDIR.
# East-west oriented axis: The insolation per m2 of collector is FDIR/cos(zenith) * cos(azimuth).

# ADJUST LAND USE BY LATITUDE!!!!!
# GET OPTIMUM PV TILT ANGLES FROM Jacobson 2018:
# https://web.stanford.edu/group/efmh/jacobson/Articles/I/TiltAngles.pdf
# http://www.ftexploring.com/solar-energy/tilt-angle1.htm
# http://www.ftexploring.com/solar-energy/tilt-angle2.htm
# http://www.ftexploring.com/solar-energy/tilt-angle3.htm
# WE DON'T CORRECT EFFICIENCY FOR TEMPERATURES (or wind speed!) see Jacobson.

# Here zenith is the solar zenith angle = π/2 - elevation. Azimuth is the direction of the sun measured from
# south with positive direction towards west.
#
# Note that when the solar elevation is close to 0, both FDIR and cos(zenith) will also be near 0, and all
# calculated solar insolations will have a "0/0" quotient. If the underlying 

