export makesolarera5, makemonthlysolarera5, clearvars_era5

# TO DO: save FDIR to disk (masked by land cover) for parabolic troughs (oriented north-south).

# Can optionally zero water cells in the landcover dataset to save a lot of disk space.
# add options to save all three CSP variables: tower, trough-NS, trough-EW.
function makesolarera5(; year=2018, land_cells_only=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")

    filename = joinpath(datafolder, "era5solar$year.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    land = imresize(JLD.load(joinpath(datafolder, "landcover.jld"), "landcover"), gridsize)

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        # create GTI and DNI variables (Global Tilted Irradiance and Direct Normal Irradiance)
        dataset_GTI = create_dataset(group, "GTI", datatype(Float32), dataspace(hours,gridsize...), chunk=(hours,16,16), blosc=3)
        dataset_DNI = create_dataset(group, "DNI", datatype(Float32), dataspace(hours,gridsize...), chunk=(hours,16,16), blosc=3)
        dataset_meanGTI = create_dataset(group, "meanGTI", datatype(Float32), dataspace(gridsize...), chunk=gridsize, blosc=3)
        dataset_meanDNI = create_dataset(group, "meanDNI", datatype(Float32), dataspace(gridsize...), chunk=gridsize, blosc=3)

        totalGTI = zeros(gridsize)
        totalDNI = zeros(gridsize)
        hour = 1

        count = 0
        for month = 1:12, monthhalf = 1:2
            if monthhalf == 1
                firstday, lastday = "01", "15"
            else
                firstday = "16"
                lastday = Dates.daysinmonth(Date("$year-$month"))
            end
            monthstr = lpad(month,2,'0')
            date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
            erafile = joinpath(downloadsfolder, "solar$year-$monthstr$firstday-$monthstr$lastday.nc")

            count += 1
            println("\nFile $count of 24:")
            println("Reading solar diffuse and direct components from $erafile...")
            ncdataset = Dataset(erafile)
            # GHI = replace(ncdataset["ssrd"][:,:,:], missing => 0.0) .* (land .> 0) ./ (3600*1000)
            # DHI = GHI - replace(ncdataset["fdir"][:,:,:], missing => 0.0) .* (land .> 0) ./ (3600*1000)
            ssrd = nomissing(ncdataset["ssrd"][:,:,:], 0.0)
            fdir = nomissing(ncdataset["fdir"][:,:,:], 0.0)
            datetime = nomissing(ncdataset["time"][:], DateTime(0))
            GTI, DNI = @time transform_solar_vars(ssrd, fdir, datetime, land, land_cells_only)

            totalGTI += sumdrop(GTI, dims=1)
            totalDNI += sumdrop(DNI, dims=1)
            len = size(GTI,1)
            println("Writing to $filename...")
            dataset_GTI[hour:hour+len-1,:,:] = GTI
            dataset_DNI[hour:hour+len-1,:,:] = DNI
            hour += len
        end
        println("\nWriting annual mean solar variables to $filename...")
        dataset_meanGTI[:,:] = totalGTI/hours
        dataset_meanDNI[:,:] = totalDNI/hours
    end
    nothing
end

function makemonthlysolarera5(; land_cells_only=true)
    years = 1979:2019
    nyears = length(years)
    nmonths = nyears*12
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")
    
    filename = joinpath(datafolder, "era5monthlysolar.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    land = imresize(JLD.load(joinpath(datafolder, "landcover.jld"), "landcover"), gridsize)

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        # create GTI and DNI variables (Global Tilted Irradiance and Direct Normal Irradiance)
        dataset_ssrd = create_dataset(group, "monthlyssrd", datatype(Float32), dataspace(nmonths,gridsize...), chunk=(nmonths,16,16), blosc=3)
        dataset_fdir = create_dataset(group, "monthlyfdir", datatype(Float32), dataspace(nmonths,gridsize...), chunk=(nmonths,16,16), blosc=3)
        dataset_annualssrd = create_dataset(group, "annualssrd", datatype(Float32), dataspace(nyears,gridsize...), chunk=(nyears,16,16), blosc=3)
        dataset_annualfdir = create_dataset(group, "annualfdir", datatype(Float32), dataspace(nyears,gridsize...), chunk=(nyears,16,16), blosc=3)
   
        erafile = in_datafolder("downloads", "monthlysolar_$(years[1])-$(years[end]).nc")

        println("Reading solar diffuse and direct components from $erafile...")
        # Permute dimensions to get hours as dimension 1 (for efficient iteration in GISwind())
        ncdataset = Dataset(erafile)
        ssrd = permutedims(nomissing(ncdataset["ssrd"][:,:,:], 0.0) .* (land .> 0), [3,1,2])
        fdir = permutedims(nomissing(ncdataset["fdir"][:,:,:], 0.0) .* (land .> 0), [3,1,2])

        println("Writing to $filename...")
        # For these monthly average insolations we skip the sun position calculations
        # made for the hourly dataset and just assign SSRD & FDIR directly.
        dataset_ssrd[:,:,:] = ssrd
        dataset_fdir[:,:,:] = fdir
        for y = 1:nyears
            dataset_annualssrd[y,:,:] = sum(ssrd[12*(y-1) .+ (1:12),:,:], dims=1) ./ 12
            dataset_annualfdir[y,:,:] = sum(fdir[12*(y-1) .+ (1:12),:,:], dims=1) ./ 12
        end
    end
    nothing
end

function transform_solar_vars(ssrd, fdir, datetime, land, land_cells_only)
    println("Calculating GTI and DNI using solar positions...")

    erares = 0.28125
    eralons = -180+erares/2:erares:180-erares/2     # longitude values (pixel center)
    eralats = 90-erares/2:-erares:-90+erares/2      # latitude values (pixel center)

    # Permute dimensions to get hours as dimension 1 (for efficient iteration in GISsolar())
    nlons, nlats, nhours = size(ssrd)
    GTI = zeros(Float32, (nhours, nlons, nlats))
    DNI = zeros(Float32, (nhours, nlons, nlats))

    albedo = 0.2                        # standard average value for ground albedo (maybe use landcover data later)
    azimuthPV = π * (eralats .< 0)      # assume equator-facing PV panels (azimuth = 0 if lat>0, π if lat<0)
    tiltPV = optimalPVtilt.(eralats)    # degrees
    # println("Check if PV tilt sign is correct in the southern hemisphere!")
    cos_tiltPVs = cosd.(tiltPV)
    sin_tiltPVs = sind.(tiltPV)

    almostzero = eps(Float32)

    for h = 1:nhours
        # Note that while the solar position calculations are instantaneous positions, ERA5 radiation variables
        # represent accumulated radiation over the hour *ending* at the indicated time. Therefore, solar positions
        # must be shifted 30 minutes BACK to correspond to the midpoint time of the ERA5 accumulations.
        # Source: ERA5 "accumulations are over the hour ending at the forecast step"
        # https://confluence.ecmwf.int//display/CKB/ERA5+data+documentation#ERA5datadocumentation-Meanratesandaccumulations
        dt = datetime[h] - Minute(30)
        TSI = solarinsolation(dt)/1000          # Total Solar Irradiance (top of atmosphere, perpendicular to sun)
        for lon = 1:nlons
            eralon = eralons[lon]
            δ, H = solarposition(dt, eralon)    # absolute solar position (declination, hour angle)
            solarpos = sines_and_cosines(δ,H)
            for lat = 1:nlats
                eralat = eralats[lat]
                land_cells_only && land[lon,lat] == 0 && continue

                zenith, azimuth = zenith_azimuth(eralat, solarpos...)     # relative solar position
                cos_zen = max(almostzero, cos(zenith))

                cos_tiltPV = cos_tiltPVs[lat]
                sin_tiltPV = sin_tiltPVs[lat]

                # ERA5 radiations are in J/m2/period, so for hourly data divide by 3600*1000 to get kW/m2
                GHI = ssrd[lon,lat,h]/(3600*1000)       # Global Horizontal Irradiance
                FDIR = fdir[lon,lat,h]/(3600*1000)
                DHI = GHI - FDIR                        # Diffuse Horizontal Irradiance

                # When the solar elevation is close to 0, both FDIR and cos(zenith) will also be near 0, and
                # calculated DNI will approach "0/0". So we'll clamp DNI to avoid artifacts.
                # That wasn't enough, so we'll add an artificial term to increase the denominator near the horizon. 
                dni = clamp(FDIR / (cos_zen + horizoncorrection(zenith)), 0, TSI)  # Direct Normal Irradiance

                # Cosine of angle of incidence
                cos_AOI = max(0, cos_zen*cos_tiltPV + sin(zenith)*sin_tiltPV*cos(azimuth-azimuthPV[lat]))
                Rb = cos_AOI/max(0.017, cos_zen)                    # cos_zen clamped to cos(89) = 0.017

                # AI is the "anisotropy index", i.e. beam radiation transmittance through atmosphere.
                # We need it to consider circumsolar diffuse radiation (the bright area of the sky near the sun). 
                AI = dni / TSI

                beamradiation = dni*cos_AOI                         # direct beam radiation from the sun
                groundradiation = GHI*albedo*(1-cos_tiltPV)/2       # diffuse reflected radiation from the ground
                skyradiation = DHI*AI*Rb + DHI*(1-AI)*(1+cos_tiltPV)/2   # diffuse radiation from the sky

                GTI[h,lon,lat] = beamradiation + groundradiation + skyradiation
                DNI[h,lon,lat] = dni          
            end
        end
    end
    return GTI, DNI
end

# Optimal tilt angle (degrees) of PV modules as a function of latitude
# Jacobson 2018, https://web.stanford.edu/group/efmh/jacobson/Articles/I/TiltAngles.pdf
# (Motivation: when the sky is overcast, a shallower tilt captures more diffuse radiation.
# Higher latitutes tend to have cloudier weather. Also, you can play with the tilt to
# trade-off between summer and winter insolation.)
function optimalPVtilt(ϕ)
    if ϕ > 0
        clamp(1.3793 + ϕ*(1.2011 + ϕ*(-0.014404 + ϕ*0.000080509)), 0, 90)
    else
        clamp(-0.41657 + ϕ*(1.4216 + ϕ*(0.024051 + ϕ*0.00021828)), -90, 0)
    end
end

# correction term to reduce horizon artifacts
horizoncorrection(zen) = 1/50 * max(0, rad2deg(zen)-85).^2

function clearvars_era5(; year=2018, datasets=["wind", "solar", "temp"])
    for dataset in datasets, month = 1:12, monthhalf = 1:2
        if monthhalf == 1
            firstday, lastday = "01", "15"
        else
            firstday = "16"
            lastday = Dates.daysinmonth(Date("$year-$month"))
        end
        monthstr = lpad(month,2,'0')
        date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
        erafile = in_datafolder("downloads", "$dataset$year-$monthstr$firstday-$monthstr$lastday.nc")
        if isfile(erafile)
            println("Deleting $erafile...")
            rm(erafile)
        end
    end
end


# SOLAR INSOLATION CALCULATIONS FOR PV
#
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
# β = PV panel tilt angle
# z = solar zenith angle = 90 - solar elevation angle
# ρ = ground albedo ≈ 0.2   (we can try estimating albedo from land cover data later)
# Rb = cos_AOI / cosz    (= I{bT}/I{b}, ratio of beam irradiance on a tilted surface and a horizontal surface)
# AZsolar is the azimuth (direction of the sun) measured from south with positive direction towards west.
# AZpv 
#
# We have:
# FDIR = DNI * cosz
# GHI = SSRD = FDIR + DHI = DNI*cosz + DHI
# cos_AOI = cosz*cosβ + sinz*sinβ*cos(AZsolar-AZpv)
#
# For solar PV, we need to estimate Global Tilted Irradiance (GTI) from GHI.
# There are three components to the radiation that hits a tilted PV panel:
# 1. direct beam radiation from the sun
# 2. diffuse radiation from the sky
# 3. diffuse reflected radiation from the ground
#
# Components 1 and 3 are easy (see any textbook or paper on solar radiation, e.g. Loutzenhiser 2007):
# 1.  FDIR * Rb = FDIR * cos_AOI / cosz = DNI * cos_AOI
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
# (Hay-Davies validation with measurements: Mubarak et al 2017)
#
# 2. DHI*AI*Rb + DHI*(1-AI)*(1+cosβ)/2, with AI = FDIR / TISR = DNI / solarconstant
#
# AI is the "anisotropy index", i.e. beam radiation transmittance through atmosphere.
# (TISR, like FDIR, is also defined as radiation hitting a horizontal plane)
#
# See also textbook by Petros Axaopoulos (chapter 4)
# http://www.labri.fr/perso/billaud/Helios2/docs/page-all-tiles.php
#
#
# SOLAR INSOLATION CALCULATIONS FOR CSP
#
# Ihe insolation per m2 on a 2-axis solar tower collector is the direct normal irradiance (DNI).
#
# For a 1-axis parabolic trough, the axis is usually (but not always) oriented north-south, and the rotation
# of the trough follows the sun east-west. The collector captures the same insolation as a horizontal surface.
# A trough with its axis oriented east-west will track the elevation of the sun, but not its azimuth.
#
# North-south oriented axis: The insolation per m2 of collector is FDIR = DNI * cos(zenith).
# East-west oriented axis: The insolation per m2 of collector is DNI * cos(azimuth).

# ADJUST LAND USE BY LATITUDE!!!!!
# WE DON'T CORRECT EFFICIENCY FOR TEMPERATURES (or wind speed!) see Jacobson.
