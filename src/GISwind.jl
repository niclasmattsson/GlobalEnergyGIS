function GISwind()

    GISREGION = "Europe8"     # 'global', 'europe8/10/12/15/17', 'china6', 'eurochine14', 'eurasia38/21' or 'mena'

    ONSHORE_DENSITY = 5        # about 30# of existing farms have at least 5 W/m2, will become more common
    OFFSHORE_DENSITY = 8       # varies a lot in existing parks (4-18 W/m2)
    AREA_ONSHORE = .05         # area available for onshore wind power after the masks have been applied
    AREA_OFFSHORE = .33        # area available for offshore wind power after the masks have been applied

    #AVAILABILITY = 1          # now included in the turbine power curve (all code for this has been removed)

    # For reference:  10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2

    # (The datasets have a 5 minute resolution, i.e. each pixel is approximately 9.3 km at the equator.)
    ELEC_ACCESS_PIXELS = 16    # Access to electricity < 16*9.3 km = 150 km  (for wind classes of category B and offshore)
    PERSONS_PER_KM2 = 75       # Not too crowded, max x persons/km2
    MAX_DEPTH = 40             # max depth for offshore wind
    MIN_PIXELS_TO_SHORE = 1    # minimum distance to shore for offshore wind > 1*9.3 km
    EXCLUDE_LANDTYPES = [0, 11, 13]         # Not water, wetlands or urban. See table below. (But wind power over forests is coming ....)
    PROTECTED_CODES = [3, 4, 5, 6, 7, 8]    # These IUCN codes are regarded as protected areas. See table below.

    SCENARIO = "ssp2_2050"
    ERA_YEAR = 2018                    # which year of the ERA Interim time series to use 
    RESCALE_ERA_TO_WIND_ATLAS = true   # Rescale the ERA Interim time series to fit annual wind speed averages from the Global Wind Atlas.
    SHOW_FIGURES = 0                   # 0, 1 or 2. Show no figures, a few, or all figures.

    #ONSHORECLASSES_min = 0:0.25:12.25
    #ONSHORECLASSES_max = 0.25:0.25:12.5
    ONSHORECLASSES_min = [2, 5, 6, 7, 8]
    ONSHORECLASSES_max = [5, 6, 7, 8, 99]
    ##ONSHORECLASSES_min = [4 5 6 7 8]
    ##ONSHORECLASSES_max = [5 6 7 8 99]

    #OFFSHORECLASSES_min = 0:0.25:12.25
    #OFFSHORECLASSES_max = 0.25:0.25:12.5
    OFFSHORECLASSES_min = [3, 6, 7, 8, 9]
    OFFSHORECLASSES_max = [6, 7, 8, 9, 99]
    ##OFFSHORECLASSES_min = [5 6 7 8 9]
    ##OFFSHORECLASSES_max = [6 7 8 9 99]

    # Land types
    #     0      'Water'                       
    #     1      'Evergreen Needleleaf Forests'
    #     2      'Evergreen Broadleaf Forests' 
    #     3      'Deciduous Needleleaf Forests'
    #     4      'Deciduous Broadleaf Forests' 
    #     5      'Mixed Forests'               
    #     6      'Closed Shrublands'           
    #     7      'Open Shrublands'             
    #     8      'Woody Savannas'              
    #     9      'Savannas'                    
    #    10      'Grasslands'                  
    #    11      'Permanent Wetlands'          
    #    12      'Croplands'                   
    #    13      'Urban'                       
    #    14      'Cropland/Natural'            
    #    15      'Snow/Ice'                    
    #    16      'Barren'                      

    # Protected areas (IUCN codes from the WDPA)
    #    1      'Not Applicable'    'Not Applicable'                 
    #    2      'Not Assigned'      'Not Assigned'                   
    #    3      'Not Reported'      'Not Reported'                   
    #    4      'Ia'                'Strict Nature Reserve'          
    #    5      'Ib'                'Wilderness Area'                
    #    6      'II'                'National Park'                  
    #    7      'III'               'Natural Monument'               
    #    8      'IV'                'Habitat/Species Management'     
    #    9      'V'                 'Protected Landscape/Seascape'   
    #    10     'VI'                'Managed Resource Protected Area'


    println("\nReading auxiliary datasets...")

    # bboxglobal = [-90 -180; 90 180]
    # if GISREGION == "MENA"
    #     bbox = [10 -15; 45 70]
    # elseif GISREGION[1:6] == "Europe"
    #     bbox = [34 -11; 72 32]
    # else
    #     bbox = bboxglobal
    # end
    # latrange, lonrange = bbox2ranges(bbox, 12)

    res = 0.01
    erares = 0.28125
    rasterdensity = round(Int, 1/res)   # pixels per degree

    regions, offshoreregions, regionlist = loadregions(GISREGION)

    # get indexes of the bounding box containing onshore region data with 1 degree of padding
    lonrange, latrange = getbboxranges(regions, rasterdensity)

    lons = (-180+res/2:res:180-res/2)[lonrange]         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = cosd.(lats) * (2*6371*pi/(360/res))^2    # area of a grid cell in km2   
    eralonranges, eralatrange = eraranges(lonrange, latrange)
    eralonindexes, eralatindexes = eraindexlookup(lons,lats,eralonranges, eralatrange)

    regions = regions[lonrange,latrange]
    offshoreregions = offshoreregions[lonrange,latrange]

    # path = joinpath(dirname(@__FILE__), "..")
    gridaccess = JLD.load("gridaccess_$(SCENARIO).jld", "gridaccess")[lonrange,latrange]
    pop = JLD.load("population_$(SCENARIO).jld", "population")[lonrange,latrange]
    topo = JLD.load("topography.jld", "topography")[lonrange,latrange]
    land = JLD.load("landcover.jld", "landcover")[lonrange,latrange]
    protected = JLD.load("protected.jld", "protected")[lonrange,latrange]
    windatlas = getwindatlas()[lonrange,latrange]

    yearlength = 8760 + 24*leapyear(ERA_YEAR)

    println("Reading ERA5 wind capacity factor timeseries...")

    @time meanwind, wind = h5open("D:/era5wind$ERA_YEAR.h5", "r") do file
        if length(eralonranges) == 1
            file["meanwind"][eralonranges[1], eralatrange],
            file["wind"][eralonranges[1], eralatrange, :]
        else
            [file["meanwind"][eralonranges[1], eralatrange]; file["meanwind"][eralonranges[2], eralatrange]],
            [file["wind"][eralonranges[1], eralatrange, :]; file["wind"][eralonranges[2], eralatrange, :]]
        end
    end




    println("Creating masks...")

    goodland = (regions .> 0)
    for i in EXCLUDE_LANDTYPES
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in PROTECTED_CODES
        protected_area[protected .== i] .= true
    end

    # Onshore wind A elec access within 9 km, hardcoded by data resolution
    gridA = (gridaccess .> 0.1)

    # Onshore wind B elec access > 9 km and < 9*ELEC_ACCESS_PIXELS km
    disk = diskfilterkernel(ELEC_ACCESS_PIXELS)
    gridB = (imfilter(gridaccess, disk) .> 0.1)

    # Onshore wind A elec access within 2*9 km
    # disk = diskfilterkernel(2)
    # gridA = (imfilter(gridaccess, disk) .> 0.1)

    # Onshore wind B elec access > 2*9 km and < 9*ELEC_ACCESS_PIXELS km
    # disk = diskfilterkernel(ELEC_ACCESS_PIXELS)
    # gridB = (imfilter(gridaccess, disk) .> 0.1)

    println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")
    # all mask conditions
    mask_onshoreA = gridA .& (pop .< PERSONS_PER_KM2) .& goodland .& .!protected_area
    mask_onshoreB = (gridB .& .!gridA) .& (pop .< PERSONS_PER_KM2) .& goodland .& .!protected_area

    # shoreline mask for offshore wind
    disk2 = diskfilterkernel(MIN_PIXELS_TO_SHORE)
    shore = (imfilter(regions .> 0, disk2) .> 0)

    # Offshore wind elec access same as onshore B (elec access < 9*ELEC_ACCESS_PIXELS km)
    # all mask conditions
    mask_offshore = gridB .& .!shore .& (topo .> -MAX_DEPTH) .& (offshoreregions .> 0) .& .!protected_area



    println("Allocating pixels to classes using the Global Wind Atlas...")

    nclasses = length(ONSHORECLASSES_min)
    onshoreclass = zeros(UInt8, size(windatlas))
    offshoreclass = zeros(UInt8, size(windatlas))
    for c = 1:nclasses
        onshoreclass[(windatlas .>= ONSHORECLASSES_min[c]) .& (windatlas .< ONSHORECLASSES_max[c])] .= c
        offshoreclass[(windatlas .>= OFFSHORECLASSES_min[c]) .& (windatlas .< OFFSHORECLASSES_max[c])] .= c
    end



    println("Calculating GW potential in wind classes...")

    nlons = size(regions,1)
    nlats = size(regions,2)
    numreg = length(regionlist)

    capacity_onshoreA = zeros(numreg,nclasses)
    capacity_onshoreB = zeros(numreg,nclasses)
    capacity_offshore = zeros(numreg,nclasses)

    for j = 1:nlats
        for i = 1:nlons
            reg = regions[i,j]
            area = cellarea[j]
            class = onshoreclass[i,j]
            offreg = offshoreregions[i,j]
            offclass = offshoreclass[i,j]
            if reg > 0 && class > 0 && mask_onshoreA[i,j] > 0
                capacity_onshoreA[reg,class] += 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * area
            elseif reg > 0 && class > 0 && mask_onshoreB[i,j] > 0
                capacity_onshoreB[reg,class] += 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * area
            elseif offreg > 0 && offclass > 0 && mask_offshore[i,j] > 0
                capacity_offshore[offreg,offclass] += 1/1000 * OFFSHORE_DENSITY * AREA_OFFSHORE * area
            end
        end
    end



    println("Calculating regional and global capacity factors of wind classes and representative timeseries...")
    println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")

    CFtime_windonshoreA, CFtime_windonshoreB, CFtime_windoffshore =
        # CF_windonshoreA_global, CF_windonshoreB_global, CF_windoffshore_global, CF_wind_time_agg
        calc_capacityfactors(wind, numreg, nclasses, regions, offshoreregions, onshoreclass, offshoreclass,
                eralonindexes, eralatindexes, mask_onshoreA, mask_onshoreB, mask_offshore)



    # firsttime1 = ncread(['D:/datasets/era5/era5wind' quartername(ERA_YEAR,1) '.nc'], 'time', 1, 1)
    # firsttime2 = ncread(['D:/datasets/era5/era5solar' quartername(ERA_YEAR,1) '.nc'], 'time', 1, 1)
    println("\n\nWind and solar ERA5 data not synced!!!!! Time diff 7 hours.")  #, [firsttime1 firsttime2])

    matopen("GISdata_wind$(ERA_YEAR)_$GISREGION.mat", "w") do file
        write(file, "CFtime_windonshoreA", CFtime_windonshoreA)
        write(file, "CFtime_windonshoreB", CFtime_windonshoreB)
        write(file, "CFtime_windoffshore", CFtime_windoffshore)
        write(file, "capacity_onshoreA", capacity_onshoreA)
        write(file, "capacity_onshoreB", capacity_onshoreB)
        write(file, "capacity_offshore", capacity_offshore)
    end

    return capacity_onshoreA, capacity_onshoreB, capacity_offshore
end

# 0 - 29 m/s
const windparkcurve = [
    0.0, 0.0014, 0.0071, 0.0229, 0.0545, 0.1067, 0.1831, 0.2850, 0.4085, 0.5434,
    0.6744, 0.7847, 0.8614, 0.9048, 0.9266, 0.9353, 0.9373, 0.9375, 0.9375, 0.9375,
    0.9375, 0.9375, 0.9375, 0.9311, 0.8683, 0.6416, 0.2948, 0.0688, 0.0063, 0.0
]

function speed2capacityfactor(windspeed)
    if windspeed >= 29 || windspeed < 0
        return 0.0
    end
    fw = floor(Int, windspeed)
    frac = windspeed - fw
    return (1-frac).*windparkcurve[fw+1] + frac.*windparkcurve[ceil(Int, windspeed)+1]
end

function calc_capacityfactors(wind, numreg, nclasses, regions, offshoreregions, onshoreclass, offshoreclass,
                eralonindexes, eralatindexes, mask_onshoreA, mask_onshoreB, mask_offshore)
    nlons, nlats = size(regions)
    yearlength = size(wind,3)

    CFtime_windonshoreA = zeros(yearlength,numreg,nclasses)
    CFtime_windonshoreB = zeros(yearlength,numreg,nclasses)
    CFtime_windoffshore = zeros(yearlength,numreg,nclasses)
    count_onshoreA = zeros(Int,numreg,nclasses)
    count_onshoreB = zeros(Int,numreg,nclasses)
    count_offshore = zeros(Int,numreg,nclasses)

    updateprogress = Progress(nlats, 1)
    for j = 1:nlats
        for i = 1:nlons
            reg = regions[i,j]
            class = onshoreclass[i,j]
            offreg = offshoreregions[i,j]
            offclass = offshoreclass[i,j]
            eralon, eralat = eralonindexes[i], eralatindexes[j]
            # can't use elseif here, probably some overlap in the masks
            if reg > 0 && class > 0 && mask_onshoreA[i,j] > 0
                CFtime_windonshoreA[:,reg,class] += speed2capacityfactor.(wind[eralon, eralat, :])
                count_onshoreA[reg,class] += 1
            elseif reg > 0 && class > 0 && mask_onshoreB[i,j] > 0
                CFtime_windonshoreB[:,reg,class] += speed2capacityfactor.(wind[eralon, eralat, :])
                count_onshoreB[reg,class] += 1
            elseif offreg > 0 && offclass > 0 && mask_offshore[i,j] > 0
                CFtime_windoffshore[:,offreg,offclass] += speed2capacityfactor.(wind[eralon, eralat, :])
                count_offshore[offreg,offclass] += 1
            end
        end
        next!(updateprogress)
    end

    # CFtime_windonshoreA_global = sumdrop(CFtime_windonshoreA, dims=2)
    # CFtime_windonshoreB_global = sumdrop(CFtime_windonshoreB, dims=2)
    # CFtime_windoffshore_global = sumdrop(CFtime_windoffshore, dims=2)
    # count_onshoreA_global = sumdrop(count_onshoreA, dims=1)
    # count_onshoreB_global = sumdrop(count_onshoreB, dims=1)
    # count_offshore_global = sumdrop(count_offshore, dims=1)
    # CF_wind_time_agg = sumdrop(CFtime_windonshoreA[:,:,3:5], dims=3) + CFtime_windonshoreA[:,:,5]
    # count_wind_agg = sumdrop(count_onshoreA[:,3:5], dims=2) + count_onshoreB[:,5]

    for y = 1:yearlength
        CFtime_windonshoreA[y,:,:] ./= count_onshoreA
        CFtime_windonshoreB[y,:,:] ./= count_onshoreB
        CFtime_windoffshore[y,:,:] ./= count_offshore
        # CFtime_windonshoreA_global[y,:] ./= count_onshoreA_global
        # CFtime_windonshoreB_global[y,:] ./= count_onshoreB_global
        # CFtime_windoffshore_global[y,:] ./= count_offshore_global
        # CF_wind_time_agg[y,:] ./= count_wind_agg
    end

    # CF_windonshoreA = meandrop(CFtime_windonshoreA, dims=1)
    # CF_windonshoreB = meandrop(CFtime_windonshoreB, dims=1)
    # CF_windoffshore = meandrop(CFtime_windoffshore, dims=1)
    # CF_windonshoreA_global = meandrop(CFtime_windonshoreA_global, dims=1)
    # CF_windonshoreB_global = meandrop(CFtime_windonshoreB_global, dims=1)
    # CF_windoffshore_global = meandrop(CFtime_windoffshore_global, dims=1)

    return CFtime_windonshoreA, CFtime_windonshoreB, CFtime_windoffshore 
            # CF_windonshoreA_global, CF_windonshoreB_global, CF_windoffshore_global, CF_wind_time_agg
end