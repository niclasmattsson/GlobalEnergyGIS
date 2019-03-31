function GISwind()

    GISREGION = "Eurasia38"     # 'global', 'europe8/10/12/15/17', 'china6', 'eurochine14', 'eurasia38/21' or 'mena'

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
    ERA_YEAR = 2016                    # which year of the ERA Interim time series to use 
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
    rasterdensity = round(Int, 1/res)   # pixels per degree

    regions, offshoreregions, regionlist = loadregions(GISREGION)
    lonrange, latrange = getbbox(regions, rasterdensity)  # get indexes of the bounding box containing region data with 1 degree of padding
    regions = regions[lonrange,latrange]
    offshoreregions = offshoreregions[lonrange,latrange]

    # path = joinpath(dirname(@__FILE__), "..")
    gridaccess = JLD.load("gridaccess_$(SCENARIO).jld", "gridaccess")[lonrange,latrange]
    pop = JLD.load("population_$(SCENARIO).jld", "population")[lonrange,latrange]
    topo = JLD.load("topography.jld", "topography")[lonrange,latrange]
    land = JLD.load("landcover.jld", "landcover")[lonrange,latrange]
    protected = JLD.load("protected.jld", "protected")[lonrange,latrange]
    windatlas = getwindatlas()[lonrange,latrange]

    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = cosd.(lats) * (2*6371*pi/(360/res))^2    # area of a grid cell in km2

    yearlength = 8760 + 24*leapyear(ERA_YEAR)

    println("\nReading ERA5 wind capacity factor timeseries...\n")

    CFtype = RESCALE_ERA_TO_WIND_ATLAS ? "RCF" : "CF"
    filename = "D:/datasets/era5/era5wind$CFtype$ERA_YEAR.h5"
    dataset = "wind$CFtype"
    if bbox == bboxglobal
        @time windCF = h5read_fast(filename, dataset)
    else
        @time windCF = h5read(filename, dataset, (1:yearlength, latrangesmall, lonrangesmall))
    end
    meanwind = h5read("D:/datasets/era5/era5wind$ERA_YEAR.h5", "meanwind", (latrangesmall, lonrangesmall))

    println("\nInterpolating ERA5 posting data to cell data...\n")

    meanwind = shifthalfcell(meanwind)
    shiftallcells!(windCF)




    # ***** Create masks for onshore wind ******

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

    # all mask conditions
    mask_onshoreA = gridA .& (pop .< PERSONS_PER_KM2) .& goodland .& .!protected_area
    mask_onshoreB = (gridB .& .!gridA) .& (pop .< PERSONS_PER_KM2) .& goodland .& .!protected_area

    # ***** Create masks for offshore wind ******

    disk2 = diskfilterkernel(MIN_PIXELS_TO_SHORE)
    shore = (imfilter(regions .> 0, disk2) .> 0)

    # Offshore wind elec access same as onshore B (elec access < 9*ELEC_ACCESS_PIXELS km)
    # all mask conditions
    mask_offshore = gridB .& .!shore .& (topo .> -MAX_DEPTH) .& (offshoreregions .> 0) .& .!protected_area

    # resize the high resolution masks to fit the ERA data
    smallmask_onshoreA = imresize(mask_onshoreA, size(meanwind))
    smallmask_onshoreB = imresize(mask_onshoreB, size(meanwind))
    smallmask_offshore = imresize(mask_offshore, size(meanwind))



    # Allocate pixels to classes using the Global Wind Atlas

    nclasses = length(ONSHORECLASSES_min)
    onshoreclass = zeros(UInt8, size(windatlas))
    offshoreclass = zeros(UInt8, size(windatlas))
    onshoreclassERA = zeros(UInt8, size(windatlassmall))
    offshoreclassERA = zeros(UInt8, size(windatlassmall))
    for c = 1:nclasses
        onshoreclass[(windatlas .>= ONSHORECLASSES_min[c]) .& (windatlas .< ONSHORECLASSES_max[c])] .= c
        onshoreclassERA[(windatlassmall .>= ONSHORECLASSES_min[c]) .& (windatlassmall .< ONSHORECLASSES_max[c])] .= c
        offshoreclass[(windatlas .>= OFFSHORECLASSES_min[c]) .& (windatlas .< OFFSHORECLASSES_max[c])] .= c
        offshoreclassERA[(windatlassmall .>= OFFSHORECLASSES_min[c]) .& (windatlassmall .< OFFSHORECLASSES_max[c])] .= c
    end



    println("\nCalculating GW potential in wind classes...\n")

    nlats = size(regions,1)
    nlons = size(regions,2)
    nlatssmall = size(smallregions,1)
    nlonssmall = size(smallregions,2)

    capacity_onshoreA = zeros(numreg,nclasses)
    capacity_onshoreB = zeros(numreg,nclasses)
    capacity_offshore = zeros(numreg,nclasses)

    # Working at high resolution here ...
    for i = 1:nlats
        for j = 1:nlons
            reg = regions[i,j]
            area = cellarea[i]
            class = onshoreclass[i,j]
            offreg = offshoreregions[i,j]
            offclass = offshoreclass[i,j]
            if reg > 0 && class > 0 && mask_onshoreA[i,j] > 0
                capacity_onshoreA[reg,class] += area
            elseif reg > 0 && class > 0 && mask_onshoreB[i,j] > 0
                capacity_onshoreB[reg,class] += area
            elseif offreg > 0 && offclass > 0 && mask_offshore[i,j] > 0
                capacity_offshore[offreg,offclass] += area
            end
        end
    end
    capacity_onshoreA = 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * capacity_onshoreA
    capacity_onshoreB = 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * capacity_onshoreB
    capacity_offshore = 1/1000 * OFFSHORE_DENSITY * AREA_OFFSHORE * capacity_offshore




    println("\nCalculating regional and global capacity factors of wind classes...\n")

    CFtime_windonshoreA = zeros(yearlength,numreg,nclasses)
    CFtime_windonshoreB = zeros(yearlength,numreg,nclasses)
    CFtime_windoffshore = zeros(yearlength,numreg,nclasses)
    CFtime_windonshoreA_global = zeros(yearlength,nclasses)
    CFtime_windonshoreB_global = zeros(yearlength,nclasses)
    CFtime_windoffshore_global = zeros(yearlength,nclasses)
    count_onshoreA = zeros(Int,numreg,nclasses)
    count_onshoreB = zeros(Int,numreg,nclasses)
    count_offshore = zeros(Int,numreg,nclasses)
    count_onshoreA_global = zeros(Int,nclasses)
    count_onshoreB_global = zeros(Int,nclasses)
    count_offshore_global = zeros(Int,nclasses)

    # Working at low resolution here ...
    updateprogress = Progress(nlatssmall, 1)
    for i = 1:nlatssmall
        for j = 1:nlonssmall
            reg = smallregions[i,j]
            class = onshoreclassERA[i,j]
            offreg = smalloffshoreregions[i,j]
            offclass = offshoreclassERA[i,j]
            # can't use elseif here, probably some overlap in the masks
            if reg > 0 && class > 0 && smallmask_onshoreA[i,j] > 0
                CFtime_windonshoreA_global[:,class] += windCF[:,i,j]
                count_onshoreA_global[class] += 1
                CFtime_windonshoreA[:,reg,class] += windCF[:,i,j]
                count_onshoreA[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_onshoreB[i,j] > 0
                CFtime_windonshoreB_global[:,class] += windCF[:,i,j]
                count_onshoreB_global[class] += 1
                CFtime_windonshoreB[:,reg,class] += windCF[:,i,j]
                count_onshoreB[reg,class] += 1
            end
            if offreg > 0 && offclass > 0 && smallmask_offshore[i,j] > 0
                CFtime_windoffshore_global[:,offclass] += windCF[:,i,j]
                count_offshore_global[offclass] += 1
                CFtime_windoffshore[:,offreg,offclass] += windCF[:,i,j]
                count_offshore[offreg,offclass] += 1
            end
        end
        next!(updateprogress)
    end
    for y = 1:yearlength
        CFtime_windonshoreA[y,:,:] ./= count_onshoreA
        CFtime_windonshoreB[y,:,:] ./= count_onshoreB
        CFtime_windoffshore[y,:,:] ./= count_offshore
        CFtime_windonshoreA_global[y,:] ./= count_onshoreA_global
        CFtime_windonshoreB_global[y,:] ./= count_onshoreB_global
        CFtime_windoffshore_global[y,:] ./= count_offshore_global
    end

    CF_windonshoreA = meandrop(CFtime_windonshoreA, dims=1)
    CF_windonshoreB = meandrop(CFtime_windonshoreB, dims=1)
    CF_windoffshore = meandrop(CFtime_windoffshore, dims=1)
    CF_windonshoreA_global = meandrop(CFtime_windonshoreA_global, dims=1)
    CF_windonshoreB_global = meandrop(CFtime_windonshoreB_global, dims=1)
    CF_windoffshore_global = meandrop(CFtime_windoffshore_global, dims=1)




    println("\nCalculating representative timeseries for wind...\n")
    CF_wind_time_agg = zeros(yearlength,numreg)
    count_wind_agg = zeros(numreg)

    # Calculate representative timeseries for each technology and region.
    # In every region, choose all sites in onshore classes A3-A5 # B5 and offshore class 5,
    # and average power output(t) over the entire region.
    updateprogress = Progress(nlatssmall, 1)
    # Working at low resolution here ...
    for i = 1:nlatssmall
        for j = 1:nlonssmall
            reg = smallregions[i,j]
            class = onshoreclassERA[i,j]
            offreg = smalloffshoreregions[i,j]
            offclass = offshoreclassERA[i,j]

            if (reg > 0 && smallmask_onshoreA[i,j] > 0 && class >= 3) ||
                (reg > 0 && smallmask_onshoreB[i,j] > 0 && class >= 5)
                    CF_wind_time_agg[:,reg] += windCF[:,i,j]
                    count_wind_agg[reg] += 1
            # add offshore class 4 for PAS since everything else is empty
            # elseif smallmask_offshore[i,j] > 0 && (offreg > 0 && offclass >= 5 || offreg == 9 && offclass >= 4)
            #     CF_wind_time_agg[:,offreg] += windCF[:,i,j]
            #     count_wind_agg[offreg] += 1
            end        
        end
        next!(updateprogress)
    end
    for r = 1:numreg
        CF_wind_time_agg[:,r] ./= count_wind_agg[r]
    end


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
