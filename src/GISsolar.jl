function GISsolar()

    GISREGION = "Europe10"     # 'global', 'europe8/10/12/15/17', 'china6', 'eurochine14', 'eurasia38/21' or 'mena'

    use_FDIR = true            # Whether to use the ERA5 FDIR dataset for direct insolation (or the clearness model)

    PV_DENSITY = 45            # Solar PV land use 45 Wp/m2 = 45 MWp/km2
    CSP_DENSITY = 35           # CSP land use 35 W/m2
    PVROOF_AREA = .05          # area available for rooftop PV after the masks have been applied
    PLANT_AREA = .03           # area available for PV or CSP plants after the masks have been applied

    # (The datasets have a 5 minute resolution, i.e. each pixel is approximately 9.3 km at the equator.)
    # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
    # roughly half the people of the world live at density > 300 ppl/km2
    ELEC_ACCESS_PIXELS = 16         # Access to electricity < 16*9.3 km = 150 km  (for solar classes of category B)
    PVROOF_PERSONS_PER_KM2 = 100    # Only in populated areas, so AT LEAST x persons/km2
    PLANT_PERSONS_PER_KM2 = 75      # Not too crowded, max x persons/km2 (both PV and CSP plants)
    EXCLUDE_LANDTYPES = [0, 12]     # Not water or croplands.
    PROTECTED_CODES = [3, 4, 5, 6, 7, 8]    # These IUCN codes are regarded as protected areas. See table below.

    CSP_SOLAR_MULTIPLE = 2.5        # Ratio of collector peak power to generator power (see IRENA CSP cost page 8)
    CSP_STORAGE = 12                # 3, 6 or 9 hours of thermal storage.

    ERA_YEAR = 2016                    # which year of the ERA Interim time series to use 
    SHOW_FIGURES = 0                   # 0, 1 or 2. Show no figures, a few, or all figures.

    #PVCLASSES_min = 0:0.01:0.49
    #PVCLASSES_max = 0.01:0.01:0.50
    PVCLASSES_min = [0.05, 0.15, 0.20, 0.24, 0.28]
    PVCLASSES_max = [0.15, 0.20, 0.24, 0.28, 1.00]
    ##PVCLASSES_min = [0.10, 0.15, 0.20, 0.24, 0.28]
    ##PVCLASSES_max = [0.15, 0.20, 0.24, 0.28, 1.00]

    #CSPCLASSES.min = 0:0.01:0.49
    #CSPCLASSES.max = 0.01:0.01:0.50
    CSPCLASSES_min = [0.05, 0.15, 0.20, 0.24, 0.28]
    CSPCLASSES_max = [0.15, 0.20, 0.24, 0.28, 1.00]
    ##CSPCLASSES_min = [0.15, 0.20, 0.24, 0.28, 0.32]
    ##CSPCLASSES_max = [0.20, 0.24, 0.28, 0.32, 1.00]

    # IMPORTANT!! Our solar data is DNI, which means we currently should assume
    # CSP of type solar tower. If we want parabolic trough then we need to
    # recalculate the solar data (capacity factors will then become lower).

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

    # NOTE ON SOLAR UNIT: the solar irradiance sets are in kW/m2. Since the
    # irradiance value used to represent "standard testing conditions" for PV
    # is 1000 W/m2, the solar datasets also directly give the capacity factor.
    # Actual insolation can occasionally go above 1000 W/m2.

    # Ideally, we should make direct assumptions of PV module efficiency as a
    # function of air temperature (see e.g. Bett & Thornton appendix A2), but
    # for now efficiency is included in our assumption of PV_DENSITY. The
    # uncertainties in PV_DENSITY and PLANT_AREA are so large that efficiency
    # variations as a function of temperature don't matter.

    bboxglobal = [-90 -180; 90 180]
    if GISREGION == "MENA"
        bbox = [10 -15; 45 70]
    elseif GISREGION[1:6] == "Europe"
        bbox = [34 -11; 72 32]
    else
        bbox = bboxglobal
    end
    # bbox = bboxglobal
    bboxsmall = roundbbox(bbox, 32/9)
    latrange, lonrange = bbox2ranges(bbox, 12)
    latrangesmall, lonrangesmall = bbox2ranges(bboxsmall, 32/9)

    matlabvars = matread("C:/Stuff/NewGIS/testregions_$(GISREGION)_12.mat")
    regionlist = Symbol.(vec(matlabvars["regionlist"]))
    numreg = length(regionlist)
    regions = matlabvars["regions"][latrange,lonrange]
    regions[regions .== numreg+1] .= 0
    # R = matlabvars["R"]

    matlabvars = matread("C:/Stuff/NewGIS/testregions_$(GISREGION)_3.5556.mat")
    smallregions = matlabvars["regions"][latrangesmall,lonrangesmall]
    smallregions[smallregions .== numreg+1] .= 0

    matlabvars = matread("C:/Stuff/NewGIS/testgrid12.mat")
    gridaccess = matlabvars["gridaccess"][latrange,lonrange]

    matlabvars = matread("C:/Stuff/NewGIS/testpop12.mat")
    pop = matlabvars["pop"][latrange,lonrange]

    matlabvars = matread("C:/Stuff/NewGIS/testtopo12.mat")
    topo = matlabvars["topo"][latrange,lonrange]

    matlabvars = matread("C:/Stuff/NewGIS/testland12.mat")
    land = matlabvars["land"][latrange,lonrange]

    matlabvars = matread("C:/Stuff/NewGIS/testprotected12.mat")
    protected = matlabvars["protected"][latrange,lonrange]

    lats = getlats(bboxglobal, 12, true)[latrange]
    pixelarea = (2*6371*pi/(360*60/5))^2        # area in km2 of 5 min pixel at equator
    areamatkm = cosd.(lats) * pixelarea         # area of a 5 min grid cell in km2

    yearlength = 8760 + 24*leapyear(ERA_YEAR)

    println("\nReading ERA5 DNI capacity factor time series...\n")
    filename = "D:/datasets/era5/era5solarCF$(use_FDIR ? "_FDIR" : "")$ERA_YEAR.h5"
    if bbox == bboxglobal
        @time solarDNI = h5read_fast(filename, "solarDNI")
    else
        @time solarDNI = h5read(filename, "solarDNI", (1:yearlength, latrangesmall, lonrangesmall))
    end

    println("\nReading ERA5 GTI capacity factor time series...\n")
    if bbox == bboxglobal
        @time solarGTI = h5read_fast(filename, "solarGTI")
    else
        @time solarGTI = h5read(filename, "solarGTI", (1:yearlength, latrangesmall, lonrangesmall))
    end

    meanDNI = h5read(filename, "meanDNI", (latrangesmall, lonrangesmall))
    meanGTI = h5read(filename, "meanGTI", (latrangesmall, lonrangesmall))

    println("\nInterpolating ERA5 posting data to cell data...\n")

    meanDNI = shifthalfcell(meanDNI)
    meanGTI = shifthalfcell(meanGTI)
    shiftallcells!(solarDNI)
    shiftallcells!(solarGTI)




    # ***** Create masks ******

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
    mask_rooftop = gridA .& (pop .> PVROOF_PERSONS_PER_KM2) .& goodland .& .!protected_area
    mask_plantA = gridA .& (pop .< PLANT_PERSONS_PER_KM2) .& goodland .& .!protected_area
    mask_plantB = (gridB .& .!gridA) .& (pop .< PLANT_PERSONS_PER_KM2) .& goodland .& .!protected_area

    # resize the high resolution masks to fit the ERA data
    smallmask_rooftop = imresize(mask_rooftop, size(meanGTI))
    smallmask_plantA = imresize(mask_plantA, size(meanGTI))
    smallmask_plantB = imresize(mask_plantB, size(meanGTI))



    CF_csp_time = solarDNI
    CF_csp_annual = meandrop(CF_csp_time, dims=1)

    # upscale the low resolution solar data to fit the other 5 minute data
    # use the new CSP capacity factors with storage
    solarpv5 = max.(0, imresize(meanGTI, size(regions)))
    solarcsp5 = max.(0, imresize(CF_csp_annual, size(regions)))

    regsize = size(regions)
    pvrooftop = zeros(regsize)
    pvplantA = zeros(regsize)
    pvplantB = zeros(regsize)
    cspplantA = zeros(regsize)
    cspplantB = zeros(regsize)

    # make high resolution global maps for the main technologies
    pvrooftop[mask_rooftop.>0] = solarpv5[mask_rooftop.>0]
    pvplantA[mask_plantA.>0] = solarpv5[mask_plantA.>0]
    pvplantB[mask_plantB.>0] = solarpv5[mask_plantB.>0]
    cspplantA[mask_plantA.>0] = solarcsp5[mask_plantA.>0]
    cspplantB[mask_plantB.>0] = solarcsp5[mask_plantB.>0]



    # Allocate pixels to classes

    nclasses = length(PVCLASSES_min)
    pvclass = zeros(UInt8, regsize)
    cspclass = zeros(UInt8, regsize)
    pvclassERA = zeros(UInt8, size(smallregions))
    cspclassERA = zeros(UInt8, size(smallregions))
    for c = 1:nclasses
        pvclass[(solarpv5 .>= PVCLASSES_min[c]) .& (solarpv5 .< PVCLASSES_max[c])] .= c
        pvclassERA[(meanGTI .>= PVCLASSES_min[c]) .& (meanGTI .< PVCLASSES_max[c])] .= c
        cspclass[(solarcsp5 .>= CSPCLASSES_min[c]) .& (solarcsp5 .< CSPCLASSES_max[c])] .= c
        cspclassERA[(CF_csp_annual .>= CSPCLASSES_min[c]) .& (CF_csp_annual .< CSPCLASSES_max[c])] .= c
    end

    println("\nCalculating GW potential and capacity factors in PV and CSP classes...\n")

    nlats = size(regions,1)
    nlons = size(regions,2)
    nlatssmall = size(smallregions,1)
    nlonssmall = size(smallregions,2)

    capacity_pvrooftop = zeros(numreg,nclasses)
    capacity_pvplantA = zeros(numreg,nclasses)
    capacity_pvplantB = zeros(numreg,nclasses)
    capacity_cspplantA = zeros(numreg,nclasses)
    capacity_cspplantB = zeros(numreg,nclasses)
    CF_pvrooftop = zeros(numreg,nclasses)
    CF_pvplantA = zeros(numreg,nclasses)
    CF_pvplantB = zeros(numreg,nclasses)
    CF_cspplantA = zeros(numreg,nclasses)
    CF_cspplantB = zeros(numreg,nclasses) 
    CF_pvrooftop_global = zeros(nclasses)
    CF_pvplantA_global = zeros(nclasses)
    CF_pvplantB_global = zeros(nclasses)
    CF_cspplantA_global = zeros(nclasses)
    CF_cspplantB_global = zeros(nclasses) 
    count_pvrooftop = zeros(Int,numreg,nclasses)
    count_pvplantA = zeros(Int,numreg,nclasses)
    count_pvplantB = zeros(Int,numreg,nclasses)
    count_cspplantA = zeros(Int,numreg,nclasses)
    count_cspplantB = zeros(Int,numreg,nclasses)
    count_pvrooftop_global = zeros(Int,nclasses)
    count_pvplantA_global = zeros(Int,nclasses)
    count_pvplantB_global = zeros(Int,nclasses)
    count_cspplantA_global = zeros(Int,nclasses)
    count_cspplantB_global = zeros(Int,nclasses) 

    # Working at high resolution here ...
    for i = 1:nlats
        for j = 1:nlons
            reg = regions[i,j]
            area = areamatkm[i]
            class = pvclass[i,j]
            if reg > 0 && class > 0 && mask_rooftop[i,j] > 0
                capacity_pvrooftop[reg,class] += area
                CF_pvrooftop[reg,class] += solarpv5[i,j]
                count_pvrooftop[reg,class] += 1
                CF_pvrooftop_global[class] += solarpv5[i,j]
                count_pvrooftop_global[class] += 1
            elseif reg > 0 && class > 0 && mask_plantA[i,j] > 0
                capacity_pvplantA[reg,class] += area
                CF_pvplantA[reg,class] += solarpv5[i,j]
                count_pvplantA[reg,class] += 1
                CF_pvplantA[class] += solarpv5[i,j]
                count_pvplantA_global[class] += 1
            elseif reg > 0 && class > 0 && mask_plantB[i,j] > 0
                capacity_pvplantB[reg,class] += area
                CF_pvplantB[reg,class] += solarpv5[i,j]
                count_pvplantB[reg,class] += 1
                CF_pvplantB[class] += solarpv5[i,j]
                count_pvplantB_global[class] += 1
            end

            class = cspclass[i,j]
            if reg > 0 && class > 0 && mask_plantA[i,j] > 0
                capacity_cspplantA[reg,class] += area
                CF_cspplantA[reg,class] += solarcsp5[i,j]
                count_cspplantA[reg,class] += 1
                CF_cspplantA[class] += solarcsp5[i,j]
                count_cspplantA_global[class] += 1
            elseif reg > 0 && class > 0 && mask_plantB[i,j] > 0
                capacity_cspplantB[reg,class] += area
                CF_cspplantB[reg,class] += solarcsp5[i,j]
                count_cspplantB[reg,class] += 1
                CF_cspplantB[class] += solarcsp5[i,j]
                count_cspplantB_global[class] += 1
            end
        end
    end
    capacity_pvrooftop = 1/1000 * PV_DENSITY * PVROOF_AREA * capacity_pvrooftop
    capacity_pvplantA = 1/1000 * PV_DENSITY * PLANT_AREA * capacity_pvplantA
    capacity_pvplantB = 1/1000 * PV_DENSITY * PLANT_AREA * capacity_pvplantB
    capacity_cspplantA = 1/1000 * CSP_DENSITY * PLANT_AREA * capacity_cspplantA
    capacity_cspplantB = 1/1000 * CSP_DENSITY * PLANT_AREA * capacity_cspplantB
    CF_pvrooftop = CF_pvrooftop ./ count_pvrooftop
    CF_pvplantA = CF_pvplantA ./ count_pvplantA
    CF_pvplantB = CF_pvplantB ./ count_pvplantB
    CF_cspplantA = CF_cspplantA ./ count_cspplantA
    CF_cspplantB = CF_cspplantB ./ count_cspplantB
    CF_pvrooftop_global = CF_pvrooftop_global ./ count_pvrooftop_global
    CF_pvplantA_global = CF_pvplantA_global ./ count_pvplantA_global
    CF_pvplantB_global = CF_pvplantB_global ./ count_pvplantB_global
    CF_cspplantA_global = CF_cspplantA_global ./ count_cspplantA_global
    CF_cspplantB_global = CF_cspplantB_global ./ count_cspplantB_global



    println("\nCalculating capacity factors of solar PV & CSP for each region & class...\n")

    CFtime_pvrooftop = zeros(yearlength,numreg,nclasses)
    CFtime_pvplantA = zeros(yearlength,numreg,nclasses)
    CFtime_pvplantB = zeros(yearlength,numreg,nclasses)
    CFtime_cspplantA = zeros(yearlength,numreg,nclasses)
    CFtime_cspplantB = zeros(yearlength,numreg,nclasses)
    count_pvrooftop = zeros(Int,numreg,nclasses)
    count_pvplantA = zeros(Int,numreg,nclasses)
    count_pvplantB = zeros(Int,numreg,nclasses)
    count_cspplantA = zeros(Int,numreg,nclasses)
    count_cspplantB = zeros(Int,numreg,nclasses)

    # Working at low resolution here ...
    updateprogress = Progress(nlatssmall, 1)
    for i = 1:nlatssmall
        for j = 1:nlonssmall
            reg = smallregions[i,j]

            class = pvclassERA[i,j]      
            # can't use elseif here, probably some overlap in the masks
            if reg > 0 && class > 0 && smallmask_rooftop[i,j] > 0
                CFtime_pvrooftop[:,reg,class] += solarGTI[:,i,j]
                count_pvrooftop[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantA[i,j] > 0
                CFtime_pvplantA[:,reg,class] += solarGTI[:,i,j]
                count_pvplantA[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantB[i,j] > 0
                CFtime_pvplantB[:,reg,class] += solarGTI[:,i,j]
                count_pvplantB[reg,class] += 1
            end

            class = cspclassERA[i,j]
            if reg > 0 && class > 0 && smallmask_plantA[i,j] > 0
                CFtime_cspplantA[:,reg,class] += CF_csp_time[:,i,j]
                count_cspplantA[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantB[i,j] > 0
                CFtime_cspplantB[:,reg,class] += CF_csp_time[:,i,j]
                count_cspplantB[reg,class] += 1
            end            
        end
        next!(updateprogress)
    end
    for y = 1:yearlength
        CFtime_pvrooftop[y,:,:] ./= count_pvrooftop
        CFtime_pvplantA[y,:,:] ./= count_pvplantA
        CFtime_pvplantB[y,:,:] ./= count_pvplantB
        CFtime_cspplantA[y,:,:] ./= count_cspplantA
        CFtime_cspplantB[y,:,:] ./= count_cspplantB
    end



    println("\nCalculating representative timeseries for solar...\n")
    CF_pv_time_agg = zeros(yearlength,numreg)
    CF_csp_time_agg = zeros(yearlength,numreg)
    count_pv_agg = zeros(numreg)
    count_csp_agg = zeros(numreg)

    # Calculate representative timeseries for each technology and region.
    # In every region, choose all sites in onshore classes A3-A5 # B5 and offshore class 5,
    # and average power output(t) over the entire region.
    updateprogress = Progress(nlatssmall, 1)
    # Working at low resolution here ...
    for i = 1:nlatssmall
        for j = 1:nlonssmall
            reg = smallregions[i,j]

            class = pvclassERA[i,j]
            if (reg > 0 && smallmask_rooftop[i,j] > 0 && class >= 3) ||
                (reg > 0 && smallmask_plantA[i,j] > 0 && class >= 3) ||
                (reg > 0 && smallmask_plantB[i,j] > 0 && class >= 5)
                    CF_pv_time_agg[:,reg] += solarGTI[:,i,j]
                    count_pv_agg[reg] += 1
            end

            class = cspclassERA[i,j]
            if (reg > 0 && smallmask_plantA[i,j] > 0 && class >= 3) ||
                (reg > 0 && smallmask_plantB[i,j] > 0 && class >= 5)
                    CF_csp_time_agg[:,reg] += CF_csp_time[:,i,j]
                    count_csp_agg[reg] += 1
            end
        end
        next!(updateprogress)
    end
    for r = 1:numreg
        CF_pv_time_agg[:,r] ./= count_pv_agg[r]
        CF_csp_time_agg[:,r] ./= count_csp_agg[r]
    end

    matopen("GISdata_solar$(ERA_YEAR)_$GISREGION.mat", "w") do file
        write(file, "CFtime_pvrooftop", CFtime_pvrooftop)
        write(file, "CFtime_pvplantA", CFtime_pvplantA)
        write(file, "CFtime_pvplantB", CFtime_pvplantB)
        write(file, "CFtime_cspplantA", CFtime_cspplantA)
        write(file, "CFtime_cspplantB", CFtime_cspplantB)
        write(file, "capacity_pvrooftop", capacity_pvrooftop)
        write(file, "capacity_pvplantA", capacity_pvplantA)
        write(file, "capacity_pvplantB", capacity_pvplantB)
        write(file, "capacity_cspplantA", capacity_cspplantA)
        write(file, "capacity_cspplantB", capacity_cspplantB)
    end

    return capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end
