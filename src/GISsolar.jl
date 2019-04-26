solaroptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"

    :pv_density => 45,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2
    :csp_density => 35,                 # CSP land use 35 W/m2

    :pvroof_area => .05,                # area available for rooftop PV after the masks have been applied
    :plant_area => .03,                 # area available for PV or CSP plants after the masks have been applied

    :distance_elec_access => 150,       # max distance to grid [km] (for solar classes of category B)
    :plant_persons_per_km2 => 75,       # not too crowded, max X persons/km2 (both PV and CSP plants)
    :pvroof_persons_per_km2 => 100,     # only in populated areas, so AT LEAST x persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :exclude_landtypes => [0,12],       # exclude water and croplands. See codes in table below.
    :protected_codes => [3,4,5,6,7,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :csp_solar_multiple => 2.5,         # ratio of collector peak power to generator power (see IRENA CSP cost page 8)
    :csp_storage_hours => 12,           # 3, 6 or 9 hours of thermal storage.

    :scenario => "ssp2_2050",           # default scenario for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.05,0.15,0.20,0.24,0.28],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49]
    :pvclasses_max => [0.15,0.20,0.24,0.28,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50]
    :cspclasses_min => [0.05,0.15,0.20,0.24,0.28],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.15,0.20,0.24,0.28,1.00]  # upper bound on annual CSP capacity factor for class X
)
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

mutable struct SolarOptions
    gisregion               ::String
    pv_density              ::Float64           # W/m2
    csp_density             ::Float64           # W/m2
    pvroof_area             ::Float64           # share [0-1]
    plant_area              ::Float64           # share [0-1]
    distance_elec_access    ::Float64           # km
    plant_persons_per_km2   ::Float64           # persons/km2
    pvroof_persons_per_km2  ::Float64           # persons/km2
    exclude_landtypes       ::Vector{Int}
    protected_codes         ::Vector{Int}
    csp_solar_multiple      ::Float64
    csp_storage_hours       ::Float64           # h
    scenario                ::String
    era_year                ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    pvclasses_min           ::Vector{Float64}
    pvclasses_max           ::Vector{Float64}
    cspclasses_min          ::Vector{Float64}
    cspclasses_max          ::Vector{Float64}
end

SolarOptions() = SolarOptions("",0,0,0,0,0,0,0,[],[],0,0,"",0,0,0,[],[],[],[])

function SolarOptions(d::Dict{Symbol,Any})
    options = SolarOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function GISsolar(; optionlist...)

    # IMPORTANT!! The function makesolarera5() uses ERA5 solar datasets to
    # calculate Global Tilted Irradiance (GTI) for solar PV and Direct Normal
    # Irradiance (DNI) for CSP solar towers. If we want CSP parabolic troughs
    # then we need to add that dataset in makesolarera5() (and capacity factors
    # will become somewhat lower).

    # NOTE ON SOLAR UNIT: the solar irradiance sets are in kW/m2. Since the
    # irradiance value used to represent "standard testing conditions" for PV
    # is 1000 W/m2, the solar datasets also directly give the capacity factor.
    # Actual insolation can occasionally go above 1000 W/m2.

    # Ideally, we should make direct assumptions of PV module efficiency as a
    # function of air temperature (see e.g. Bett & Thornton appendix A2), but
    # for now efficiency is included in our assumption of :pv_density. Wind
    # speed also affects PV module temperature and efficiency. However, the
    # uncertainties in :pv_density and :plant_area are so large that efficiency
    # variations as a function of temperature don't matter.

    options = SolarOptions(merge(solaroptions(), optionlist))

    regions, offshoreregions, regionlist, gridaccess, pop, topo, land, protected, lonrange, latrange = read_datasets(options)
    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(options, lonrange, latrange)

return meanGTI, solarGTI, meanDNI, solarDNI
nclasses, pvclass, cspclass = makesolarclasses(options, meanGTI, meanDNI)

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(options, regions, gridaccess, pop, land, protected)

    windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore =
        calc_solar_vars(options, solarGHI, solarDHI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    matopen("GISdata_solar$(ERA_YEAR)_$GISREGION.mat", "w") do file
        write(file, "CFtime_pvrooftop", CF_pvrooftop)
        write(file, "CFtime_pvplantA", CF_pvplantA)
        write(file, "CFtime_pvplantB", CF_pvplantB)
        write(file, "CFtime_cspplantA", CF_cspplantA)
        write(file, "CFtime_cspplantB", CF_cspplantB)
        write(file, "capacity_pvrooftop", capacity_pvrooftop)
        write(file, "capacity_pvplantA", capacity_pvplantA)
        write(file, "capacity_pvplantB", capacity_pvplantB)
        write(file, "capacity_cspplantA", capacity_cspplantA)
        write(file, "capacity_cspplantB", capacity_cspplantB)
    end

    return capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function read_solar_datasets(options, lonrange, latrange)
    @unpack res, erares, era_year = options

    println("Reading ERA5 solar datasets...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    @time meanGTI, solarGTI, meanDNI, solarDNI = h5open("D:/era5solar$era_year.h5", "r") do file
        if length(eralonranges) == 1
            file["meanGTI"][eralonranges[1], eralatrange],
                file["GTI"][:,eralonranges[1], eralatrange],
                file["meanDNI"][eralonranges[1], eralatrange],
                file["DNI"][:,eralonranges[1], eralatrange]
        else
            [file["meanGTI"][eralonranges[1], eralatrange]; file["meanGTI"][eralonranges[2], eralatrange]],
                [file["GTI"][:, eralonranges[1], eralatrange]; file["GTI"][:, eralonranges[2], eralatrange]],
                [file["meanDNI"][eralonranges[1], eralatrange]; file["meanDNI"][eralonranges[2], eralatrange]],
                [file["DNI"][:, eralonranges[1], eralatrange]; file["DNI"][:, eralonranges[2], eralatrange]]
        end
    end
    return meanGTI, solarGTI, meanDNI, solarDNI
end

function create_solar_masks(options, regions, gridaccess, pop, land, protected)
    @unpack res, exclude_landtypes, protected_codes, distance_elec_access, plant_persons_per_km2, pvroof_persons_per_km2 = options

    println("Creating masks...")

    goodland = (regions .> 0)
    for i in exclude_landtypes
        goodland[land .== i] .= false
    end
    protected_area = zeros(Bool, size(protected))
    for i in protected_codes
        protected_area[protected .== i] .= true
    end

    # Pixels with electricity access for onshore wind A 
    gridA = (gridaccess .> 0.1)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    gridB = (imfilter(gridaccess, disk) .> 0.1)

    println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    mask_rooftop = gridA .& (pop .> pvroof_persons_per_km2) .& goodland .& .!protected_area
    mask_plantA = gridA .& (pop .< plant_persons_per_km2) .& goodland .& .!protected_area
    mask_plantB = (gridB .& .!gridA) .& (pop .< plant_persons_per_km2) .& goodland .& .!protected_area

    return mask_rooftop, mask_plantA, mask_plantB
end


function calc_everything_5min(options, windatlas, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    # INIT

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

    CF_pvrooftop = zeros(yearlength,numreg,nclasses)
    CF_pvplantA = zeros(yearlength,numreg,nclasses)
    CF_pvplantB = zeros(yearlength,numreg,nclasses)
    CF_cspplantA = zeros(yearlength,numreg,nclasses)
    CF_cspplantB = zeros(yearlength,numreg,nclasses)
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
                CF_pvrooftop[:,reg,class] += solarGTI[:,i,j]
                count_pvrooftop[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantA[i,j] > 0
                CF_pvplantA[:,reg,class] += solarGTI[:,i,j]
                count_pvplantA[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantB[i,j] > 0
                CF_pvplantB[:,reg,class] += solarGTI[:,i,j]
                count_pvplantB[reg,class] += 1
            end

            class = cspclassERA[i,j]
            if reg > 0 && class > 0 && smallmask_plantA[i,j] > 0
                CF_cspplantA[:,reg,class] += CF_csp_time[:,i,j]
                count_cspplantA[reg,class] += 1
            end
            if reg > 0 && class > 0 && smallmask_plantB[i,j] > 0
                CF_cspplantB[:,reg,class] += CF_csp_time[:,i,j]
                count_cspplantB[reg,class] += 1
            end            
        end
        next!(updateprogress)
    end
    for y = 1:yearlength
        CF_pvrooftop[y,:,:] ./= count_pvrooftop
        CF_pvplantA[y,:,:] ./= count_pvplantA
        CF_pvplantB[y,:,:] ./= count_pvplantB
        CF_cspplantA[y,:,:] ./= count_cspplantA
        CF_cspplantB[y,:,:] ./= count_cspplantB
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
end

function increment_solarCF!(cf::AbstractVector{<:AbstractFloat}, speed_or_cf::AbstractVector{<:AbstractFloat}, factor, rescale::Bool)
    if rescale
        @inbounds for i = 1:length(cf)
            cf[i] += speed2capacityfactor(speed_or_cf[i]*factor)
        end
    else
        @inbounds for i = 1:length(cf)
            cf[i] += speed_or_cf[i]
        end
    end
end

function makesolarclasses(options, meanGTI, meanDNI)
    println("Allocating pixels to classes using ERA5 annual means...")

    @unpack pvclasses_min, pvclasses_max, cspclasses_min, cspclasses_max = options

    nclasses = length(onshoreclasses_min)
    pvclass = zeros(UInt8, size(meanGTI))
    cspclass = zeros(UInt8, size(meanDNI))
    for c = 1:nclasses
        pvclass[(meanGTI .>= pvclasses_min[c]) .& (meanGTI .< pvclasses_max[c])] .= c
        cspclass[(meanDNI .>= cspclasses_min[c]) .& (meanDNI .< cspclasses_max[c])] .= c
    end

    return nclasses, pvclass, cspclass
end

function calc_solar_vars(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and class...")
    println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    nclasses, pvclass, cspclass = makesolarclasses(options, meanGTI, meanDNI)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(options, lonrange, latrange)

    @unpack era_year, res, erares, pv_density, csp_density, pvroof_area, plant_area = options

    numreg = length(regionlist)
    yearlength, nlons, nlats = size(windspeed)
    firsttime = DateTime(era_year, 1, 1)

    capacity_pvrooftop = zeros(numreg,nclasses)
    capacity_pvplantA = zeros(numreg,nclasses)
    capacity_pvplantB = zeros(numreg,nclasses)
    capacity_cspplantA = zeros(numreg,nclasses)
    capacity_cspplantB = zeros(numreg,nclasses)
    CF_pvrooftop = zeros(yearlength,numreg,nclasses)
    CF_pvplantA = zeros(yearlength,numreg,nclasses)
    CF_pvplantB = zeros(yearlength,numreg,nclasses)
    CF_cspplantA = zeros(yearlength,numreg,nclasses)
    CF_cspplantB = zeros(yearlength,numreg,nclasses)
    count_pvrooftop = zeros(Int,numreg,nclasses)
    count_pvplantA = zeros(Int,numreg,nclasses)
    count_pvplantB = zeros(Int,numreg,nclasses)
    count_cspplantA = zeros(Int,numreg,nclasses)
    count_cspplantB = zeros(Int,numreg,nclasses)

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for i = 1:nlons  #in randperm(nlons)
        eralon = eralons[i]
        rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2-res/5, res)]
        for j = 1:nlats
            meanwind[i,j] == 0 && continue
            wind = rescale_to_wind_atlas ? windspeed[:, i, j] : speed2capacityfactor.(windspeed[:, i, j])
            eralat = eralats[j]
            # get all high resolution row and column indexes within this ERA5 cell 
            colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2+res/5, res)]   
            for t = 1:yearlength
                datetime = firsttime + Hour(t)
                zenith, azimuth = solarposition_grena1(datetime, eralat, eralon)
                for c in colrange, r in rowrange
                    (c == 0 || r == 0) && continue
                    reg = regions[r,c]
                    area = cellarea[c]
                    class = onshoreclass[r,c]
                    offreg = offshoreregions[r,c]
                    offclass = offshoreclass[r,c]
                    
                    # can't use elseif here, probably some overlap in the masks
                    if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                        capacity_onshoreA[reg,class] += 1/1000 * onshore_density * area_onshore * area
                        incrementCF!(windCF_onshoreA[:,reg,class], wind, windatlas[r,c] / meanwind[i,j], rescale_to_wind_atlas)
                        count_onshoreA[reg,class] += 1
                    elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                        capacity_onshoreB[reg,class] += 1/1000 * onshore_density * area_onshore * area
                        incrementCF!(windCF_onshoreB[:,reg,class], wind, windatlas[r,c] / meanwind[i,j], rescale_to_wind_atlas)
                        count_onshoreB[reg,class] += 1
                    elseif offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                        capacity_offshore[offreg,offclass] += 1/1000 * offshore_density * area_offshore * area
                        incrementCF!(windCF_offshore[:,offreg,offclass], wind, windatlas[r,c] / meanwind[i,j], rescale_to_wind_atlas)
                        count_offshore[offreg,offclass] += 1
                    end
                end
            end
        end
        next!(updateprogress)
    end

    for y = 1:yearlength
        windCF_onshoreA[y,:,:] ./= count_onshoreA
        windCF_onshoreB[y,:,:] ./= count_onshoreB
        windCF_offshore[y,:,:] ./= count_offshore
    end

    return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

