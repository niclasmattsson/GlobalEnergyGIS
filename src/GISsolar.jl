solaroptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat" 

    :pv_density => 45,                  # Solar PV land use 45 Wp/m2 = 45 MWp/km2 (includes PV efficiency & module spacing, add latitude dependency later)
    :csp_density => 35,                 # CSP land use 35 W/m2

    :pvroof_area => .05,                # area available for rooftop PV after the masks have been applied
    :plant_area => .05,                 # area available for PV or CSP plants after the masks have been applied

    :distance_elec_access => 300,       # max distance to grid [km] (for solar classes of category B)
    :plant_persons_per_km2 => 150,      # not too crowded, max X persons/km2 (both PV and CSP plants)
    :pvroof_persons_per_km2 => 200,     # only in populated areas, so AT LEAST x persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :exclude_landtypes => [0,1,2,3,4,5,8],       # exclude water and forests. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.08,0.14,0.18,0.22,0.26],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49;]
    :pvclasses_max => [0.14,0.18,0.22,0.26,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50;]
    :cspclasses_min => [0.10,0.18,0.24,0.28,0.32],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.18,0.24,0.28,0.32,1.00],  # upper bound on annual CSP capacity factor for class X

    :downsample_masks => 1              # set to 2 or higher to scale down mask sizes to avoid GPU errors in Makie plots for large regions 
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
    #    1      'Ia'                'Strict Nature Reserve'          
    #    2      'Ib'                'Wilderness Area'                
    #    3      'II'                'National Park'                  
    #    4      'III'               'Natural Monument'               
    #    5      'IV'                'Habitat/Species Management'     
    #    6      'V'                 'Protected Landscape/Seascape'   
    #    7      'VI'                'Managed Resource Protected Area'
    #    8      'Not Reported'      'Not Reported'                   
    #    9      'Not Applicable'    'Not Applicable'                 
    #    10     'Not Assigned'      'Not Assigned'                   


mutable struct SolarOptions
    gisregion               ::String
    filenamesuffix          ::String
    pv_density              ::Float64           # W/m2
    csp_density             ::Float64           # W/m2
    pvroof_area             ::Float64           # share [0-1]
    plant_area              ::Float64           # share [0-1]
    distance_elec_access    ::Float64           # km
    plant_persons_per_km2   ::Float64           # persons/km2
    pvroof_persons_per_km2  ::Float64           # persons/km2
    exclude_landtypes       ::Vector{Int}
    protected_codes         ::Vector{Int}
    scenarioyear            ::String
    era_year                ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    pvclasses_min           ::Vector{Float64}
    pvclasses_max           ::Vector{Float64}
    cspclasses_min          ::Vector{Float64}
    cspclasses_max          ::Vector{Float64}
    downsample_masks        ::Int
end

SolarOptions() = SolarOptions("","",0,0,0,0,0,0,0,[],[],"",0,0,0,[],[],[],[],0)

function SolarOptions(d::Dict{Symbol,Any})
    options = SolarOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function GISsolar(; savetodisk=true, plotmasks=false, optionlist...)

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
    @unpack gisregion, era_year, filenamesuffix, pv_density, csp_density, downsample_masks = options

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(options)

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(options, regions, gridaccess, popdens, land, protected, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)

    return nothing  # uncomment to terminate after plotting masks
    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(options, lonrange, latrange)

    CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB, solar_overlap_areaA, solar_overlap_areaB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB =
        calc_solar_vars(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    if savetodisk
        mkpath(in_datafolder("output"))
        matopen(in_datafolder("output", "GISdata_solar$(era_year)_$gisregion$filenamesuffix.mat"), "w") do file
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
            write(file, "solar_overlap_areaA", solar_overlap_areaA)
            write(file, "solar_overlap_areaB", solar_overlap_areaB)
            write(file, "pv_density", pv_density)
            write(file, "csp_density", csp_density)
        end
    end

    nothing
    # return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
    #         capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

function read_solar_datasets(options, lonrange, latrange)
    @unpack res, erares, era_year = options

    println("Reading ERA5 solar datasets...")
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)

    @time meanGTI, solarGTI, meanDNI, solarDNI = h5open(in_datafolder("era5solar$era_year.h5"), "r") do file
        if length(eralonranges) == 1
            file["meanGTI"][eralonranges[1], eralatrange],
                file["GTI"][:,eralonranges[1], eralatrange],
                file["meanDNI"][eralonranges[1], eralatrange],
                file["DNI"][:,eralonranges[1], eralatrange]
        else
            [file["meanGTI"][eralonranges[1], eralatrange]; file["meanGTI"][eralonranges[2], eralatrange]],
                [file["GTI"][:, eralonranges[1], eralatrange] file["GTI"][:, eralonranges[2], eralatrange]],
                [file["meanDNI"][eralonranges[1], eralatrange]; file["meanDNI"][eralonranges[2], eralatrange]],
                [file["DNI"][:, eralonranges[1], eralatrange] file["DNI"][:, eralonranges[2], eralatrange]]
        end
    end
    return meanGTI, solarGTI, meanDNI, solarDNI
end

function create_solar_masks(options, regions, gridaccess, popdens, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, plant_persons_per_km2, pvroof_persons_per_km2 = options

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
    gridB = (imfilter(gridaccess, disk) .> 10^(-1.5))

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    mask_rooftop = gridA .& (popdens .> pvroof_persons_per_km2) .& .!protected_area
    mask_plantA = gridA .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area
    mask_plantB = (gridB .& .!gridA) .& (popdens .< plant_persons_per_km2) .& goodland .& .!protected_area

    if plotmasks
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> plant_persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        masks[(masks .== 0) .& .!gridA .& gridB] .= 6
        masks[(masks .== 0) .& isregion] .= 5
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "no grid", "solar plant A", "solar plant B", "", ""]
        maskmap("$(gisregion)_masks_solar", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    end

    return mask_rooftop, mask_plantA, mask_plantB
end

function increment_solarCF!(cf::AbstractVector{<:AbstractFloat}, solardata::AbstractVector{<:AbstractFloat})
    @inbounds for i = 1:length(cf)
        cf[i] += solardata[i]
    end
end

function makesolarclasses(options, meanGTI, meanDNI)
    println("Allocating pixels to classes using ERA5 annual means...")

    @unpack pvclasses_min, pvclasses_max, cspclasses_min, cspclasses_max = options

    pvclass = getclasses(meanGTI, pvclasses_min, pvclasses_max)
    cspclass = getclasses(meanDNI, cspclasses_min, cspclasses_max)

    return pvclass, cspclass
end

function calc_solar_vars(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    pvclass, cspclass = makesolarclasses(options, meanGTI, meanDNI)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(options, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and class...")
    println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, pvclasses_min, cspclasses_min, res, erares, pv_density, csp_density, pvroof_area, plant_area = options

    numreg = length(regionlist)
    npvclasses, ncspclasses = length(pvclasses_min), length(cspclasses_min)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    capacity_pvrooftop = zeros(numreg,npvclasses)
    capacity_pvplantA = zeros(numreg,npvclasses)
    capacity_pvplantB = zeros(numreg,npvclasses)
    capacity_cspplantA = zeros(numreg,ncspclasses)
    capacity_cspplantB = zeros(numreg,ncspclasses)
    CF_pvrooftop = zeros(yearlength,numreg,npvclasses)
    CF_pvplantA = zeros(yearlength,numreg,npvclasses)
    CF_pvplantB = zeros(yearlength,numreg,npvclasses)
    CF_cspplantA = zeros(yearlength,numreg,ncspclasses)
    CF_cspplantB = zeros(yearlength,numreg,ncspclasses)
    count_pvrooftop = zeros(Int,numreg,npvclasses)
    count_pvplantA = zeros(Int,numreg,npvclasses)
    count_pvplantB = zeros(Int,numreg,npvclasses)
    count_cspplantA = zeros(Int,numreg,ncspclasses)
    count_cspplantB = zeros(Int,numreg,ncspclasses)
    solar_overlap_areaA = zeros(numreg,npvclasses,ncspclasses)
    solar_overlap_areaB = zeros(numreg,npvclasses,ncspclasses)

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanGTI[i,j] == 0 && meanDNI[i,j] == 0 && continue
            GTI = solarGTI[:, i, j]
            DNI = solarDNI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                (reg == 0 || reg == NOREGION) && continue 

                area = cellarea[c]
                class = pvclass[i,j]
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if class > 0
                    if mask_rooftop[r,c] > 0
                        capacity_pvrooftop[reg,class] += 1/1000 * pv_density * pvroof_area * area
                        increment_solarCF!(CF_pvrooftop[:,reg,class], GTI)
                        count_pvrooftop[reg,class] += 1
                    elseif mask_plantA[r,c] > 0
                        capacity_pvplantA[reg,class] += 1/1000 * pv_density * plant_area * area
                        increment_solarCF!(CF_pvplantA[:,reg,class], GTI)
                        count_pvplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_pvplantB[reg,class] += 1/1000 * pv_density * 2 * plant_area * area
                        increment_solarCF!(CF_pvplantB[:,reg,class], GTI)
                        count_pvplantB[reg,class] += 1
                    end
                end

                class_pv = class
                class = cspclass[i,j]
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if class > 0
                    if mask_plantA[r,c] > 0
                        capacity_cspplantA[reg,class] += 1/1000 * csp_density * plant_area * area
                        increment_solarCF!(CF_cspplantA[:,reg,class], DNI)
                        count_cspplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_cspplantB[reg,class] += 1/1000 * csp_density * 2 * plant_area * area
                        increment_solarCF!(CF_cspplantB[:,reg,class], DNI)
                        count_cspplantB[reg,class] += 1
                    end
                end

                if class_pv > 0 && class > 0
                    if mask_plantA[r,c] > 0
                        solar_overlap_areaA[reg,class_pv,class] += 1/1000 * plant_area * area
                    elseif mask_plantB[r,c] > 0
                        solar_overlap_areaB[reg,class_pv,class] += 1/1000 * 2 * plant_area * area
                    end
                end
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

    return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB, solar_overlap_areaA, solar_overlap_areaB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end





# Quick and ugly copy/paste hack to create resource maps for solar classes combined with masks.
function GISsolarmap(; optionlist...)
    options = SolarOptions(merge(solaroptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix = options

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(options)
    meanGTI, solarGTI, meanDNI, solarDNI = read_solar_datasets(options, lonrange, latrange)

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(options, regions, gridaccess, popdens, land, protected, lonrange, latrange, plotmasks=true)

    pvmap, pvrooftopmap, cspmap =
        calc_solar_map(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    return pvmap, pvrooftopmap, cspmap
end

function calc_solar_map(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    pvclass, cspclass = makesolarclasses(options, meanGTI, meanDNI)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(options, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and class...")
    println("Interpolate ERA5 insolation later (maybe 4x runtime).")

    @unpack era_year, pvclasses_min, cspclasses_min, res, erares, pv_density, csp_density, pvroof_area, plant_area = options

    numreg = length(regionlist)
    npvclasses, ncspclasses = length(pvclasses_min), length(cspclasses_min)
    yearlength, nlons, nlats = size(solarGTI)
    firsttime = DateTime(era_year, 1, 1)

    # pvmap = zeros(size(regions))
    # pvrooftopmap = zeros(size(regions))
    # cspmap = zeros(size(regions))
    pvmap = zeros(Int16, size(regions))
    pvrooftopmap = zeros(Int16, size(regions))
    cspmap = zeros(Int16, size(regions))

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanGTI[i,j] == 0 && meanDNI[i,j] == 0 && continue
            GTI = solarGTI[:, i, j]
            DNI = solarDNI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                area = cellarea[c]

                class = pvclass[i,j]
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0
                    if mask_rooftop[r,c] > 0
                        # pvrooftopmap[r,c] = meanGTI[i,j]
                        pvrooftopmap[r,c] = class
                    elseif mask_plantA[r,c] > 0
                        # pvmap[r,c] = meanGTI[i,j]
                        pvmap[r,c] = class
                    elseif mask_plantB[r,c] > 0
                        # pvmap[r,c] = class
                    end
                end

                class = cspclass[i,j]
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0
                    if mask_plantA[r,c] > 0
                        # cspmap[r,c] = meanGTI[i,j]
                        cspmap[r,c] = class
                    elseif mask_plantB[r,c] > 0
                        # cspmap[r,c] = class
                    end
                end
            end
        end
        next!(updateprogress)
    end

    return pvmap, pvrooftopmap, cspmap
end