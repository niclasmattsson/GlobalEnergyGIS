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

    :scenario => "ssp2_2050",           # default scenario for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.08,0.14,0.18,0.22,0.26],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49;]
    :pvclasses_max => [0.14,0.18,0.22,0.26,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50;]
    :cspclasses_min => [0.10,0.18,0.24,0.28,0.32],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.18,0.24,0.28,0.32,1.00]  # upper bound on annual CSP capacity factor for class X
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
    scenario                ::String
    era_year                ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    pvclasses_min           ::Vector{Float64}
    pvclasses_max           ::Vector{Float64}
    cspclasses_min          ::Vector{Float64}
    cspclasses_max          ::Vector{Float64}
end

SolarOptions() = SolarOptions("",0,0,0,0,0,0,0,[],[],"",0,0,0,[],[],[],[])

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

    mask_rooftop, mask_plantA, mask_plantB =
        create_solar_masks(options, regions, gridaccess, pop, land, protected)

    CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB =
        calc_solar_vars(options, meanGTI, solarGTI, meanDNI, solarDNI, regions, offshoreregions, regionlist,
                mask_rooftop, mask_plantA, mask_plantB, lonrange, latrange)

    matopen("GISdata_solar$(options.era_year)_$(options.gisregion).mat", "w") do file
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

    return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
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
                [file["GTI"][:, eralonranges[1], eralatrange] file["GTI"][:, eralonranges[2], eralatrange]],
                [file["meanDNI"][eralonranges[1], eralatrange]; file["meanDNI"][eralonranges[2], eralatrange]],
                [file["DNI"][:, eralonranges[1], eralatrange] file["DNI"][:, eralonranges[2], eralatrange]]
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

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2+res/5, res)]
        for i = 1:nlons
            meanGTI[i,j] == 0 && meanDNI[i,j] == 0 && continue
            GTI = solarGTI[:, i, j]
            DNI = solarDNI[:, i, j]
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2-res/5, res)]

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
                        capacity_pvrooftop[reg,class] += 1/1000 * pv_density * pvroof_area * area
                        increment_solarCF!(CF_pvrooftop[:,reg,class], GTI)
                        count_pvrooftop[reg,class] += 1
                    elseif mask_plantA[r,c] > 0
                        capacity_pvplantA[reg,class] += 1/1000 * pv_density * plant_area * area
                        increment_solarCF!(CF_pvplantA[:,reg,class], GTI)
                        count_pvplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_pvplantB[reg,class] += 1/1000 * pv_density * plant_area * area
                        increment_solarCF!(CF_pvplantB[:,reg,class], GTI)
                        count_pvplantB[reg,class] += 1
                    end
                end

                class = cspclass[i,j]
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0
                    if mask_plantA[r,c] > 0
                        capacity_cspplantA[reg,class] += 1/1000 * csp_density * plant_area * area
                        increment_solarCF!(CF_cspplantA[:,reg,class], DNI)
                        count_cspplantA[reg,class] += 1
                    elseif mask_plantB[r,c] > 0
                        capacity_cspplantB[reg,class] += 1/1000 * csp_density * plant_area * area
                        increment_solarCF!(CF_cspplantB[:,reg,class], DNI)
                        count_cspplantB[reg,class] += 1
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

    return CF_pvrooftop, CF_pvplantA, CF_pvplantB, CF_cspplantA, CF_cspplantB,
            capacity_pvrooftop, capacity_pvplantA, capacity_pvplantB, capacity_cspplantA, capacity_cspplantB
end

