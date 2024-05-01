windoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat" 

    :onshore_density => 5,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 8,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => .08,               # area available for onshore wind power after the masks have been applied
    :area_offshore => .33,              # area available for offshore wind power after the masks have been applied

    :distance_elec_access => 150,       # max distance to grid [km] (for wind classes of category B and offshore)
    :persons_per_km2 => 150,            # not too crowded, max X persons/km2
                                        # US census bureau requires 1000 ppl/mile^2 = 386 ppl/km2 for "urban" (half in Australia)
                                        # roughly half the people of the world live at density > 300 ppl/km2
    :max_depth => 40,                   # max depth for offshore wind [m]
    :min_shore_distance => 5,           # minimum distance to shore for offshore wind [km]
    :exclude_landtypes => [0,11,13],    # exclude water, wetlands and urban areas. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 
    :rescale_to_wind_atlas => true,     # rescale the ERA5 time series to fit annual wind speed averages from the Global Wind Atlas
    :wind_speed_altitude => 100,        # (m) altitude of GWA wind speeds that determine hourly power output & potentials 
    :wind_class_altitude => 100,        # (m) altitude of GWA wind speeds for selection of wind classes 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :onshoreclasses_min => [2,5,6,7,8],     # lower bound on annual onshore wind speeds for class X    [0:0.25:12.25;]
    :onshoreclasses_max => [5,6,7,8,99],    # upper bound on annual onshore wind speeds for class X    [0.25:0.25:12.5;]
    :offshoreclasses_min => [3,6,7,8,9],    # lower bound on annual offshore wind speeds for class X
    :offshoreclasses_max => [6,7,8,9,99],   # upper bound on annual offshore wind speeds for class X

    :grid_everywhere => false,  # set to true to assume all pixels have grid access
    :downsample_masks => 1,     # set to 2 or higher to scale down mask sizes to avoid GPU errors in Makie plots for large regions 
    :classB_threshold => 0.001, # minimum share of pixels within distance_elec_access km that must have grid access
                                # for a pixel to be considered for wind class B.

    :climate_scenario => ""     # e.g. "HCLIM_EC-EARTH_100m_rcp85_2050", "CORDEX_ictp_EC-EARTH_100m_rcp85_2050"
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

mutable struct WindOptions
    gisregion               ::String
    filenamesuffix          ::String
    onshore_density         ::Float64           # W/m2
    offshore_density        ::Float64           # W/m2
    area_onshore            ::Float64           # share [0-1]
    area_offshore           ::Float64           # share [0-1]
    distance_elec_access    ::Float64           # km
    persons_per_km2         ::Float64           # persons/km2
    max_depth               ::Float64           # m
    min_shore_distance      ::Float64           # km
    exclude_landtypes       ::Vector{Int}
    protected_codes         ::Vector{Int}
    scenarioyear            ::String
    era_year                ::Int
    rescale_to_wind_atlas   ::Bool
    wind_speed_altitude     ::Int
    wind_class_altitude     ::Int
    res                     ::Float64           # degrees/pixel
    erares                  ::Float64           # degrees/pixel
    onshoreclasses_min      ::Vector{Float64}
    onshoreclasses_max      ::Vector{Float64}
    offshoreclasses_min     ::Vector{Float64}
    offshoreclasses_max     ::Vector{Float64}
    grid_everywhere         ::Bool
    downsample_masks        ::Int
    classB_threshold        ::Float64
    climate_scenario        ::String
end

WindOptions() = WindOptions("","",0,0,0,0,0,0,0,0,[],[],"",0,false,100,100,0,0,[],[],[],[],false,0,0.0,"")

function WindOptions(d::Dict{Symbol,Any})
    options = WindOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function GISwind(; savetodisk=true, plotmasks=false, optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix, downsample_masks, climate_scenario = options

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(options)

    mask_onshoreA, mask_onshoreB, mask_offshore =
        create_wind_masks(options, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange,
                            plotmasks=plotmasks, downsample=downsample_masks)

    plotmasks == :onlymasks && return nothing

    windatlas, windatlas_class, meanwind, windspeed, meanwind_allyears = read_wind_datasets(options, lonrange, latrange)

    windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore =
        calc_wind_vars(options, windatlas, windatlas_class, meanwind, windspeed, meanwind_allyears, regions,
                offshoreregions, regionlist, mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    if savetodisk
        mkpath(in_datafolder("output"))
        suffix = isempty(climate_scenario) ? filenamesuffix : "_$climate_scenario$filenamesuffix"
        matopen(in_datafolder("output", "GISdata_wind$(era_year)_$gisregion$suffix.mat"), "w", compress=true) do file
            write(file, "CFtime_windonshoreA", windCF_onshoreA)
            write(file, "CFtime_windonshoreB", windCF_onshoreB)
            write(file, "CFtime_windoffshore", windCF_offshore)
            write(file, "capacity_onshoreA", capacity_onshoreA)
            write(file, "capacity_onshoreB", capacity_onshoreB)
            write(file, "capacity_offshore", capacity_offshore)
        end
    end

    nothing
    # return windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore
end

function read_datasets(options)
    @unpack res, gisregion, scenarioyear = options

    println("\nReading auxiliary datasets...")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    cellarea = rastercellarea.(lats, res)

    gridaccess = JLD.load(in_datafolder("gridaccess_$scenarioyear.jld"), "gridaccess")[lonrange,latrange]
    pop = JLD.load(in_datafolder("population_$scenarioyear.jld"), "population")[lonrange,latrange]
    topo = JLD.load(in_datafolder("topography.jld"), "topography")[lonrange,latrange]
    land = JLD.load(in_datafolder("landcover.jld"), "landcover")[lonrange,latrange]
    protected = JLD.load(in_datafolder("protected.jld"), "protected")[lonrange,latrange]

    popdens = pop ./ cellarea'

    # drawmap(log.(pop))

    return regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange
end

function read_wind_datasets(options, lonrange, latrange)
    @unpack res, erares, era_year, wind_speed_altitude, wind_class_altitude, climate_scenario = options

    println("Reading ERA5 wind speeds and calculating capacity factors...")
    if isempty(climate_scenario)
        eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
        @time meanwind, windspeed = h5open(in_datafolder("era5wind$era_year.h5"), "r") do file
            if length(eralonranges) == 1
                file["meanwind"][eralonranges[1], eralatrange],
                    file["wind"][:, eralonranges[1], eralatrange]
            else
                [file["meanwind"][eralonranges[1], eralatrange]; file["meanwind"][eralonranges[2], eralatrange]],
                    [file["wind"][:, eralonranges[1], eralatrange] file["wind"][:, eralonranges[2], eralatrange]]
            end
        end
        annualwind = getmonthlywind(:annual, :wind, eralonranges, eralatrange, "")
        meanwind_allyears = meandrop(annualwind, dims=1)
    else
        climatefolder = getconfig("climatefolder")
        @time meanwind, windspeed = h5open("$climatefolder/wind_$(climate_scenario).h5", "r") do file
            file["meanwind"][:,:], permutedims(file["wind"][:,:,:], (3,1,2))
        end
        climyear = Base.parse(Int, climate_scenario[end-3:end])
        meanyears = climyear < 2030 ? "1996_2005" : climyear < 2070 ? "2046_2055" : "2091_2100"
        datasource = Base.contains(climate_scenario, "HCLIM") ? "HCLIM" : "CORDEX"
        meanwind_allyears = h5read("$climatefolder/meanwind_$(meanyears)_$(datasource)_EC-EARTH_100m.h5", "/meanwind")
    end

    windatlas = getwindatlas(wind_speed_altitude)[lonrange,latrange]

    if wind_speed_altitude == wind_class_altitude
        return windatlas, windatlas, meanwind, windspeed, meanwind_allyears
    else
        windatlas_class = getwindatlas(wind_class_altitude)[lonrange,latrange]
        return windatlas, windatlas_class, meanwind, windspeed, meanwind_allyears
    end
end

function create_wind_masks(options, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange; plotmasks=false, downsample=1)
    @unpack res, gisregion, exclude_landtypes, protected_codes, distance_elec_access, persons_per_km2,
                min_shore_distance, max_depth, classB_threshold, filenamesuffix, grid_everywhere = options

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
    gridA = grid_everywhere ? fill(true, size(gridaccess)) : (gridaccess .> 0)

    # Pixels with electricity access for onshore wind B and offshore wind
    km_per_degree = π*2*6371/360
    disk = diskfilterkernel(distance_elec_access/km_per_degree/res)
    gridB = (imfilter(gridaccess, disk) .> max(1e-9, classB_threshold)) # avoid artifacts if classB_threshold == 0

    # println("MAKE SURE MASKS DON'T OVERLAP! (regions & offshoreregions, mask_*)")

    # all mask conditions
    mask_onshoreA = gridA .& (popdens .< persons_per_km2) .& goodland .& .!protected_area
    mask_onshoreB = (gridB .& .!gridA) .& (popdens .< persons_per_km2) .& goodland .& .!protected_area

    # shoreline mask for offshore wind
    disk = diskfilterkernel(min_shore_distance/km_per_degree/res)
    shore = (imfilter(regions .> 0, disk) .> 1e-6)

    # all mask conditions
    mask_offshore = gridB .& .!shore .& (topo .> -max_depth) .& (offshoreregions .> 0) .& .!protected_area

    if plotmasks != false   # can == :onlymasks as well
        # drawmap(land)
        isregion = (regions .> 0) .& (regions .!= NOREGION)

        # mask values refer to colors in ColorBrewer Set2_7:
        # https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#colorbrewer-1
        masks = zeros(Int16, size(regions))
        masks[(masks .== 0) .& (popdens .> persons_per_km2)] .= 2
        masks[(masks .== 0) .& protected_area] .= 3
        masks[(masks .== 0) .& .!gridA .& .!gridB] .= 4
        masks[(masks .== 0) .& .!goodland] .= 1
        masks[(masks .== 0) .& .!gridA .& gridB] .= 8
        masks[(masks .== 0) .& isregion] .= 7
        masks[regions .== 0] .= 0
        masks[regions .== NOREGION] .= NOREGION
        legendtext = ["bad land type", "high population", "protected area", "no grid", "", "", "wind plant A", "wind plant B"]
        maskmap("$(gisregion)_masks_wind$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)

        isregion = (offshoreregions .> 0) .& (offshoreregions .!= NOREGION)
        masks = zeros(Int16, size(offshoreregions))
        masks[(masks .== 0) .& shore .& isregion] .= 1
        masks[(masks .== 0) .& protected_area] .= 3
        masks[(masks .== 0) .& (topo .<= -max_depth)] .= 8
        masks[(masks .== 0) .& .!gridB] .= 4
        masks[(masks .== 0) .& isregion] .= 7
        masks[offshoreregions .== 0] .= NOREGION
        masks[offshoreregions .== NOREGION] .= 0
        legendtext = ["near shore", "", "protected area", "no grid", "", "", "wind offshore", "too deep water"]
        maskmap("$(gisregion)_masks_windoffshore$filenamesuffix", masks, legendtext, lonrange, latrange; legend=true, downsample=downsample)
    end

    return mask_onshoreA, mask_onshoreB, mask_offshore
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

function increment_windCF!(cf::AbstractVector{<:AbstractFloat}, speed_or_cf::AbstractVector{<:AbstractFloat}, factor, rescale::Bool)
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

function getclasses0(windatlas, classes_min, classes_max)
    class = zeros(Int16, size(windatlas))
    for c = 1:length(classes_min)
        class[(windatlas .>= classes_min[c]) .& (windatlas .< classes_max[c])] .= c
    end
    return class
end

function getclasses(windatlas, classes_min, classes_max)
    class = similar(windatlas, Int16)
    nclasses = length(classes_min)
    @inbounds for i = 1:length(windatlas)
        wa = windatlas[i]
        if wa < classes_min[1] || wa > classes_max[end]
            class[i] = 0
            continue
        end
        for c = 1:nclasses
            if wa < classes_max[c]
                class[i] = c
                break
            end
        end
    end
    return class
end

function makewindclasses(options, windatlas)
    # println("CHANGE TO CAPACITY FACTOR LATER!")

    @unpack onshoreclasses_min, onshoreclasses_max, offshoreclasses_min, offshoreclasses_max = options

    onshoreclass = getclasses(windatlas, onshoreclasses_min, onshoreclasses_max)
    offshoreclass = getclasses(windatlas, offshoreclasses_min, offshoreclasses_max)

    return onshoreclass, offshoreclass
end

function calc_wind_vars(options, windatlas, windatlas_class, meanwind, windspeed, meanwind_allyears, regions,
                offshoreregions, regionlist, mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore, climate_scenario = options
    
    numreg = length(regionlist)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    capacity_onshoreA = zeros(numreg,nonshoreclasses)
    capacity_onshoreB = zeros(numreg,nonshoreclasses)
    capacity_offshore = zeros(numreg,noffshoreclasses)
    windCF_onshoreA = zeros(yearlength,numreg,nonshoreclasses)
    windCF_onshoreB = zeros(yearlength,numreg,nonshoreclasses)
    windCF_offshore = zeros(yearlength,numreg,noffshoreclasses)
    count_onshoreA = zeros(Int,numreg,nonshoreclasses)
    count_onshoreB = zeros(Int,numreg,nonshoreclasses)
    count_offshore = zeros(Int,numreg,noffshoreclasses)

    if Base.contains(climate_scenario, "HCLIM")
        erares = 0.03
        extent = [2, 50.5, 32, 71.5]
        # meanwind_allyears = zeros(size(meanwind))
        GWAclasses = false
        windatlas_class = meanwind_allyears
    elseif Base.contains(climate_scenario, "CORDEX")
        erares = 0.1
        extent = [-10.5, 34.5, 32, 71.5]
        # meanwind_allyears = zeros(size(meanwind))
        GWAclasses = false
        windatlas_class = meanwind_allyears
    else
        erares = 0.28125
        eralons, eralats, _, _, _ = eralonlat(options, lonrange, latrange)
        extent = [eralons[1]-erares/2, eralats[end]-erares/2, eralons[end]+erares/2, eralats[1]+erares/2]
        GWAclasses = true
        if !rescale_to_wind_atlas
            GWAclasses = false
            windatlas_class = meanwind_allyears
        end
    end

    println("Calculating GW potential and hourly capacity factors for each region and wind class...")
    println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")
    println("Allocating pixels to classes using the Global Wind Atlas...")
    onshoreclass, offshoreclass = makewindclasses(options, windatlas_class)

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        # println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    @assert size(meanwind_allyears) == size(meanwind)
    meanwindGeo = GeoArray(meanwind, erares, extent)

    lons = (-180+res/2:res:180-res/2)[lonrange]         # longitude values (pixel center)
    lats = (90-res/2:-res:-90+res/2)[latrange]          # latitude values (pixel center)
    lonlim = (lons[1]-res/2, lons[end]+res/2)
    latlim = (lats[end]-res/2, lats[1]+res/2)
    regionsGeo = GeoArray(regions, res, lonlim, latlim)

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    for j in randperm(nlats)
        eralat = getlat(meanwindGeo, j)
        colrange = lat_indices_within(regionsGeo, eralat - erares/2, eralat + erares/2)
        for i = 1:nlons
            meanwind[i,j] == 0 && continue
            wind = rescale_to_wind_atlas ? windspeed[:, i, j] : speed2capacityfactor.(windspeed[:, i, j])
            eralon = getlon(meanwindGeo, i)
            rowrange = lon_indices_within(regionsGeo, eralon - erares/2, eralon + erares/2)
            
            # for all high resolution row and column indexes within this ERA5 cell         
            for c in colrange, r in rowrange
                reg = regions[r,c]
                offreg = offshoreregions[r,c]
                area = rastercellarea(getlat(regionsGeo, c), res)

                # @views is needed to make sure increment_windCF!() works with matrix
                # slice. It's also faster since it avoids making copies.
                # We divide by meanwind_allyears instead of meanwind in order to
                # account for more and less windy years, while still rescaling 
                # hourly ERA5 wind speeds to match annual averages from the Global
                # Wind Atlas.
                if reg > 0 && reg != NOREGION
                    class = GWAclasses ? onshoreclass[r,c] : onshoreclass[i,j]
                    @views if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                        capacity_onshoreA[reg,class] += 1/1000 * onshore_density * area_onshore * area
                        increment_windCF!(windCF_onshoreA[:,reg,class], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                        count_onshoreA[reg,class] += 1
                    elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                        capacity_onshoreB[reg,class] += 1/1000 * onshore_density * 2 * area_onshore * area
                        increment_windCF!(windCF_onshoreB[:,reg,class], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
                        count_onshoreB[reg,class] += 1
                    end
                end
                if offreg > 0 && offreg != NOREGION
                    offclass = GWAclasses ? offshoreclass[r,c] : offshoreclass[i,j]
                    @views if offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                        capacity_offshore[offreg,offclass] += 1/1000 * offshore_density * area_offshore * area
                        increment_windCF!(windCF_offshore[:,offreg,offclass], wind, windatlas[r,c] / meanwind_allyears[i,j], rescale_to_wind_atlas)
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

    return windCF_onshoreA, windCF_onshoreB, windCF_offshore, capacity_onshoreA, capacity_onshoreB, capacity_offshore
end

# returns a wind speed index for 1979-2019, i.e. normalized average wind speed over
# each pixel in each region for each year
# also return a matrix of annual wind speeds in each ERA5 pixel over the region
function annualwindindex(options; resource=:wind, sites_quantile=1.0, aggregateregions=[])
    @unpack gisregion, era_year, res, erares, filenamesuffix = options
    regions, _, regionlist, lonrange, latrange = loadregions(gisregion)
    smallregions = resize_categorical(regions, regionlist, lonrange, latrange;
                skipNOREGION=true)
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    varname = (resource == :wind) ? "annualwind" : "annualssrd"
    annualwind = getmonthlywind(:annual, resource, eralonranges, eralatrange, filenamesuffix)
    meanwind = meandrop(annualwind, dims=1)
    nyears, nreg = size(annualwind, 1), length(regionlist)
    aggregateregions = isempty(aggregateregions) ? [[r] for r = 1:nreg] : aggregateregions
    regionmasks = [in.(smallregions, Ref(aggreg)) for aggreg in aggregateregions]
    cutoffs = [quantile(vec(meanwind[(meanwind.>0) .& regmask]), 1-sites_quantile)
            for regmask in regionmasks]
    masks = [(meanwind .> cutoffs[i]) .& regionmasks[i] for i = 1:length(regionmasks)]
    index = [sum(annualwind[y,:,:] .* mask) ./ sum(mask) for y = 1:nyears, mask in masks]
    return index, index ./ mean(index, dims=1), meanwind, annualwind[era_year-1979+1,:,:]
end

function annualwindindex(; resource=:wind, sites_quantile=1.0, aggregateregions=[], optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    return annualwindindex(options, resource=resource, sites_quantile=sites_quantile,
                    aggregateregions=aggregateregions)
end

function seasonalwindprofile(; resource=:wind, sites_quantile=1.0, aggregateregions=[], optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, res, erares, filenamesuffix = options
    regions, _, regionlist, lonrange, latrange = loadregions(gisregion)
    smallregions = resize_categorical(regions, regionlist, lonrange, latrange;
                skipNOREGION=true)
    eralonranges, eralatrange = eraranges(lonrange, latrange, res, erares)
    annualwind = getmonthlywind(:annual, resource, eralonranges, eralatrange, filenamesuffix)
    monthlywind = getmonthlywind(:monthly, resource, eralonranges, eralatrange, filenamesuffix)
    meanwind = meandrop(annualwind, dims=1)
    cutoff = quantile(vec(meanwind[meanwind.>0]), 1-sites_quantile)
    nyears, nreg = size(monthlywind, 1), length(regionlist)
    aggregateregions = isempty(aggregateregions) ? [[r] for r = 1:nreg] : aggregateregions
    profile = [sum(meandrop(monthlywind[m:12:end,:,:], dims=1) .*
        (mask = (meanwind .> cutoff) .& in.(smallregions, Ref(aggreg)))) ./ sum(mask)
                    for m = 1:12, aggreg in aggregateregions]
end

function getmonthlywind(time, resource, eralonranges, eralatrange, filenamesuffix)
    h5open(in_datafolder("era5monthly$resource$filenamesuffix.h5"), "r") do file
        varname = (resource == :wind) ? "$(time)wind" : "$(time)ssrd"
        if length(eralonranges) == 1
            file[varname][:, eralonranges[1], eralatrange]
        else
            [file[varname][:, eralonranges[1], eralatrange] file[varname][:, eralonranges[2], eralatrange]]
        end
    end
end

# Calculate annual mean wind speeds and capacity factors for a list of locations,
# given vectors of lons and lats.
function calcCF(lons, lats; optionlist...)
    len = length(lons)
    @assert len == length(lats)
    options = WindOptions(merge(windoptions(), optionlist))
    windatlas, _, meanwind, windspeed, _ = read_wind_datasets(options, 1:36000, 1:18000)
    annualwind = getmonthlywind(:annual, :wind, [1:1280], 1:640, "")
    meanwind_allyears = meandrop(annualwind, dims=1)

    res, erares = 0.01, 0.28125
    windatlasGeo = GeoArray(windatlas, res)
    meanwindGeo = GeoArray(meanwind, erares)

    speeds, cf = zeros(len), zeros(len)
    for i = 1:len
        lon, lat = lons[i], lats[i]
        index = lonlat_index(windatlasGeo, lon, lat)
        eraindex = lonlat_index(meanwindGeo, lon, lat)
        factor = windatlas[index] / meanwind_allyears[eraindex]
        speed = windspeed[:, eraindex] * factor
        speeds[i] = mean(speed)
        cf[i] = mean(speed2capacityfactor.(speed))
    end
    return speeds, cf
end

function calcCF_all(lons, lats, years=2008:2020; optionlist...)
    len = length(lons)
    speeds, cf = zeros(len), zeros(len)
    for y in years
        println("\nYear $y:")
        _speeds, _cf = calcCF(lons, lats; optionlist..., era_year=y)
        speeds += _speeds
        cf += _cf
    end
    speeds /= length(years)
    cf /= length(years)
    return speeds, cf
end

# Quick and ugly copy/paste hack to create resource maps for wind classes combined with masks.
function GISwindmap(; optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    @unpack gisregion, era_year, filenamesuffix = options

    regions, offshoreregions, regionlist, gridaccess, popdens, topo, land, protected, lonrange, latrange =
                read_datasets(options)
    _, windatlas_class, meanwind, windspeed = read_wind_datasets(options, lonrange, latrange)

    mask_onshoreA, mask_onshoreB, mask_offshore =
        create_wind_masks(options, regions, offshoreregions, gridaccess, popdens, topo, land, protected, lonrange, latrange, plotmasks=true)

    onshoremap, offshoremap = calc_wind_map(options, windatlas_class, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    return onshoremap, offshoremap
end

function calc_wind_map(options, windatlas_class, meanwind, windspeed, regions, offshoreregions, regionlist,
                mask_onshoreA, mask_onshoreB, mask_offshore, lonrange, latrange)

    println("Calculating GW potential and hourly capacity factors for each region and wind class...")
    # println("Interpolate ERA5 wind speeds later (maybe 4x runtime).")

    onshoreclass, offshoreclass = makewindclasses(options, windatlas_class)
    eralons, eralats, lonmap, latmap, cellarea = eralonlat(options, lonrange, latrange)

    @unpack onshoreclasses_min, offshoreclasses_min, rescale_to_wind_atlas, res, erares,
                onshore_density, area_onshore, offshore_density, area_offshore = options

    numreg = length(regionlist)
    nonshoreclasses, noffshoreclasses = length(onshoreclasses_min), length(offshoreclasses_min)
    yearlength, nlons, nlats = size(windspeed)

    onshoremap = zeros(Int16, size(regions))
    offshoremap = zeros(Int16, size(regions))

    if rescale_to_wind_atlas
        println("\nRescaling ERA5 wind speeds to match annual wind speeds from the Global Wind Atlas.")
        println("This will increase run times by an order of magnitude (since GWA has very high spatial resolution).")
    end

    # Run times vary wildly depending on geographical area (because of far offshore regions with mostly zero wind speeds).
    # To improve the estimated time of completing the progress bar, iterate over latitudes in random order.
    Random.seed!(1)
    updateprogress = Progress(nlats, 1)
    @inbounds for j in randperm(nlats)
        eralat = eralats[j]
        colrange = latmap[lat2col(eralat+erares/2, res):lat2col(eralat-erares/2, res)-1]
        for i = 1:nlons
            meanwind[i,j] == 0 && continue
            wind = rescale_to_wind_atlas ? windspeed[:, i, j] : speed2capacityfactor.(windspeed[:, i, j])
            eralon = eralons[i]
            # get all high resolution row and column indexes within this ERA5 cell         
            rowrange = lonmap[lon2row(eralon-erares/2, res):lon2row(eralon+erares/2, res)-1]

            for c in colrange, r in rowrange
                (c == 0 || r == 0) && continue
                reg = regions[r,c]
                area = cellarea[c]
                class = onshoreclass[r,c]
                offreg = offshoreregions[r,c]
                offclass = offshoreclass[r,c]
                
                # can't use elseif here, probably some overlap in the masks
                # @views is needed to make sure increment_windCF!() works with matrix slices
                # also faster since it avoids making copies
                @views if reg > 0 && class > 0 && mask_onshoreA[r,c] > 0
                    onshoremap[r,c] = class
                elseif reg > 0 && class > 0 && mask_onshoreB[r,c] > 0
                    # onshoremap[r,c] = class
                elseif offreg > 0 && offclass > 0 && mask_offshore[r,c] > 0
                    offshoremap[r,c] = class
                end
            end
        end
        next!(updateprogress)
    end

    return onshoremap, offshoremap
end
