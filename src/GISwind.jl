function GISwind()

    GISREGION = "Europe10"     # 'global', 'europe8/10/12/15/17', 'china6', 'eurochine14', 'eurasia38/21' or 'mena'

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

    path = joinpath(dirname(@__FILE__), "..")
    file = h5open("C:/Stuff/NewGIS/windatlas.h5", "r")
    windatlas = read(file, "windatlas5min")
    close(file)

    matlabvars = matread("C:/Stuff/NewGIS/testregions_$(GISREGION)_12.mat")
    regionlist = Symbol.(vec(matlabvars["regionlist"]))
    numreg = length(regionlist)
    regions = matlabvars["regions"]
    regions[regions .== numreg+1] .= 0
    # R = matlabvars["R"]

    matlabvars = matread("C:/Stuff/NewGIS/offshoreregions_$(GISREGION)_12.mat")
    offshoreregions = matlabvars["offshoreregions"]
    offshoreregions[offshoreregions .== numreg+1] .= 0

    matlabvars = matread("C:/Stuff/NewGIS/testregions_$(GISREGION)_3.5556.mat")
    smallregions = matlabvars["regions"]
    smallregions[smallregions .== numreg+1] .= 0

    matlabvars = matread("C:/Stuff/NewGIS/offshoreregions_$(GISREGION)_3.5556.mat")
    smalloffshoreregions = matlabvars["offshoreregions"]
    smalloffshoreregions[smalloffshoreregions .== numreg+1] .= 0

    bboxglobal = [-90 -180; 90 180]
    if GISREGION[1:6] == "Europe"
        bbox = [34 -11; 72 32]
    elseif GISREGION == "MENA"
        bbox = [10 -15; 45 70]
    else
        bbox = bboxglobal
    end
    bbox = bboxglobal
    bboxsmall = roundbbox(bbox, 32/9)

    latindex, lonindex = bbox2ranges(bbox, 12)
    regions = regions[latindex, lonindex]
    offshoreregions = offshoreregions[latindex, lonindex]
    latindex, lonindex = bbox2ranges(bboxsmall, 32/9)
    smallregions = smallregions[latindex, lonindex]
    smalloffshoreregions = smalloffshoreregions[latindex, lonindex]

    matlabvars = matread("C:/Stuff/NewGIS/testgrid12.mat")
    gridaccess = matlabvars["gridaccess"]

    matlabvars = matread("C:/Stuff/NewGIS/testpop12.mat")
    pop = matlabvars["pop"]

    matlabvars = matread("C:/Stuff/NewGIS/testtopo12.mat")
    topo = matlabvars["topo"]

    matlabvars = matread("C:/Stuff/NewGIS/testland12.mat")
    land = matlabvars["land"]

    matlabvars = matread("C:/Stuff/NewGIS/testprotected12.mat")
    protected = matlabvars["protected"]

    lats = getlats(bboxglobal, 12, true)
    pixelarea = (2*6371*pi/(360*60/5))^2        # area in km2 of 5 min pixel at equator
    areamatkm = cosd.(lats) * pixelarea         # area of a 5 min grid cell in km2

    yearlength = 8760 + 24*leapyear(ERA_YEAR)

    println("\nReading ERA5 wind capacity factor timeseries...\n")

    latrange, lonrange = bbox2ranges(bboxsmall, 32/9)
    CFtype = RESCALE_ERA_TO_WIND_ATLAS ? "RCF" : "CF"
    filename = "D:/datasets/era5/era5wind$CFtype$ERA_YEAR.h5"
    dataset = "wind$CFtype"
    # if bbox == bboxglobal
    #     @time windCF = h5read_fast(filename, dataset)
    # else
    #     @time windCF = h5read(filename, dataset, (1:yearlength, latrange, lonrange))
    # end
    @time windCF = h5read(filename, dataset, (1:yearlength, latrange, lonrange))
    meanwind = h5read("D:/datasets/era5/era5wind$ERA_YEAR.h5", "meanwind", (latrange, lonrange))

    println("\nInterpolating ERA5 posting data to cell data...\n")

    lat0 = filterrange(getlats(bboxglobal, 32/9, false), bboxsmall[:,1])[1:end-1]
    lon0 = filterrange(getlons(bboxglobal, 32/9, false), bboxsmall[:,2])[1:end-1]
    latrangesmall = getlats(bboxsmall, 32/9, true)
    lonrangesmall = getlons(bboxsmall, 32/9, true)

    meanwind2 = shifthalfcell(meanwind);
    shiftallcells!(windCF)

    return windCF, lat0, lon0, meanwind, meanwind2, latrangesmall, lonrangesmall, latrange, lonrange

end




shifthalfcell(data; kernel=centered([.25 .25; .25 .25])) = imfilter(data, kernel)

function shiftallcells!(windCF)
    yearlength = size(windCF, 1)
    updateprogress = Progress(yearlength, 1)
    kernel = centered([.25 .25; .25 .25])
    for i=1:yearlength
        windCF[i,:,:] = imfilter(windCF[i,:,:], kernel) #shifthalfcell(windCF[i,:,:], kernel)
        next!(updateprogress)
    end
end

function h5read_fast(filename::String, varname::String)
    h5open(filename, "r") do file
        return read_fast(file[varname], Float32)
    end
end

function read_fast(dset::HDF5.HDF5Dataset, T::DataType)
    filetype = HDF5.datatype(dset) # packed layout on disk
    memtype_id = HDF5.h5t_get_native_type(filetype.id) # padded layout in memory
    @assert sizeof(T) == HDF5.h5t_get_size(memtype_id) "Type sizes don't match!"
    out = Array{Float64, length(size(dset))}(undef, size(dset))
    HDF5.h5d_read(dset.id, memtype_id, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, out)
    HDF5.h5t_close(memtype_id)
    out
end


#=

[latindex,lonindex] = bbox2ranges(bbox, 12)
windatlassmall = imresize(windatlas(latindex,lonindex), size(meanwind), 'bicubic')
#windatlassmall = imresize(windatlas, size(meanwind), 'bicubic')
toctic

windatlas = windatlas(latindex,lonindex)
topo = topo(latindex,lonindex)
protected = protected(latindex,lonindex)
pop = pop(latindex,lonindex)
land = land(latindex,lonindex)
gridaccess = gridaccess(latindex,lonindex)
areamatkm = areamatkm(latindex,lonindex)

nlats = size(regions,1)
nlons = size(regions,2)
nlatssmall = size(smallregions,1)
nlonssmall = size(smallregions,2)

# ***** Create masks for onshore wind ******

goodland = single(regions > 0)
for i = EXCLUDE_LANDTYPES
    goodland(land == i) = 0
end
protected_area = single(zeros(size(protected)))
for i = PROTECTED_CODES
    protected_area(protected == i) = 1
end


# Onshore wind A elec access within 9 km, hardcoded by data resolution
gridA = single(gridaccess > 0.1)

# Onshore wind B elec access > 9 km and < 9*ELEC_ACCESS_PIXELS km
disk = fspecial('disk', ELEC_ACCESS_PIXELS)
gridB = single(filter2(disk, gridaccess) > 0.1)

# Onshore wind A elec access within 2*9 km
#disk = fspecial('disk', 2)
#gridA = single(filter2(disk, gridaccess) > 0.1)

# Onshore wind B elec access > 2*9 km and < 9*ELEC_ACCESS_PIXELS km
#disk = fspecial('disk', ELEC_ACCESS_PIXELS)
#gridB = single(filter2(disk, gridaccess) > 0.1)

# all mask conditions
mask_onshoreA = gridA & (pop < PERSONS_PER_KM2) & goodland & ~protected_area
mask_onshoreB = (gridB & ~gridA) & (pop < PERSONS_PER_KM2) & goodland & ~protected_area

if SHOW_FIGURES == 2
    newmap meshm(gridA, Rgrid) title('Direct access to electricity grid')
    newmap meshm(gridB, Rgrid) title('Access to electricity within 150 km')
    gridBonly = single(gridB & ~gridA)
    newmap meshm(gridBonly, Rgrid) title('Access to electricity within 150 km')
    newmap meshm(single(mask_onshoreA), Rgrid) title('OnshoreA mask: direct access to electricity, not too crowded and land type OK')
    newmap meshm(single(mask_onshoreB), Rgrid) title('OnshoreB mask: longer distance to electricity, not too crowded and land type OK')
end

# ***** Create masks for offshore wind ******

disk2 = fspecial('disk', MIN_PIXELS_TO_SHORE)
shore = single(filter2(disk2, single(regions>0)) > 0)

# Offshore wind elec access same as onshore B (elec access < 9*ELEC_ACCESS_PIXELS km)
# all mask conditions
mask_offshore = gridB & ~shore & (topo > -MAX_DEPTH) & (offshoreregions > 0) & ~protected_area

if SHOW_FIGURES == 2
    newmap meshm(single(shore), Rgrid) plotm(borders.lat, borders.lon,'k') title('No offshore wind within 1 pixel outside the shoreline')
    newmap meshm(single(mask_offshore), Rgrid) title('Offshore mask: longer distance to electricity, not too close to shore and water not too deep')
    plotm(borders.lat, borders.lon,'k')
end

# resize the high resolution masks to fit the ERA data
smallmask_onshoreA = imresize(mask_onshoreA, size(meanwind), 'bilinear')
smallmask_onshoreB = imresize(mask_onshoreB, size(meanwind), 'bilinear')
smallmask_offshore = imresize(mask_offshore, size(meanwind), 'bilinear')

#clear gridA gridB shore goodland gridaccess lat lon topo land





[windonshoreA, windonshoreB, windoffshore] = deal(single(zeros(size(windatlas))))

# make high resolution global maps for the main types of wind resources
windonshoreA(mask_onshoreA>0) = windatlas(mask_onshoreA>0)
windonshoreB(mask_onshoreB>0) = windatlas(mask_onshoreB>0)
windoffshore(mask_offshore>0) = windatlas(mask_offshore>0)

if SHOW_FIGURES >= 1
    newmap meshm(windonshoreA, Rwind) title('OnshoreA wind speeds: direct access to electricity, not too crowded and land type OK')
    plotm(borders.lat, borders.lon,'k')
    newmap meshm(windonshoreB, Rwind) title('OnshoreB wind speeds: longer distance to electricity, not too crowded and land type OK')
    plotm(borders.lat, borders.lon,'k')
    newmap meshm(windoffshore, Rwind) title('Offshore wind speeds: Longer distance to electricity, not too close to shore and water not too deep')
    plotm(borders.lat, borders.lon,'k')
    drawnow
end
if SHOW_FIGURES == 2
    figure histogram(windonshoreA(:),1000) axis([0 12 0 2500]) title('Onshore A annual average wind speeds (spatial variation)')
    figure histogram(windonshoreB(:),1000) axis([0 12 0 2500]) title('Onshore B annual average wind speeds (spatial variation)')
    figure histogram(windoffshore(:),1000) axis([0 12 0 150]) title('Offshore annual average wind speeds (spatial variation)')
    
    smallmask = smallmask_onshoreA & (windERA >= ONSHORECLASSES.min(5)) & (windERA < ONSHORECLASSES.max(5)) & (smallregions == 3)
    smallmask = permute(repmat(smallmask, [1 1 length(windERAtime)]), [3 1 2])
    figure histogram(windERAtime(smallmask), 1000)
    title('Onshore A wind speed distribution of class 5 winds in EUR (mostly temporal variation)')

    smallmask = smallmask_onshoreA & (windERA >= ONSHORECLASSES.min(3)) & (windERA < ONSHORECLASSES.max(3)) & (smallregions == 3)
    smallmask = permute(repmat(smallmask, [1 1 length(windERAtime)]), [3 1 2])
    figure histogram(windERAtime(smallmask), 1000)
    title('Onshore A wind speed distribution of class 3 winds in EUR (mostly temporal variation)')

    smallmask = smallmask_offshore & (windERA >= ONSHORECLASSES.min(5)) & (windERA < ONSHORECLASSES.max(5)) & (smallregions == 3)
    smallmask = permute(repmat(smallmask, [1 1 length(windERAtime)]), [3 1 2])
    figure histogram(windERAtime(smallmask), 1000)
    title('Offshore wind speed distribution of class 5 winds in EUR (mostly temporal variation)')
    
    drawnow
end






onshoreclass = zeros(size(windatlas))
offshoreclass = zeros(size(windatlas))
onshoreclassERA = zeros(size(windatlassmall))
offshoreclassERA = zeros(size(windatlassmall))
for c = 1:length(ONSHORECLASSES.min)
    f = windatlas >= ONSHORECLASSES.min(c) & windatlas < ONSHORECLASSES.max(c)
    onshoreclass(f) = c
    f = windatlassmall >= ONSHORECLASSES.min(c) & windatlassmall < ONSHORECLASSES.max(c)
    onshoreclassERA(f) = c
end
for c = 1:length(OFFSHORECLASSES.min)
    f = windatlas >= OFFSHORECLASSES.min(c) & windatlas < OFFSHORECLASSES.max(c)
    offshoreclass(f) = c
    f = windatlassmall >= OFFSHORECLASSES.min(c) & windatlassmall < OFFSHORECLASSES.max(c)
    offshoreclassERA(f) = c
end

toctic
mdisp('\Calculating GW potential in wind classes...\')

nclasses = numel(ONSHORECLASSES.min)
[capacity_onshoreA, capacity_onshoreB, capacity_offshore] = deal(zeros(numreg,nclasses))

updateprogress = textprogressbar(nlats)
# Working at high resolution here ...
for i = 1:nlats
    for j = 1:nlons
        reg = regions(i,j)
        area = areamatkm(i,j)
        class = onshoreclass(i,j)
        offreg = offshoreregions(i,j)
        offclass = offshoreclass(i,j)
        # none of the ">0" comparisons are needed here, but the loop runs 20x faster with them
        if reg > 0 && class > 0 && mask_onshoreA(i,j) > 0
            capacity_onshoreA(reg,class) = capacity_onshoreA(reg,class) + area
        elseif reg > 0 && class > 0 && mask_onshoreB(i,j) > 0
            capacity_onshoreB(reg,class) = capacity_onshoreB(reg,class) + area
        elseif offreg > 0 && offclass > 0 && mask_offshore(i,j) > 0
            capacity_offshore(offreg,offclass) = capacity_offshore(offreg,offclass) + area
        end
    end
    updateprogress()
end
capacity_onshoreA = 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * capacity_onshoreA
capacity_onshoreB = 1/1000 * ONSHORE_DENSITY * AREA_ONSHORE * capacity_onshoreB
capacity_offshore = 1/1000 * OFFSHORE_DENSITY * AREA_OFFSHORE * capacity_offshore




toctic
mdisp('\Calculating regional and global capacity factors of wind classes...\')

[CFtime_windonshoreA, CFtime_windonshoreB, CFtime_windoffshore] = deal(zeros(yearlength,numreg,nclasses))
[CFtime_windonshoreA_global, CFtime_windonshoreB_global, CFtime_windoffshore_global] = deal(zeros(yearlength,nclasses))
[count_onshoreA, count_onshoreB, count_offshore] = deal(zeros(numreg,nclasses))
[count_onshoreA_global, count_onshoreB_global, count_offshore_global] = deal(zeros(1,nclasses))

updateprogress = textprogressbar(nlatssmall)
# Working at low resolution here ...
for i = 1:nlatssmall
    for j = 1:nlonssmall
        reg = smallregions(i,j)
        class = onshoreclassERA(i,j)
        offreg = smalloffshoreregions(i,j)
        offclass = offshoreclassERA(i,j)
        # none of the ">0" comparisons are needed here, but the loop runs 20x faster with them
        # can't use elseif here, probably some overlap in the masks
        if reg > 0 && class > 0 && smallmask_onshoreA(i,j) > 0
            CFtime_windonshoreA_global(:,class) = CFtime_windonshoreA_global(:,class) + windCF(:,i,j)
            count_onshoreA_global(class) = count_onshoreA_global(class) + 1
            CFtime_windonshoreA(:,reg,class) = CFtime_windonshoreA(:,reg,class) + windCF(:,i,j)
            count_onshoreA(reg,class) = count_onshoreA(reg,class) + 1
        end
        if reg > 0 && class > 0 && smallmask_onshoreB(i,j) > 0
            CFtime_windonshoreB_global(:,class) = CFtime_windonshoreB_global(:,class) + windCF(:,i,j)
            count_onshoreB_global(class) = count_onshoreB_global(class) + 1
            CFtime_windonshoreB(:,reg,class) = CFtime_windonshoreB(:,reg,class) + windCF(:,i,j)
            count_onshoreB(reg,class) = count_onshoreB(reg,class) + 1
        end
        if offreg > 0 && offclass > 0 && smallmask_offshore(i,j) > 0
            CFtime_windoffshore_global(:,offclass) = CFtime_windoffshore_global(:,offclass) + windCF(:,i,j)
            count_offshore_global(offclass) = count_offshore_global(offclass) + 1
            CFtime_windoffshore(:,offreg,offclass) = CFtime_windoffshore(:,offreg,offclass) + windCF(:,i,j)
            count_offshore(offreg,offclass) = count_offshore(offreg,offclass) + 1
        end
    end
    updateprogress()
end
CFtime_windonshoreA = CFtime_windonshoreA ./ permute(count_onshoreA, [3 1 2])
CFtime_windonshoreB = CFtime_windonshoreB ./ permute(count_onshoreB, [3 1 2])
CFtime_windoffshore = CFtime_windoffshore ./ permute(count_offshore, [3 1 2])
CFtime_windonshoreA_global = CFtime_windonshoreA_global ./ count_onshoreA_global
CFtime_windonshoreB_global = CFtime_windonshoreB_global ./ count_onshoreB_global
CFtime_windoffshore_global = CFtime_windoffshore_global ./ count_offshore_global

CF_windonshoreA = squeeze(mean(CFtime_windonshoreA))
CF_windonshoreB = squeeze(mean(CFtime_windonshoreB))
CF_windoffshore = squeeze(mean(CFtime_windoffshore))
CF_windonshoreA_global = mean(CFtime_windonshoreA_global)
CF_windonshoreB_global = mean(CFtime_windonshoreB_global)
CF_windoffshore_global = mean(CFtime_windoffshore_global)

#capacity_onshoreA, capacity_onshoreB, capacity_offshore
#CF_windonshoreA, CF_windonshoreB, CF_windoffshore
#CF_windonshoreA_global, CF_windonshoreB_global, CF_windoffshore_global





toctic
mdisp('\Calculating representative timeseries for wind...\')
CF_wind_time_agg = single(zeros(yearlength,numreg))
count_wind_agg = zeros(1,numreg)

# Calculate representative timeseries for each technology and region.
# In every region, choose all sites in onshore classes A3-A5 # B5 and offshore class 5,
# and average power output(t) over the entire region.
updateprogress = textprogressbar(nlatssmall)
# Working at low resolution here ...
for i = 1:nlatssmall
    for j = 1:nlonssmall
        reg = smallregions(i,j)
        class = onshoreclassERA(i,j)
        offreg = smalloffshoreregions(i,j)
        offclass = offshoreclassERA(i,j)

        if reg > 0 && smallmask_onshoreA(i,j) > 0 && class >= 3 || ...
            reg > 0 && smallmask_onshoreB(i,j) > 0 && class >= 5
                CF_wind_time_agg(:,reg) = CF_wind_time_agg(:,reg) + windCF(:,i,j)
                count_wind_agg(reg) = count_wind_agg(reg) + 1
        # add offshore class 4 for PAS since everything else is empty
        elseif smallmask_offshore(i,j) > 0 && (offreg > 0 && offclass >= 5 || offreg == 9 && offclass >= 4)
            CF_wind_time_agg(:,offreg) = CF_wind_time_agg(:,offreg) + windCF(:,i,j)
            count_wind_agg(offreg) = count_wind_agg(offreg) + 1
        end        
    end
    updateprogress()
end
CF_wind_time_agg = CF_wind_time_agg ./ count_wind_agg

toc

#CF_wind_time_agg(1:20,:)

if SHOW_FIGURES == 2
    figure plot(CF_wind_time_agg)
end

firsttime1 = ncread(['D:/datasets/era5/era5wind' quartername(ERA_YEAR,1) '.nc'], 'time', 1, 1)
firsttime2 = ncread(['D:/datasets/era5/era5solar' quartername(ERA_YEAR,1) '.nc'], 'time', 1, 1)
mdisp('\\Wind and solar ERA5 data not synced!!!!! Time diff 7 hours: ', [firsttime1 firsttime2])

if strcmp(GISREGION, 'global')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_global']
elseif strcmp(GISREGION, 'europe8')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_europe8']
elseif strcmp(GISREGION, 'europe10')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_europe10']
elseif strcmp(GISREGION, 'europe12')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_europe12']
elseif strcmp(GISREGION, 'europe15')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_europe15']
elseif strcmp(GISREGION, 'europe17')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_europe17']
elseif strcmp(GISREGION, 'china6')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_china6']
elseif strcmp(GISREGION, 'eurochine14')
    #filename = ['GISdata_wind' num2str(ERA_YEAR) '_eurochine14']
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_1000km_eurochine14']
elseif strcmp(GISREGION, 'eurasia38')
    #filename = ['GISdata_wind' num2str(ERA_YEAR) '_eurasia38']
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_1000km_eurasia38']
elseif strcmp(GISREGION, 'eurasia21')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_eurasia21']
    #filename = ['GISdata_wind' num2str(ERA_YEAR) '_1000km_eurasia21']
elseif strcmp(GISREGION, 'mena')
    filename = ['GISdata_wind' num2str(ERA_YEAR) '_mena']
end
filename = [filename '_50classes']
save(filename, 'CFtime_windonshoreA', 'CFtime_windonshoreB', 'CFtime_windoffshore', ...
               'capacity_onshoreA', 'capacity_onshoreB', 'capacity_offshore')


#{
# Export to Excel.
mdisp('\Writing Excel spreadsheet...\')
tic
file = ['GISoutput' num2str(ERA_YEAR) '.xlsx']
classnames = {'c1' 'c2' 'c3' 'c4' 'c5'}'

xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [capacity_onshoreA' sum(capacity_onshoreA)'], 'onshoreA capacity (GW)', 'Wind classes', 'A1')
xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [capacity_onshoreB' sum(capacity_onshoreB)'], 'onshoreB capacity (GW)', 'Wind classes', 'O1')
xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [capacity_offshore' sum(capacity_offshore)'], 'offshore capacity (GW)', 'Wind classes', 'A18')

xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [CF_windonshoreA' CF_windonshoreA_global'], 'onshoreA capacity factors (#)', 'Wind classes', 'A9')
xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [CF_windonshoreB' CF_windonshoreB_global'], 'onshoreB capacity factors (#)', 'Wind classes', 'O9')
xlswritetable(file, classnames, {regionlist{:} 'GLOBAL'}, [CF_windoffshore' CF_windoffshore_global'], 'offshore capacity factors (#)', 'Wind classes', 'A26')

time = (1:yearlength)'
xlswrite(file, mean(CF_wind_time_agg), 'Time series', 'B1')
xlswritetable(file, time, regionlist, CF_wind_time_agg, 'Wind capacity factors', 'Time series', 'A3')

parameters = {  'ONSHORE_DENSITY' 'OFFSHORE_DENSITY' 'AREA_ONSHORE' 'AREA_OFFSHORE' 'ELEC_ACCESS_PIXELS' 'PERSONS_PER_KM2' ...
                'MAX_DEPTH' 'MIN_PIXELS_TO_SHORE' 'EXCLUDE_LANDTYPES' 'PROTECTED_CODES' 'ERA_YEAR' 'RESCALE_ERA_TO_WIND_ATLAS' ...
                'ONSHORECLASSES.min' 'ONSHORECLASSES.max' 'OFFSHORECLASSES.min' 'OFFSHORECLASSES.max' }'

parametervalues = {ONSHORE_DENSITY OFFSHORE_DENSITY AREA_ONSHORE AREA_OFFSHORE ELEC_ACCESS_PIXELS PERSONS_PER_KM2 ...
                MAX_DEPTH MIN_PIXELS_TO_SHORE mat2str(EXCLUDE_LANDTYPES) mat2str(PROTECTED_CODES) ERA_YEAR RESCALE_ERA_TO_WIND_ATLAS ...
                mat2str(ONSHORECLASSES.min) mat2str(ONSHORECLASSES.max) mat2str(OFFSHORECLASSES.min) mat2str(OFFSHORECLASSES.max) }'

xlswritetable(file, parameters, '', parametervalues, 'Wind parameters', 'Parameters', 'A1')
toc
#}


=#






# @btime posting2cell($meanwind, $lat0, $lon0, $latrangesmall, $lonrangesmall);
# meanwind2 = posting2cell(meanwind, lat0, lon0, latrangesmall, lonrangesmall)
# convertallpostings!(windCF, lat0, lon0, latrangesmall, lonrangesmall)

function convertallpostings!(windCF, lat0, lon0, lat, lon)
    yearlength = size(windCF, 1)
    updateprogress = Progress(yearlength, 1)
    for i=1:yearlength
        windCF[i,:,:] = posting2cell(windCF[i,:,:], lat0, lon0, lat, lon)
        next!(updateprogress)
    end
end

function posting2cell(data, lat0, lon0, lat, lon)
    # I thought these should both be reversed to do the shift in the right direction,
    # but I get the same result as the Matlab code without the reverses. Also, the
    # spatial correlation with the wind atlas seems better without reverses.
    interp = interpolate((lat0, lon0), data, Gridded(Linear()))       # reverse(data, dims=1)
    extrap = extrapolate(interp, Line())
    testxy(x,y) = y >= lat0[1] && y <= lat0[end] && x >= lon0[1] && x <= lon0[end]
    return [testxy(x,y) ? interp(y,x) : extrap(y,x) for y in lat, x in lon]     # reverse(lat)
end

function posting2cell_slower1(data, lat0, lon0, lat, lon)
    interp = interpolate((lat0, lon0), data, Gridded(Linear()))       # reverse(data, dims=1)
    extrap = extrapolate(interp, Line())
    return [extrap(y,x) for y in lat, x in lon]                                 # reverse(lat)
end

function posting2cell_slower2(data, lat0, lon0, lat, lon)
    interp = LinearInterpolation((lat0, lon0), data, extrapolation_bc = Line())       # reverse(data, dims=1)
    return [interp(y,x) for y in lat, x in lon]                                 # reverse(lat)
end