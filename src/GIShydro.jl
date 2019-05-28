hydrooptions() = Dict(
    :gisregion => "Eurasia21",                  # "Europe8", "Eurasia38", "Scand3"

    :costclasses_min => [ 0,  50, 100],         # US $/MWh
    :costclasses_max => [50, 100, 999],

    :storageclasses_min => [   0, 1e-6,  12],   # weeks (discharge time)
    :storageclasses_max => [1e-6,   12, 9e9]
)

mutable struct HydroOptions
    gisregion               ::String
    costclasses_min         ::Vector{Float64}
    costclasses_max         ::Vector{Float64}
    storageclasses_min      ::Vector{Float64}
    storageclasses_max      ::Vector{Float64}
end

HydroOptions() = HydroOptions("",[],[],[],[])

function HydroOptions(d::Dict{Symbol,Any})
    options = HydroOptions()
    for (key,val) in d
        setproperty!(options, key, val)
    end
    return options
end

function dms2deg(dms)
    sg = sign(dms)
    dms = abs(dms)
    d = sg*floor(dms/10000)
    dms = dms-10000*d
    m = floor(dms/100)
    return d + m/60 + (dms-100*m)/3600
end

nonmissingtype(x::Type{Union{Missing,T}}) where T = T
nonmissing_eltype(x) = isa(eltype(x), Union) ? nonmissingtype(eltype(x)) : eltype(x)

# Replace all missing values in a DataFrame with simple defaults
function unmissify!(df)
    replacements = Dict(String => "", Int => 0, Float64 => 0.0)
    for col in names(df)
        df[col] = coalesce.(df[col], replacements[nonmissing_eltype(df[col])])
    end
end

# function GIShydro(; optionlist...)

println("\nReading hydro databases...")

# lat,lon,COE,Production_GWh,Lake_surface_m2,Lake_volume_m3,
#   Qm1,Qm2,Qm3,Qm4,Qm5,Qm6,Qm7,Qm8,Qm9,Qm10,Qm11,Qm12,Qm13,ContID,BasinID,SysID,CapCost
# original headers:
# lat,lon,COE ($/kWh),Production (GWh),Lake surface (m2),Lake volume (m3),
#   Qm1,Qm2,Qm3,Qm4,Qm5,Qm6,Qm7,Qm8,Qm9,Qm10,Qm11,Qm12,Qm13,ContID,BasinID,
#   SysID (1=DiversionalCanalPower/2=RiverPower), CapCost ($/kW)
potential = readdlm("C:/Stuff/Datasets/Hydro database (Gernaat) - potential.csv", ',', skipstart=1)

# GrandID,lat,lon,Production_kWh,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13
# original headers:
# GrandID,lat,lon,Production_GWh,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13
existing = readdlm("C:/Stuff/Datasets/Hydro database (Gernaat) - existing (GRanD).csv", ',', skipstart=1)

# country,country_long,name,gppd_idnr,capacity_mw,latitude,longitude,fuel1,fuel2,fuel3,fuel4,commissioning_year,owner,source,url,geolocation_source,year_of_capacity_data,generation_gwh_2013,generation_gwh_2014,generation_gwh_2015,generation_gwh_2016,estimated_generation_gwh
elecplants = CSV.read("C:/Stuff/Datasets/WRI - Global Power Plant Database v1.10/global_power_plant_database.csv", copycols=true)
unmissify!(elecplants)

# clean up wrong coordinates
wrong = findall((elecplants.fuel1 .== "Hydro") .& (
        (elecplants.latitude.>90) .| (elecplants.latitude.<-90) .| (elecplants.longitude.>180) .| (elecplants.longitude.<-180)
    ))
# elecplants.url[wrong]     # 8 results, check urls online
# 1: DDMMSS in lat and lon
# 2-8: lat & lon flipped
w = wrong[1]
elecplants.latitude[w], elecplants.longitude[w] = dms2deg(elecplants.latitude[w]), dms2deg(elecplants.longitude[w])
w = wrong[2:8]
elecplants.latitude[w], elecplants.longitude[w] = elecplants.longitude[w], elecplants.latitude[w] 
wrong = findall((elecplants.fuel1 .== "Hydro") .& (
        (elecplants.latitude.>90) .| (elecplants.latitude.<-90) .| (elecplants.longitude.>180) .| (elecplants.longitude.<-180)
    ))
if !isempty(wrong)
    error("Cleanup of power plant database coordinates failed.")
end

# capacity_mw,latitude,longitude
hydroplants = table2array(elecplants(strcmp(elecplants.fuel1,'Hydro'), 5:7))

dayspermonth = [31 28 31 30 31 30 31 31 30 31 30 31];

% Too little data on Nordic sites (60 degree north limit), so use old monthly inflow instead.
nordic = h5read('C:/Stuff/Julia/Modules/Supergridmodel/inputdata/demand_Europe10.h5','/hydroRoR');
nordic = nordic(:,1);
lasthour = 24 * cumsum(dayspermonth);
firsthour = [1 1+lasthour(1:end-1)];
nordicprofile = zeros(1,12);
for m = 1:12
    nordicprofile(m) = mean(nordic(firsthour(m):lasthour(m)));
end

% Country, Capacity_MW, Pumped_MW, Other_MW, Generation_GWh, Generation_BP_GWh   % [0 = no data]
% original headers:
% Country, Total Hydropower Capacity (MW) in 2015, Pumped Storage Capacity (MW) in 2015, Excluding Pumped Storage (MW) in 2015, Estimated Net Hydropower Generation (GWh) in 2015, Consumption of Hydroelectricity (GWh) in 2015 (BP 2016)
WECcapacity = readtable('C:/Stuff/Datasets/WEC hydro capacity 2015.csv');

% Country, Undeveloped_GWh, Potential_GWh, Utilisation
% original headers:
% Country, Undeveloped (GWh/year), Total Potential (GWh/year), Current Utilisation (%)
WECpotential = readtable('C:/Stuff/Datasets/WEC hydro potentials.csv');

prod = potential(:,4);              % GWh
profile = potential(:, 7:18);       % m3/s
eleccost = potential(:,3);          % $/kWh
capcost = potential(:,23);          % $/kW
type = potential(:,22);
cfriver = mean(profile,2) ./ max(profile,[],2);
reservoirsize = potential(:, 6) / 1e6;        % Mm3
% reservoirsize(reservoirsize==0) = 0.1;

existingprofile = existing(:, 5:16);       % m3/s

crf = CRF(0.1,40);
cf = capcost * crf/8760 ./ eleccost;
capac = prod ./ cf * 1000/8760;    % MW
effic = 0.9*(type==1) + 0.7*(type==2);
sortedflow = sort(profile,2);
designflow = sortedflow(:,9);       % m3/s
waterdensity = 997;                 % kg/m3
grav = 9.81;                        % m/s2

% My estimate of fallheight using eq 2-3 in the Gernaat paper.
% Ask for a new version of the database with this field included.
fallheight = prod*1e9/8760 ./ cf ./ effic ./ (designflow*waterdensity*grav);    % m

% remove outliers
f = find(fallheight > 2000);
fallheight(f) = 2000;
prod(f) = fallheight(f) .* cf(f) .* effic(f) .* designflow(f) * waterdensity*grav*8760/1e9;
capac(f) = prod(f) ./ cf(f) * 1000/8760;
eleccost(f) = eleccost(f) .* potential(f,4) ./ prod(f);
capcost(f) = capcost(f) .* potential(f,4) ./ prod(f);

% J = kg m2/s2
% m * m3/s * kg/m3 * m/s2 = J/s = W
energyprofile = fallheight .* profile * waterdensity*grav .* dayspermonth*24/1e9;  % GWh/month
waterenergy = sum(energyprofile,2);     % GWh/year

% m * Mm3 * kg/m3 * m/s2 = MJ
reservoirenergy = fallheight .* reservoirsize * waterdensity*grav / 3600 / 1000;    % GWh
dischargetime = reservoirenergy ./ capac * 1000;    % h
%figure; hist(log10(reservoirenergy(reservoirenergy>0)),1000)
%hist(log10(dischargetime(reservoirenergy>0)),1000)

% Try calculating a monthly cf based on water inflow and turbine size.
% Maybe not always accurate for large reservoirs since they can store between months,
% but should work for run-of-river. Annual production and CF should still be fine.
% (But why is my production estimate sometimes much higher than Gernaat's?)
monthlyinflow = energyprofile .* effic;
monthlyprod = min(monthlyinflow, capac.*dayspermonth*24/1000);      % GWh/month
prod2 = sum(monthlyprod,2);
cf2 = prod2 ./ (capac*8760/1000);

load regions0_Global_12;
countries = regions;
Rglobal = R;
numcountries = numel(countrynames);

if strcmp(GISREGION, 'global')
    load regions0_Global_12;
    regionlist = countrynames;
    numreg = numel(regionlist);
    bbox = [-90 -180; 90 180];
elseif strcmp(GISREGION, 'europe8')
    load regions_Europe8_GADM12;
    numreg = numel(regionlist)-1;
    regions(regions == numreg+1) = 0;
    bbox = [34 -11; 72 32];
elseif strcmp(GISREGION, 'europe10')
    load regions_Europe10_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [34 -11; 72 32];
elseif strcmp(GISREGION, 'europe12')
    load regions_Europe12_GADM12;
    numreg = numel(regionlist)-1;
    regions(regions == numreg+1) = 0;
    bbox = [34 -11; 72 32];
elseif strcmp(GISREGION, 'eurochine14')
    load regions_EuroChine14_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [-90 -180; 90 180];
elseif strcmp(GISREGION, 'eurasia38')
    load regions_Eurasia38_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [-90 -180; 90 180];
elseif strcmp(GISREGION, 'eurasia21')
    load regions_Eurasia21_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [-90 -180; 90 180];
elseif strcmp(GISREGION, 'china6')
    load regions_China6_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [-90 -180; 90 180];
elseif strcmp(GISREGION, 'mena')
    load regions_MENA_12;
    numreg = numel(regionlist);
    regions(regions == numreg+1) = 0;
    bbox = [10 -15; 45 70];
end

mdisp('/Calculate monthly inflow for existing hydro...');

existingcapac = zeros(numreg,1);             % GW
existinginflow = zeros(numreg,12);           % GWh/month

capacwri = zeros(numcountries,1);    % GW
scalefactor_capacwri = zeros(numcountries,1);
annualcf = zeros(numcountries,1);

for i = 1:length(hydroplants)
    cty = getregion(Rglobal, countries, hydroplants(i,2), hydroplants(i,3));
    if cty == 0
        continue;
    end
    capacwri(cty) = capacwri(cty) + hydroplants(i, 1) / 1e3;    % GW
end

% calculate scale factor for national hydro capacity (WRI capac/WEC capac)
% calculate annual cf for each country from WEC data
for i = 1:size(WECcapacity,1)-1
    cty = find(strcmp(countrynames, WECcapacity{i,1}));
    scalefactor_capacwri(cty) = capacwri(cty) / table2array(WECcapacity(i, 'Capacity_MW')) * 1000;
    annualcf(cty) = table2array(WECcapacity(i, 'Generation_GWh')) / table2array(WECcapacity(i, 'Capacity_MW')) * 1000/8760;
end

% calculate annual generation at each site using capacity and scale factor & CF from above
% then get monthly inflow profile from nearest site in Gernaat/Grand
% finally calculate monthly inflow for existing sites
for i = 1:length(hydroplants)
    cty = getregion(Rglobal, countries, hydroplants(i,2), hydroplants(i,3));
    reg = getregion(R, regions, hydroplants(i,2), hydroplants(i,3));
    if reg == 0
        continue;
    end
    scalefactor = scalefactor_capacwri(cty);
    if scalefactor == 0
        scalefactor = 1;
    end
    capacity = hydroplants(i, 1) / 1e3 / scalefactor;      % GW
    annualgeneration = capacity * annualcf(cty) * 8760;    % GWh/year
    existingcapac(reg) = existingcapac(reg) + capacity;
    
    % approx distance in km to all existing plants in Gernaat/Grand
    dist = 111*sqrt((hydroplants(i,2)-existing(:,2)).^2 + (hydroplants(i,3)-existing(:,3)).^2);
    [d,ndx] = min(dist);
    if strcmp(regionlist(reg),'NOR')
        inflowprofile = nordicprofile;
    else
        inflowprofile = existingprofile(ndx,:);
    end
    monthlyprofile = inflowprofile .* dayspermonth;
    monthlyprofile = monthlyprofile / sum(monthlyprofile);
    
    monthlygeneration = annualgeneration * monthlyprofile;          % GWh/month
    existinginflow(reg,:) = existinginflow(reg,:) + monthlygeneration;  
end

existinginflowcf = existinginflow ./ dayspermonth/24 ./ existingcapac;



mdisp('/Calculate monthly inflow for potential hydro...');

ncostclasses = length(COSTCLASSES.min);
nstorageclasses = length(STORAGECLASSES.min);

potentialcapac = zeros(numreg,ncostclasses,nstorageclasses);        % GW
potentialinflow = zeros(numreg,ncostclasses,nstorageclasses,12);    % GWh/month

potentialmeancost = zeros(numreg,ncostclasses,nstorageclasses);             % $/kWh
potentialmeandischargetime = zeros(numreg,ncostclasses,nstorageclasses);    % hours
nobservations = zeros(numreg,ncostclasses,nstorageclasses);

regfilter = zeros(length(potential),1);

for i = 1:length(potential)
    reg = getregion(R, regions, potential(i,1), potential(i,2));
    if reg == 0
        continue;
    end
    regfilter(i) = reg;

    weeks = dischargetime(i) / 168;
    storageclass = find(weeks >= STORAGECLASSES.min & weeks < STORAGECLASSES.max);
    cost = eleccost(i)*1000;     % US $/MWh
    costclass = find(cost >= COSTCLASSES.min & cost < COSTCLASSES.max);
    
    if strcmp(regionlist(reg),'NOR')
        inflowprofile = nordicprofile * sum(monthlyinflow(i,:))/sum(nordicprofile);
    else
        inflowprofile = monthlyinflow(i,:);
    end
    
    potentialcapac(reg,costclass,storageclass) = potentialcapac(reg,costclass,storageclass) + capac(i)/1000;
    potentialinflow(reg,costclass,storageclass,:) = potentialinflow(reg,costclass,storageclass,:) + reshape(inflowprofile, [1 1 1 12]);

    potentialmeancost(reg,costclass,storageclass) = potentialmeancost(reg,costclass,storageclass) + eleccost(i);
    potentialmeandischargetime(reg,costclass,storageclass) = potentialmeandischargetime(reg,costclass,storageclass) + min(8760, dischargetime(i));
    nobservations(reg,costclass,storageclass) = nobservations(reg,costclass,storageclass) + 1;
end

potentialinflowcf = potentialinflow ./ potentialcapac ./ reshape(dayspermonth, [1 1 1 12])/24;
potentialmeancost = potentialmeancost ./ nobservations;
potentialmeandischargetime = potentialmeandischargetime ./ nobservations;

mdisp('/Saving...');

savename = "GISdata_hydro_" + GISREGION + ".mat";
save(savename, 'existingcapac', 'existinginflowcf', 'potentialcapac', 'potentialinflowcf', 'potentialmeancost', 'potentialmeandischargetime');

mdisp('//Done!/');
return


%hist(cf2(regfilter>0),1000)
%pcf2 = sum(potentialinflow,4) ./ potentialcapac / 8760;
%hist(pcf2(:), 20)





function deg = dms2deg(dms)
    sg = sign(dms);
    dms = abs(dms);
    d = sg*floor(dms/10000);
    dms = dms-10000*d;
    m = floor(dms/100);
    deg = d + m/60 + (dms-100*m)/3600;
end

function crf = CRF(r,T)
    crf = r / (1 - 1/(1+r)^T);
end

function reg = getregion(R, regmat, lat, lon)
    % contains(R,lat,lon) would be nicer but doesn't work on the boundaries
    latlim = lat > R.LatitudeLimits(1) && lat <= R.LatitudeLimits(2);
    lonlim = lon >= R.LongitudeLimits(1) && lon < R.LongitudeLimits(2);
    if latlim && lonlim
        [r,c] = latlon2pix(R,lat,lon);
        reg = regmat(round(r), round(c));
    else
        reg = 0;
    end
end

function url = mapurl(latlon1,latlon2)
    url = "https://www.google.com/maps/dir/?api=1&origin=" + latlon1(1) + "," + latlon1(2);
    url = url + "&destination=" + latlon2(1) + "," + latlon2(2) + "&travelmode=walking";
end

    