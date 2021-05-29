# GlobalEnergyGIS.jl

Automatic generation of renewable energy input data for energy models in arbitrary world regions using
public datasets. Written in Julia.

## Paper

The work here has been described in a
[scientific paper](paper/Mattsson%20et%20al.%202019%20-%20An%20autopilot%20for%20energy%20models.pdf)
submitted to Energy Strategy Reviews along with its
[supplementary material](paper/Mattsson%20et%20al.%202019%20-%20Supplementary%20-%20An%20autopilot%20for%20energy%20models.pdf).

## Disk space requirements

This package uses several large datasets and requires a lot of disk space: roughly 10 GB + 29 GB/year of
reanalysis data stored. Also, at least 50 GB of **additional** disk space will be required temporarily. Please
ensure that you have enough space available (perhaps on a secondary hard drive) before proceeding with the
data download. You also need a minimum of 8 GB of RAM memory.

## Installation

Make sure you are on Julia v1.6 or higher. Type `]` to enter Julia's package mode, then:

```
(@v1.6) pkg> add https://github.com/niclasmattsson/GlobalEnergyGIS
``` 

Grab some coffee, because installing and compiling dependencies can take quite some time to run. If you don't
yet have a Copernicus account, you can create one while you wait for the compilation to complete.

## List of datasets and terms of use

The GlobalEnergyGIS package makes use of the following datasets. By using this package, you agree to abide by
their terms of use. Please click the links to open the terms of use in your browser.

* ECMWF ERA5 reanalysis and Copernicus download service: https://apps.ecmwf.int/datasets/licences/copernicus/
* Global Wind Atlas version 1: https://globalwindatlas.info/about/TermsOfUse
* World Database of Protected Areas: https://www.protectedplanet.net/c/terms-and-conditions
* GADM (Global Administrative Areas): https://gadm.org/license.html
* Eurostat NUTS (Administrative areas in Europe): https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units
* USGS MODIS 500m Land Cover: https://www.usgs.gov/centers/eros/data-citation
* ETOPO1 Topography: https://www.ngdc.noaa.gov/mgg/global/dem_faq.html#sec-2
* Population scenarios downscaled to 1 km resolution: http://www.cgd.ucar.edu/iam/modeling/spatial-population-scenarios.html
* Global population & GDP. Original data: http://www.cger.nies.go.jp/gcp/population-and-gdp.html. Raster converted: https://github.com/Nowosad/global_population_and_gdp.

## Setup and data preparation

### 1. Create Copernicus account

GlobalEnergyGIS is based on several public datasets, most notably the
[ERA5 reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) by the ECMWF
(European Centre for Medium-Range Weather Forecasts). The reanalysis data is used to calculate hourly solar
insolation and wind speeds for any location on earth. To download ERA5 data you need to [create a free
account at the Copernicus Data Service (CDS)](https://cds.climate.copernicus.eu/user/register).

### 2. Create config files and agree to dataset terms

Now we will create two config files that will keep track of your preferred data storage path and your
Copernicus ID.

Since the datasets will require a lot of disk space (see hardware requirements above), the package allows
you to store all data on a secondary disk drive if you have one available. For example, if you have a fast
SSD as a boot drive and a larger HDD, then consider storing the data on the HDD. It may also be a good idea
to choose a storage location which doesn't get automatically backuped if your backup storage or bandwidth is
limited.

[Login to your Copernicus account](https://cds.climate.copernicus.eu/user/login?destination=%2F%23!%2Fhome)
and go to your profile page (click your name in the upper right), then scroll down to the API key section.
Here you will find your user ID (UID) and the long API key string that needs to be copy/pasted to the command
below.

Run `saveconfig(folder_path, Copernicus_UID, "your API key", agree_terms=true)` and substitute your own
Copernicus data to create a small configuration file (.cdsapirc) in your home directory. Use forward slash
`/` or double backslashes `\\` as folder delimiters in the path string. For example:

```
julia> using GlobalEnergyGIS

julia> saveconfig("D:/GISdata", 12345, "abcdef123-ab12-cd34-ef56-abcdef123456", agree_terms=true)
```

The argument `agree_terms=true` is required to continue. By including it you agree to the terms of use of all
datasets listed above. The first time you run `using GlobalEnergyGIS` there will be another delay (a minute
or two) while Julia precompiles the dependencies.

To agree to Copernicus terms, you **must** download a small amount of test data once using the web interface
(otherwise you will get errors in step 5 below). Visit [the CDS web interface](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=form),
select one checkbox in each parameter section (e.g. Product type, Variable, Pressure level, etc.),
and finally click the "Agree terms" and "Submit form" buttons as in the screenshot below. If you have any problem with this, follow
[step 2 in these detailed instructions](https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-1-Prerequisites).

![screenshot of CDS terms](https://github.com/niclasmattsson/GlobalEnergyGIS/blob/master/CDSterms.png)

### 3. Download auxiliary datasets

Now we will download the auxiliary public datasets listed above. The following command will download them all
to the data folder you supplied when you created the configuration file.

```
julia> download_datasets()
```

A total of 17 files will be downloaded and unpacked. This will probably take 30 minutes or so depending on
your internet connection.

### 4. Rasterization

Several of the datasets are vector data (shapefiles). To speed up and simplify the subsequent processing, we
will rasterize them all to a common resolution (0.01 degrees, or roughly 1 km at the equator).

```
julia> rasterize_datasets(cleanup=:all)
```

This command will automatically delete the original datasets to save disk space. Use the argument
`cleanup=:limited` to keep the original files, or `cleanup=:none` to also keep intermediate raster files.

### 5. Download and convert ERA5 data

Run `download_and_convert_era5(data_year)` to begin downloading the ERA5 wind, solar and temperature data for a given year. About 26 GB of raw data must be downloaded (split into 24 files of 1.1 GB) *per variable, per year*. So 72 downloads in total. It will likely take several hours to download, depending on your internet connection and whether or not there is congestion in the CDS download queue.

```
julia> download_and_convert_era5(2018)
```

After the raw data for a variable has been downloaded, it is immediately aggregated and converted into more suitable file format (HDF5) and recompressed to save disk space. To minimize disk usage we also throw away data we won't need (by default we discard wind direction, far offshore wind speeds and solar insolation over oceans). This will reduce disk usage to about 31 GB per year of ERA5 data in total (on average 14.7 GB for solar, 8.6 GB for wind and 7.4 GB for temperatures - solar needs roughly twice the space since we need to store both direct and diffuse solar insolation). 

## Usage (creating renewable energy input data for arbitrary model regions)

There are three main steps:

1. Create economic background scenario (run only once per combination of SSP scenario and year)
2. Create region files (run only once per set of model regions)
3. Calculate potential capacities and hourly capacity factors for solar-, wind- and hydropower.

The output of steps 2 and 3 are directly read by [Supergridmodel](https://github.com/niclasmattsson/Supergrid)
(a companion energy system model designed to work with this GIS package).

### Create economic background scenario

Run `create_scenario_datasets(SSPscenario, target_year)`, where `SSPscenario` is one of `"SSP1"`, `"SSP2"`
or `"SSP3"` and target_year is one of 2020, 2030, ..., 2100.

```
julia> create_scenario_datasets("SSP2", 2050)
```

### Create region files

Run `saveregions(regionset_name, region_definitions)`, where `regionset_name` is a string and
`region_definitions` is a matrix as specified below. Then run `makedistances(regionset_name)` to determine
which regions are connected onshore and offshore and calculate distances between population-weighted region
centers.

```
julia> saveregions("Europe13", europe13)

julia> makedistances("Europe13")

```

Here `europe13` is a region matrix defined in `regiondefinitions.jl`, but you can refer to your own region
matrices defined elsewhere. See the next section for syntax examples. To get visual confirmation of the
results, run `createmaps(regionset_name)` to create images of onshore and offshore region territories
(in `/GISdata_folder_path/output`).

```
julia> createmaps("Europe13")
```

### Region definition matrix syntax

Regions are specified using an n x 2 matrix, where the first column contains names of model regions, and the
second column contains information on which subregions are included in each main region using GADM or NUTS
subregion names. This is facilitated using the `GADM()` and `NUTS()` helper functions. Here's a simple
example of making a 4-region definition matrix of Scandiavian countries using both administrative border
datasets:

```
scandinavia4 = [
    "SWE"   GADM("Sweden")
    "NOR"   GADM("Norway")
    "DEN"   GADM("Denmark")
    "FIN"   GADM("Finland", "Åland")    # in GADM, the island of Åland is a separate level-0 region
]

scandinavia4_nuts = [
    "SWE"   NUTS("SE")
    "NOR"   NUTS("NO")
    "DEN"   NUTS("DK")
    "FIN"   NUTS("FI")
]
```

Subregions of the same level can be concatenated by listing multiple arguments to the `GADM()` or `NUTS()`
call. For example, mainland Portugal can be defined by `NUTS("PT11","PT15","PT16","PT17","PT18")`. This will
exclude islands that are otherwise included in `NUTS("PT")`.

Subregion names (codes) are unique in NUTS but not always in GADM. For this reason, GADM subregions must
indicate their parent regions. Here are two examples:

```
gadm_subregions = [
    "India_E"   GADM(["India"], "Odisha", "Jharkhand", "West Bengal", "Bihar")
    "Öland"     GADM(["Sweden","Kalmar"], "Borgholm", "Mörbylånga")
]
```

If the first argument to a `GADM()` call is a vector, then the remaining arguments are subregions to the last
vector element. Here India_E (Eastern India) is defined using four level-1 subregions of the "India" level-0
region. Next we define the Swedish island of Öland using its two municipalities that are level-2 subregions
of Kalmar, which is itself a level-1 subregion of Sweden. If the first argument to a `GADM()` call is not a
vector, then all arguments are assumed to be level-0 regions.

This vector syntax is not used in `NUTS()` calls since all subregion code names are unique.

To concatenate different levels of subregions or to mix NUTS and GADM calls in the same region, use a tuple
of GADM and NUTS calls by enclosing them in parentheses:

```
concatenation_examples = [
    "China_SC"  (GADM(["China"], "Henan","Hubei","Hunan","Guangdong","Guangxi","Hainan"), GADM("Hong Kong","Macao"))
    "France"    (NUTS("FR"), GADM("Monaco"))
]
```

Here China_SC (southcentral) is defined using six level-1 subregions of China, in addition to Hong Kong and
Macao which are level-0 regions in GADM. Next, we define a France NUTS region that includes Monaco (which is
not included in NUTS) by concatenating the GADM definition of Monaco.

For more syntax examples, see `regiondefinitions.jl` (in the GlobalEnergyGIS /src folder).
[Maps of NUTS regions can be found here.](https://ec.europa.eu/eurostat/web/nuts/nuts-maps).

### Subregion helper function

There is a simple helper function `subregion()` to facilitate finding the names of subregions. Note that the
GADM version takes multiple subregion arguments while the NUTS version only takes a single argument (and matches
the beginning of the subregion name). The function will return a vector of subregion names.

```
julia> subregions(GADM)
256-element Array{String,1}:
 "Afghanistan"
 "Akrotiri and Dhekelia"
 "Albania"
 "Algeria"
[...]

julia> subregions(GADM, "France")
13-element Array{String,1}:
 "Auvergne-Rhône-Alpes"
 "Bourgogne-Franche-Comté"
 "Bretagne"
 "Centre-Val de Loire"
 [...]

julia> subregions(GADM, "France", "Bretagne")
4-element Array{String,1}:
 "Côtes-d'Armor"
 "Finistère"
 "Ille-et-Vilaine"
 "Morbihan"

julia> subregions(NUTS)
37-element Array{String,1}:
 "AL"
 "AT"
 "BE"
 "BG"
 [...]

julia> subregions(NUTS, "UK")
179-element Array{String,1}:
 "UKC11"
 "UKC12"
 "UKC13"
 "UKC14"
 [...]

julia> subregions(NUTS, "UKN")
11-element Array{String,1}:
 "UKN06"
 "UKN07"
 "UKN08"
 "UKN09"
 [...]

```

### The actual GIS analysis

Finally we have everything we need to actually calculate potential capacities and hourly capacity factors for
solar-, wind- and hydropower. This is the basic syntax, which assumes default values for all unlisted GIS
parameters.

```
julia> GISsolar(gisregion="Europe13")

julia> GISwind(gisregion="Europe13")

julia> GIShydro(gisregion="Europe13")

```

Here is a call that changes some parameters:

```
julia> GISwind(gisregion="Europe13", scenarioyear="ssp2_2020", era_year=2016, persons_per_km2=100,
				max_depth=60, min_shore_distance=2, area_onshore=0.05, area_offshore=0.20)
```

The parameters that can be changed are listed in the section "GIS options" below. `GISwind()` and
`GISsolar()` also take an optional boolean parameter `plotmasks=true` that will generate .png images
of the dataset masks resulting from the other parameters. This will increase run times by a minute
or two. Images will be placed in `/GISdata_folder_path/output`.

```
julia> GISwind(gisregion="Europe13", ..., plotmasks=true)

julia> GISsolar(gisregion="Europe13", ..., plotmasks=true)
```

## Synthetic electricity demand

The synthetic demand module estimates the profile of hourly electricity demand in each model region given the
total annual demand (determined by current national electricity demand per capita extrapolated using the SSP
background scenario and target year). This is done using machine learning, specifically a method called
gradient boosting tree regression. This is similar to ordinary regression, except that underlying mathematical
relationships between variables are determined automatically using a black box approach.

We train the model based on real electricity demand in 44 countries for the year 2015. Regression variables
include calendar effects (e.g. hour of day and weekday/weekend indicators), temperature variables (e.g. hourly
temperature series in the most populated areas of each model region, or monthly averages as seasonality
indicators) and economic indicators, e.g. local GDP per capita or electricity demand per capita (using the
latter variable is not "cheating", since we are merely interested in predicting hourly profiles of normalized
demand, not the demand level).

### Easy version using our default parameters and regression variables

Assuming you have already downloaded the requisite temperature data for the year you want to study
(see `era5download()` and `maketempera5()`), and created population and GDP datasets for the SSP scenario
(see `create_scenario_datasets()`):

```
julia> predictdemand(gisregion="Europe8", sspscenario="ssp2-26", sspyear=2050, era_year=2018)
```

This will create a matrix (size 8760x`number_of_regions`) with the predicted electricity demand for each
model region and hour of the year. This data is saved in a new JLD file in `/GISdata_folder_path/output`. Here the full SSP scenario variant must be specified including the 2-digit code representing radiative forcing target (e.g. 19, 26, 34, 45).

### Selecting variables to train on

These are the default nine variables (so this will produce the exact same result as the previous example).

```
julia> selectedvars = [:hour, :weekend01, :temp_monthly, :ranked_month, :temp_top3,
                        :temp1_mean, :temp1_qlow, :temp1_qhigh, :demandpercapita]

julia> predictdemand(variables=selectedvars, gisregion="Europe8", sspscenario="ssp2-26", sspyear=2050, era_year=2018)
```

And here is a simpler example using seven variables:

```
julia> selectedvars = [:hour, :weekend01, :ranked_month, :temp_top3, :temp1_qlow, :temp1_qhigh, :gdppercapita]

julia> predictdemand(variables=selectedvars, gisregion="Europe8", sspscenario="ssp2-26", sspyear=2050, era_year=2018)
```

Currently we calculate data for 12 different variables. Any combination of these can be used to train on.
The full list along with brief explanations appears below near the bottom of this README.

### Selecting custom learning parameters

These are the default parameters:

```
julia> predictdemand(variables=defaultvariables, gisregion="Europe8", sspscenario="ssp2-34", sspyear=2050, era_year=2018,
            nrounds=100, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"])         
```

And here we modify them:

```
julia> predictdemand(variables=defaultvariables, gisregion="Europe8", sspscenario="ssp1-45", sspyear=2030, era_year=2018,
            nrounds=40, max_depth=8, eta=0.30, subsample=0.85, metrics=["rmse"])         
```

These parameters are explained briefly below. Additionally, any other
[XGBoost parameters](https://xgboost.readthedocs.io/en/latest/parameter.html) can be specified.

### Cross-validation of the training demand dataset (44 countries)

Cross-validation of the training dataset can help determine which variables to train on and what values of
learning parameters to use. This will predict the demand for all 44 countries in the demand dataset, but the
model built for each country will only use data from the other 43 countries.

Iterations appear significantly slower than `predictdemand()` because it trains 44 models in parallel. 

```
julia> crossvalidate(variables=defaultvariables, nrounds=100, max_depth=7, eta=0.05, subsample=0.75, metrics=["mae"])
```

Note that there is no `gisregion` argument since we are both training and predicting the same 44 country
dataset. The log will show two columns of training errors. The right column `cv-train-mae` will have lower
errors, but this results from evaluating the model on the same data it trained on (i.e "cheating"). The
column `cv-test-mae` on the left is the real (non-cheating) result. Resist the temptation of adapting
parameters to the right column.

## GIS options

### Default GISsolar() options

```
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
    :exclude_landtypes => [0,1,2,3,4,5,8,12],       # exclude water, forests and croplands. See codes in table below.
    :protected_codes => [1,2,3,4,5,8],  # IUCN codes to be excluded as protected areas. See codes in table below.

    :scenarioyear => "ssp2_2050",       # default scenario and year for population and grid access datasets
    :era_year => 2018,                  # which year of the ERA5 time series to use 

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :pvclasses_min => [0.08,0.14,0.18,0.22,0.26],   # lower bound on annual PV capacity factor for class X    [0:0.01:0.49;]
    :pvclasses_max => [0.14,0.18,0.22,0.26,1.00],   # upper bound on annual PV capacity factor for class X    [0.01:0.01:0.50;]
    :cspclasses_min => [0.10,0.18,0.24,0.28,0.32],  # lower bound on annual CSP capacity factor for class X
    :cspclasses_max => [0.18,0.24,0.28,0.32,1.00]  # upper bound on annual CSP capacity factor for class X
)
```

### Default GISwind() options

```
windoptions() = Dict(
    :gisregion => "Europe8",            # "Europe8", "Eurasia38", "Scand3"
    :filenamesuffix => "",              # e.g. "_landx2" to save high land availability data as "GISdata_solar2018_Europe8_landx2.mat" 

    :onshore_density => 5,              # about 30% of existing farms have at least 5 W/m2, will become more common
    :offshore_density => 8,             # varies a lot in existing parks (4-18 W/m2)
                                        # For reference: 10D x 5D spacing of 3 MW turbines (with 1D = 100m) is approximately 6 MW/km2 = 6 W/m2
    :area_onshore => .08,               # area available for onshore wind power after the masks have been applied
    :area_offshore => .33,              # area available for offshore wind power after the masks have been applied

    :distance_elec_access => 300,       # max distance to grid [km] (for wind classes of category B and offshore)
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

    :res => 0.01,                       # resolution of auxiliary datasets [degrees per pixel]
    :erares => 0.28125,                 # resolution of ERA5 datasets [degrees per pixel]

    :onshoreclasses_min => [2,5,6,7,8],     # lower bound on annual onshore wind speeds for class X    [0:0.25:12.25;]
    :onshoreclasses_max => [5,6,7,8,99],    # upper bound on annual onshore wind speeds for class X    [0.25:0.25:12.5;]
    :offshoreclasses_min => [3,6,7,8,9],    # lower bound on annual offshore wind speeds for class X
    :offshoreclasses_max => [6,7,8,9,99]    # upper bound on annual offshore wind speeds for class X
)
```

### Default GIShydro() options

```
hydrooptions() = Dict(
    :gisregion => "Europe8",                    # "Europe8", "Eurasia38", "Scand3"

    :costclasses_min => [ 0,  50, 100],         # US $/MWh
    :costclasses_max => [50, 100, 999],

    :storageclasses_min => [   0, 1e-6,  12],   # weeks (discharge time)
    :storageclasses_max => [1e-6,   12, 9e9]
)
```

### Land types

```
 0      'Water'                       
 1      'Evergreen Needleleaf Forests'
 2      'Evergreen Broadleaf Forests' 
 3      'Deciduous Needleleaf Forests'
 4      'Deciduous Broadleaf Forests' 
 5      'Mixed Forests'               
 6      'Closed Shrublands'           
 7      'Open Shrublands'             
 8      'Woody Savannas'              
 9      'Savannas'                    
10      'Grasslands'                  
11      'Permanent Wetlands'          
12      'Croplands'                   
13      'Urban'                       
14      'Cropland/Natural'            
15      'Snow/Ice'                    
16      'Barren'
```               

### Protected areas (IUCN codes from the WDPA)

```  
1      'Ia'                'Strict Nature Reserve'          
2      'Ib'                'Wilderness Area'                
3      'II'                'National Park'                  
4      'III'               'Natural Monument'               
5      'IV'                'Habitat/Species Management'     
6      'V'                 'Protected Landscape/Seascape'   
7      'VI'                'Managed Resource Protected Area'
8      'Not Reported'      'Not Reported'                   
9      'Not Applicable'    'Not Applicable'                 
10     'Not Assigned'      'Not Assigned'        
```  

<!-- 
## Output file format

To do.


## Include this in README later

* mission statement
* hardware requirements
* add hydro and other data to dataset list
* describe output file formats

 -->

### Synthetic demand: list of training variables

List of variables in the training dataset that can be use for the regression:

```
Calendar variables
   :hour               hour of day
   :month              month of year
   :weekend01          weekend indicator

Hourly temperatures:
   :temp1              temperature in the largest population center of each region
   :temp_top3          average temperature in the three largest population centers of each region

Monthly temperatures (season indicators):
   :temp_monthly       average monthly temperature in the largest population center of each region
   :ranked_month       rank of the average monthly temperature of each month (1-12)

Annual temperature levels and variability:
   :temp1_mean         average annual temperature in the largest population center of each region
   :temp1_qlow         low annual temperature - 5% quantile of hourly temperatures
   :temp1_qhigh        high annual temperature - 95% quantile of hourly temperatures

Economic indicators:
   :demandpercapita    level of annual average electricity demand [MWh/year/capita] in each region
   :gdppercapita       level of annual average GDP per capita [USD(2010)/capita] in each region
```

### Default learning parameters for the synthetic demand regression

Our default values are adapted to our data and differ from the default XGBoost parameters. The parameters
also need to be adapted to each other. For example, a lower `eta` value may require a higher `nrounds` value
to reach full benefit. 

```
nrounds=100      # number of rounds of learning improvements
max_depth=7      # tree depth, i.e. complexity of the underlying black box model. Increasing this may lead to overfitting.
eta=0.05         # learning rate. Higher values will improve faster, but may ultimately lead to a less efficient model.
subsample=0.75   # how much of the training data to use in each iteration. Same tradeoff as 'eta' parameter.
metrics=["mae"]  # "mae" = mean absolute error, "rmse" = root mean square error, or both. Note the brackets (it's a vector).
```

A slightly better but much more computationally demanding set of parameters: `nnrounds=1000`, `max_depth=7`, `eta=0.005`, and `subsample=0.05`. These were the parameters used to produce figures 1 and 2 in the paper above.

For more information on these and other selectable parameters, see https://xgboost.readthedocs.io/en/latest/parameter.html.
