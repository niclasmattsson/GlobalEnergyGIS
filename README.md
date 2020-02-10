# GlobalEnergyGIS.jl

Automatic generation of renewable energy input data for energy models in arbitrary world regions using public datasets. Written in Julia.

<!--
## Mission statement

Supergridmodel. Link paper.
-->

## Disk space requirements

This package uses several large datasets and requires a lot of disk space: roughly 10 GB + 21 GB/year of reanalysis data stored. Also, about 50 GB of **additional** disk space will be required temporarily. Please ensure that you have enough space available (perhaps on a secondary hard drive) before proceeding with the data download. You also need a minimum of 8 GB of RAM memory.

## Installation

If you are on Windows, first **restart Julia as administrator** (right click on the Julia shortcut in the Start Menu and select "Run as Administrator"). Administrator privileges are not required later to use the package, just for the installation step.

Type `]` to enter Julia's package mode, then:

```
(v1.1) pkg> add https://github.com/niclasmattsson/GlobalEnergyGIS
``` 

Grab some coffee, because installing and compiling dependencies can take quite some time to run. If you don't yet have a Copernicus account, you can create one while you wait for the compilation to complete.

## List of datasets and terms of use

The GlobalEnergyGIS package makes use of the following datasets. By using this package, you agree to abide by their terms of use. Please click the links to open the terms of use in your browser.

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

GlobalEnergyGIS is based on several public datasets, most notably the [ERA5 reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) by the ECMWF (European Centre for Medium-Range Weather Forecasts). The reanalysis data is used to calculate hourly solar insolation and wind speeds for any location on earth. To download ERA5 data you need to [create a free account at the Copernicus Data Service (CDS)](https://cds.climate.copernicus.eu/user/register).

### 2. Create config files

Now we will create two config files that will keep track of your preferred data storage path and your Copernicus ID.

Since the datasets will require a lot of disk space (see hardware requirements above), the package allows you to store all data on a secondary disk drive if you have one available. For example, if you have a fast SSD as a boot drive and a larger HDD, then consider storing the data on the HDD. It may also be a good idea to choose a storage location which doesn't get automatically backuped if your backup storage or bandwidth is limited.

[Login to your Copernicus account](https://cds.climate.copernicus.eu/user/login?destination=%2F%23!%2Fhome) and go to your profile page (click your name in the upper right), then scroll down to the API key section. Here you will find your user ID (UID) and the long API key string that needs to be copy/pasted to the command below.

Run `saveconfig(folder_path, Copernicus_UID, "your API key", agree_terms=true)` and substitute your own Copernicus data to create a small configuration file (.cdsapirc) in your home directory. Use forward slash `/` or double backslashes `\\` as folder delimiters in the path string. For example:

```
julia> using GlobalEnergyGIS

julia> saveconfig("D:/GISdata", 12345, "abcdef123-ab12-cd34-ef56-abcdef123456", agree_terms=true)
```

The argument `agree_terms=true` is required to continue. By including it you agree to the terms of use of all datasets listed above. The first time you run `using GlobalEnergyGIS` there will be another delay (a minute or two) while Julia precompiles the dependencies.

### 3. Download auxiliary datasets

Now we will download the auxiliary public datasets listed above. The following command will download them all to the data folder you supplied when you created the configuration file.

```
julia> download_datasets()
```

A total of 16 files will be downloaded and unpacked. This will probably take 15 minutes or so depending on your internet connection.

### 4. Rasterization

Several of the datasets are vector data (shapefiles). To speed up and simplify the subsequent processing, we will rasterize them all to a common resolution (0.01 degrees, or roughly 1 km at the equator).

```
julia> rasterize_datasets(cleanup=:all)
```

This command will automatically delete the original datasets to save disk space. Use the argument `cleanup=:limited` to keep the original files, or `cleanup=:none` to also keep intermediate raster files.

### 5. Download ERA5 data and disk usage

Run `era5download(data_year)` to begin downloading the ERA5 solar and wind data for a given year.

```
julia> era5download(2018)
```

One year of global solar and wind data consists of roughly 53 GB of raw data. It will likely take several hours or possibly a day or more to download, depending on your internet connection and whether or not there is congestion in the CDS download queue. A total of 48 files will be downloaded, so please note that the first progress bar you see will be followed by 47 more... :) 

### 6. Conversion and recompression of ERA5 data

Next we need to convert the ERA5 data to a more suitable file format (HDF5). Additionally, to save some disk space in the long run, the raw data will be reduced (by default we discard wind direction, far offshore wind speeds and solar insolation over oceans) and recompressed. This will reduce disk usage to about 22 GB per year of ERA5 data (15.5 GB for solar and 6.5 GB for wind). 

```
julia> makewindera5(year=2018)
```

Now you can optionally delete the original data. One of the upstream packages sometimes keeps a lock on the files, so you may need to restart Julia to run this command.

```
julia> clearvars_era5("wind", year=2018)
```

Now convert the solar data, and optionally delete the original data.

```
julia> makesolarera5(year=2018)
```

```
julia> clearvars_era5("solar", year=2018)
```

## Usage (creating renewable energy input data for arbitrary model regions)

There are three main steps:

1. Create economic background scenario (run only once per combination of SSP scenario and year)
2. Create region files (run only once per set of model regions)
3. Calculate potential capacities and hourly capacity factors for solar-, wind- and hydropower.

The output of steps 2 and 3 are directly read by [Supergridmodel](https://github.com/niclasmattsson/Supergrid) (a companion energy system model designed to work with this GIS package).

### Create economic background scenario

Run `create_scenario_datasets(SSPscenario, target_year)`, where `SSPscenario` is one of `"SSP1"`, `"SSP2"` or `"SSP3"` and target_year is one of 2020, 2030, ..., 2100.

```
julia> create_scenario_datasets("SSP2", 2050)
```

### Create region files

Run `saveregions(regionset_name, region_definitions)`, where `regionset_name` is a string and `region_definitions` is a matrix as specified below. Then run `makedistances(regionset_name)` to determine which regions are connected onshore and offshore and calculate distances between population-weighted region centers.

```
julia> saveregions("Europe8", europe8)

julia> makedistances("Europe8")

```

Here `europe8` is a region matrix defined in `regiondefinitions.jl`.

### Region definition matrix syntax

Regions are specified using an n x 2 matrix, where the first column contains names of model regions, and the second column contains information on which subregions are included in each main region using GADM or NUTS subregion names. This is facilitated using the `GADM()` and `NUTS()` helper functions. Here's a simple example of making a 4-region definition matrix of Scandiavian countries using both administrative border datasets:

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

Subregions of the same level can be concatenated by listing multiple arguments to the `GADM()` or `NUTS()` call. For example, mainland Portugal can be defined by `NUTS("PT11","PT15","PT16","PT17","PT18")`. This will exclude islands that are otherwise included in `NUTS("PT")`.

Subregion names (codes) are unique in NUTS but not always in GADM. For this reason, GADM subregions must indicate their parent regions. Here are two examples:

```
gadm_subregions = [
    "India_E"   GADM(["India"], "Odisha", "Jharkhand", "West Bengal", "Bihar")
    "Öland"     GADM(["Sweden","Kalmar"], "Borgholm", "Mörbylånga")
]
```

If the first argument to a `GADM()` call is a vector, then the remaining arguments are subregions to the last vector element. Here India_E (Eastern India) is defined using four level-1 subregions of the "India" level-0 region. Next we define the Swedish island of Öland using its two municipalities that are level-2 subregions of Kalmar, which is itself a level-1 subregion of Sweden. If the first argument to a `GADM()` call is not a vector, then all arguments are assumed to be level-0 regions.

This vector syntax is not used in `NUTS()` calls since all subregion code names are unique.

To concatenate different levels of subregions or to mix NUTS and GADM calls in the same region, use a tuple of GADM and NUTS calls by enclosing them in parentheses:

```
concatenation_examples = [
    "China_SC"  (GADM(["China"], "Henan","Hubei","Hunan","Guangdong","Guangxi","Hainan"), GADM("Hong Kong","Macao"))
    "France"    (NUTS("FR"), GADM("Monaco"))
]
```

Here China_SC (southcentral) is defined using six level-1 subregions of China, in addition to Hong Kong and Macao which are level-0 regions in GADM. Next, we define a France NUTS region that includes Monaco (which is not included in NUTS) by concatenating the GADM definition of Monaco.

For more syntax examples, see `regiondefinitions.jl` (in the GlobalEnergyGIS /src folder). [Maps of NUTS regions can be found here.](https://ec.europa.eu/eurostat/web/nuts/nuts-maps).

### The actual GIS analysis

Finally we have everything we need to actually calculate potential capacities and hourly capacity factors for solar-, wind- and hydropower. This is the basic syntax, which assumes default values for all unlisted GIS parameters.

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

The parameters that can be changed are listed in the section "GIS options" below.

## Synthetic demand demo

The synthetic demand code still requires a few more days of development before it can generate demand for arbitrary world regions. All components are in place except hourly temperature data for model regions. When completed, the code will sample hourly temperatures from the three most highly populated pixels of each model region. This will require downloading a temperature dataset from the ERA5 reanalysis in addition to the solar and wind data. 

This will be completed in the week of January 20-24. For now, we have hardcoded temperature data for two sample regions, MENA (Middle East and North Africa) and EU10. Below are instructions for generating synthetic demand for these regions.

### Installation of Python and geopandas

Download and install Anaconda (Python 3.7). Accept the default options in the installer.
https://www.anaconda.com/distribution

Open the Anaconda prompt (from the start menu on a Windows system) and install geopandas. Accept default options if prompted.

```
(base) C:\Users\niclas>conda install geopandas
```

### Creating input data for synthetic demand

Skip this step for now. The requisite input data for MENA and EU10 are already included.

In Julia:

```
julia> using GlobalEnergyGIS

julia> regions, regionlist, demand, pop, popdens, gdp = makesyntheticdemandinput(gisregion="Europe8");
```

### Running the synthetic demand module

In the Anaconda prompt, change to the synthetic demand folder and start the script, using the region name as a command line argument (MENA or EU10).

```
(base) C:\Users\niclas> cd D:\GISdata\syntheticdemand

(base) D:\GISdata\syntheticdemand>python synthetic_demand_demo.py MENA
```

The machine learning code runs for about 5-10 minutes. The resulting output files will be placed in `D:/GISdata/syntheticdemand/output/`. 

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
