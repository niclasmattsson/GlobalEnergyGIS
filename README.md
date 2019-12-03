# GlobalEnergyGIS.jl

Automatic generation of renewable energy input data for energy models in arbitrary world regions using public datasets. Written in Julia. 

## Work in progress

This README is a work in progress and currently only describes the basics of using the package. There are many options that are not yet documented.

## Installation

If you are on Windows, first **restart Julia as administrator** (right click on the Julia shortcut in the Start Menu and select "[More\\]Run as Administrator"). Administrator privileges are not required later to use the package, just for the installation step.

Type `]` to enter Julia's package mode, then:

```
(v1.1) pkg> add https://github.com/niclasmattsson/GlobalEnergyGIS
``` 

Grab some coffee, because installing and compiling dependencies can take quite some time to run. If you don't yet have a Copernicus account, you can create one while you wait for the compilation to complete.

## Copernicus account

GlobalEnergyGIS is based on several public datasets, most importantly hourly global solar insolation and wind speeds from the [ERA5 reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) by the ECMWF (European Centre for Medium-Range Weather Forecasts). To download ERA5 data you need to [create a free account at the Copernicus Data Service (CDS)](https://cds.climate.copernicus.eu/user/register).

Login to your account and go to your profile page (click your name in the upper right), then scroll down to the API key section. Make a note of your user ID (UID) and the long API key string, or just keep the page open so you can copy/paste them in the next step.

## List of datasets and terms of use

The GlobalEnergyGIS package makes use of the following datasets. By using this package, you agree to abide by their terms of use. Please click the links to open the terms of use in your browser.

* ECMWF ERA5 reanalysis and Copernicus download service: https://apps.ecmwf.int/datasets/licences/copernicus/
* Global Wind Atlas version 1: https://globalwindatlas.info/about/TermsOfUse
* World Database of Protected Areas: https://www.protectedplanet.net/c/terms-and-conditions
* GADM (Global Administrative Areas): https://gadm.org/license.html
* Eurostat NUTS (Administrative areas in Europe): https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units
* USGS MODIS 500m Land Cover: https://www.usgs.gov/centers/eros/data-citation
* ETOPO1 Topography: https://www.ngdc.noaa.gov/mgg/global/dem_faq.html#sec-2
* Population scenarios downscaled to 1km resolution: http://www.cgd.ucar.edu/iam/modeling/spatial-population-scenarios.html
* Global population & GDP. Original data: http://www.cger.nies.go.jp/gcp/population-and-gdp.html. Raster converted: https://github.com/Nowosad/global_population_and_gdp.

## Configuration

Run `saveconfig(folder_path, UID, "your API key")` and substitute your own Copernicus data to create a small configuration file (.cdsapirc) in your home directory.

```
julia> using GlobalEnergyGIS

julia> saveconfig("D:/GISdata", 12345, "abcdef123-ab12-cd34-ef56-abcdef123456", agree_terms=true)
```

The first time you run `using GlobalEnergyGIS` there may be another long delay (a minute or two).

## Downloading ERA5 data and disk usage

Run `era5download(data_year)` to begin downloading the ERA5 data for a given year. Since this will require a lot of disk space (see below), you must also provide a path to where to store this data. If you have a fast SSD as a boot drive and a larger HDD, then consider using a folder on the HDD. It may also be a good idea to choose a storage location which doesn't get automatically backuped. Use forward slash `/` or double backslashes `\\` as folder delimiters in the path string. 

```
julia> era5download(2018, "D:/GISdata")
```

One year of global solar and wind data consists of roughly 53 GB of raw data. It will likely take several hours or possibly a day or more to download, depending on your internet connection and whether or not there is congestion in the CDS download queue. A total of 48 files will be downloaded, so remember that the first progress bar you see will be followed by 47 more... :) 

## Conversion and recompression

Next we need to convert the data to a more suitable file format (HDF5). Additionally, to save some disk space in the long run, the raw data will be reduced (by default we discard wind direction, far offshore wind speeds and solar insolation over oceans) and recompressed. This will reduce disk usage to about 22 GB per year of ERA5 data (15.5 GB for solar and 6.5 GB for wind). When this step has finished, the downloaded data will be deleted.

Temporary disk requirements: 

## Downloading other datasets



## Output

