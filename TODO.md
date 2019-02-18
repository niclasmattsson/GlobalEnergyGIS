# To do list for porting GIS code from Matlab to Julia

## Datasets

Code to read file formats, rasterization and other preprocessing.

* Regions (GADM 3.6)
* Solar & wind time series (ERA5)
* High resolution annual wind speeds (Global wind atlas)
* Land cover (MODIS)
* Population (GHS-POP)
    - GPWv4: NASA-SEDAC Gridded Population of the World, census data for ~2010 extrapolated to 2020, 30 arc seconds = 1 km
    - Global Human Settlement Population (GHS-POP): uses remote sensing dataset of settlements to downscale GPWv4 to 250m resolution (alternatively 1 km)
    - O'Neill (2016): GPWv3 (not v4, census from ~2000) projected to 2100 consistent with SSPs, "gravity model" for population movement, 7.5 arc minutes = 14 km, "zero artifacts" and national differences clearly visible on map
    - Gao (2017): O'Neill data downscaled to 1 km resolution
    - Global Rural-Urban Mapping Project (GRUMP)

## Sooner

## Later

* Mapping
    - projections
    - basic heatmap plot (projected)
    - lines, points, text (all projected)


## Pandoc table formats

  Right     Left     Center     Default
-------     ------ ----------   -------
     12     12        12            12
    123     123       123          123
      1     1          1             1

Table:  Demonstration of simple table syntax.

-------     ------ ----------   -------
     12     12        12             12
    123     123       123           123
      1     1          1              1
-------     ------ ----------   -------

-------------------------------------------------------------
 Centered   Default           Right Left
  Header    Aligned         Aligned Aligned
----------- ------- --------------- -------------------------
   First    row                12.0 Example of a row that
                                    spans multiple lines.

   Second    row                 5.0 Here's another one. Note
                                    the blank line between
                                    rows.
-------------------------------------------------------------

Table: Here's the caption. It, too, may span
multiple lines.

