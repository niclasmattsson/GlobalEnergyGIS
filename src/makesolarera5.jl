export makesolarera5

# Can optionally zero water cells in the landcover dataset to save a lot of disk space.
function makesolarera5(year=2018, land_cells_only=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    filename = "D:/era5solar$year.h5"
    isfile(filename) && error("File $filename exists in $(pwd()), please delete or rename manually.")

    land = imresize(JLD.load("landcover.jld", "landcover"), gridsize)

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        # create GHI and DHI variables (Global Horizontal Irradiance and Diffuse Horizontal Irradiance)
        dataset_GHI = d_create(group, "GHI", datatype(Float32), dataspace(gridsize...,hours), "chunk", (gridsize...,5), "blosc", 3)
        dataset_DHI = d_create(group, "DHI", datatype(Float32), dataspace(gridsize...,hours), "chunk", (gridsize...,5), "blosc", 3)

        hour = 1
        for month = 1:12, monthhalf = 1:2
            if monthhalf == 1
                firstday, lastday = "01", "15"
            else
                firstday = "16"
                lastday = Dates.daysinmonth(Date("$year-$month"))
            end
            monthstr = lpad(month,2,'0')
            date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
            erafile = "D:/testera5/solar$year-$monthstr$firstday-$monthstr$lastday.nc"

            println("\nReading solar diffuse and direct components from $erafile...")
            ncdataset = Dataset(erafile)
            # ERA5 radiations are in J/m2/period, so for hourly data divide by 3600*1000 to get kW/m2
            GHI = replace(ncdataset["ssrd"][:,:,:], missing => 0.0) .* (land .> 0) ./ (3600*1000)
            DHI = GHI - replace(ncdataset["fdir"][:,:,:], missing => 0.0) .* (land .> 0) ./ (3600*1000)

            len = size(GHI,3)
            println("Writing to $filename...")
            dataset_GHI[:,:,hour:hour+len-1] = GHI
            dataset_DHI[:,:,hour:hour+len-1] = DHI
            hour += len
        end
    end
    nothing
end