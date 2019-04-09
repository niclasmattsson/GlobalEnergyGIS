export makewindera5

# Can optionally filter out cells that are zero in the Global Wind Atlas to save disk space.
function makewindera5(year=2018, filter_windatlas_cells=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    filename = "D:/era5wind$year.h5"
    isfile(filename) && error("File $filename exists in $(pwd()), please delete or rename manually.")

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        dataset_wind = d_create(group, "wind", datatype(Float32), dataspace(gridsize...,hours), "chunk", (gridsize...,5), "blosc", 3)
        dataset_meanwind = d_create(group, "meanwind", datatype(Float32), dataspace(gridsize...), "chunk", gridsize, "blosc", 3)

        totalwind = zeros(gridsize)
        hour = 1

        for month = 1:12, monthhalf = 1:2
            var1, var2 = "100m_u_component_of_wind", "100m_v_component_of_wind"
            if monthhalf == 1
                firstday, lastday = "01", "15"
            else
                firstday = "16"
                lastday = Dates.daysinmonth(Date("$year-$month"))
            end
            monthstr = lpad(month,2,'0')
            date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
            erafile = "D:/testera5/wind$year-$monthstr$firstday-$monthstr$lastday.nc"

            println("\nReading wind components from $erafile...")
            ncdataset = Dataset(erafile)
            u100 = ncdataset["u100"][:,:,:]
            v100 = ncdataset["v100"][:,:,:]

            println("Calculating absolute speed...")
            wind = replace(sqrt.(u100.^2 + v100.^2), missing => 0.0)

            totalwind = totalwind + sum(wind, dims=3)
            len = size(wind,3)
            println("Writing to $filename...")
            dataset_wind[:,:,hour:hour+len-1] = wind
            hour += len
        end
        println("\nWriting meanwind to $filename...")
        dataset_meanwind[:,:] = totalwind/hours
    end
end