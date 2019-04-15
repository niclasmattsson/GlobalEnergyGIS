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
        dataset_diffuse = d_create(group, "diffuse", datatype(Float32), dataspace(gridsize...,hours), "chunk", (gridsize...,5), "blosc", 3)
        dataset_direct = d_create(group, "direct", datatype(Float32), dataspace(gridsize...,hours), "chunk", (gridsize...,5), "blosc", 3)
        dataset_meandiffuse = d_create(group, "meandiffuse", datatype(Float32), dataspace(gridsize...), "chunk", gridsize, "blosc", 3)
        dataset_meandirect = d_create(group, "meandirect", datatype(Float32), dataspace(gridsize...), "chunk", gridsize, "blosc", 3)

        totaldiffuse = zeros(gridsize)
        totaldirect = zeros(gridsize)
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
            # ERA5 radiations are in J/m2/period, so for hourly data divide by 3600 to get W/m2
            diffuse = replace(ncdataset["ssrd"][:,:,:], missing => 0.0) .* (land .> 0) ./ 3600
            direct = replace(ncdataset["fdir"][:,:,:], missing => 0.0) .* (land .> 0) ./ 3600

            totaldiffuse = totaldiffuse + sum(diffuse, dims=3)
            totaldirect = totaldirect + sum(direct, dims=3)
            len = size(diffuse,3)
            println("Writing to $filename...")
            dataset_diffuse[:,:,hour:hour+len-1] = diffuse
            dataset_direct[:,:,hour:hour+len-1] = direct
            hour += len
        end
        println("\nWriting mean solar variables to $filename...")
        dataset_meandiffuse[:,:] = totaldiffuse/hours
        dataset_meandirect[:,:] = totaldirect/hours
    end
    nothing
end