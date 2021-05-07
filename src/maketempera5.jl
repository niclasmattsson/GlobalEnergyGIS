export maketempera5

# Can optionally zero cells that are zero in the Global Wind Atlas to save a lot of disk space.
function maketempera5(; year=2018, land_cells_only=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")

    filename = joinpath(datafolder, "era5temp$year.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    land = imresize(JLD.load(joinpath(datafolder, "landcover.jld"), "landcover"), gridsize)

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        dataset_temp = create_dataset(group, "temp", datatype(Float32), dataspace(hours,gridsize...), chunk=(hours,16,16), blosc=3)
        dataset_meantemp = create_dataset(group, "meantemp", datatype(Float32), dataspace(gridsize...), chunk=gridsize, blosc=3)

        totaltemp = zeros(gridsize)
        hour = 1

        count = 0
        for month = 1:12, monthhalf = 1:2
            if monthhalf == 1
                firstday, lastday = "01", "15"
            else
                firstday = "16"
                lastday = Dates.daysinmonth(Date("$year-$month"))
            end
            monthstr = lpad(month,2,'0')
            date = "$year-$monthstr-$firstday/$year-$monthstr-$lastday"
            erafile = joinpath(downloadsfolder, "temp$year-$monthstr$firstday-$monthstr$lastday.nc")

            count += 1
            println("\nFile $count of 24:")
            println("Reading temperatures from $erafile...")
            # Permute dimensions to get hours as dimension 1
            ncdataset = Dataset(erafile)
            temp = permutedims((ncdataset["t2m"][:,:,:] .- 273.15) .* (land .> 0), [3,1,2])

            totaltemp += sumdrop(temp, dims=1)
            len = size(temp,1)
            println("Writing to $filename...")
            dataset_temp[hour:hour+len-1,:,:] = temp
            hour += len
        end
        println("\nWriting meantemp to $filename...")
        dataset_meantemp[:,:] = totaltemp/hours
    end
    nothing
end
