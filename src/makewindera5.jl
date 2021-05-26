export makewindera5, makemonthlywindera5

# Can optionally zero cells that are zero in the Global Wind Atlas to save a lot of disk space.
function makewindera5(; year=2018, windatlas_only=true)
    hours = 24*Dates.daysinyear(year)
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    downloadsfolder = joinpath(datafolder, "downloads")
    
    filename = joinpath(datafolder, "era5wind$year.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    windatlas = reshape(imresize(getwindatlas(), gridsize), (1,gridsize...))

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        dataset_wind = create_dataset(group, "wind", datatype(Float32), dataspace(hours,gridsize...), chunk=(hours,16,16), blosc=3)
        dataset_meanwind = create_dataset(group, "meanwind", datatype(Float32), dataspace(gridsize...), chunk=gridsize, blosc=3)

        totalwind = zeros(gridsize)
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
            erafile = joinpath(downloadsfolder, "wind$year-$monthstr$firstday-$monthstr$lastday.nc")

            count += 1
            println("\nFile $count of 24:")
            println("Reading wind components from $erafile...")
            # Permute dimensions to get hours as dimension 1 (for efficient iteration in GISwind())
            ncdataset = Dataset(erafile)
            u100 = permutedims(ncdataset["u100"][:,:,:], [3,1,2])
            v100 = permutedims(ncdataset["v100"][:,:,:], [3,1,2])

            println("Calculating absolute speed...")
            wind = replace(sqrt.(u100.^2 + v100.^2), missing => 0.0) .* (windatlas .> 0)

            totalwind = totalwind + sumdrop(wind, dims=1)
            len = size(wind,1)
            println("Writing to $filename...")
            dataset_wind[hour:hour+len-1,:,:] = wind
            hour += len
        end
        println("\nWriting meanwind to $filename...")
        dataset_meanwind[:,:] = totalwind/hours
    end
    nothing
end

# Unlike makemonthlysolarera5(), the wind version can't use monthly ERA5 variables
# which take monthly averages of u and v components of wind speed. Doing so would
# greatly underestimate absolute wind speeds in locations where wind direction
# changes frequently. So instead we use ordinary hourly ERA5 data. This assumes
# that every year between 1979-2019 has been downloaded (for wind).
function makemonthlywindera5(; windatlas_only=true)
    years = 1979:2019
    nyears = length(years)
    nmonths = nyears*12
    gridsize = (1280,640)

    datafolder = getconfig("datafolder")
    filename = joinpath(datafolder, "era5monthlywind.h5")
    isfile(filename) && error("File $filename exists in $datafolder, please delete or rename manually.")

    println("Creating HDF5 file:  $filename")
    h5open(filename, "w") do file 
        group = file["/"]
        monthlywind = create_dataset(group, "monthlywind", datatype(Float32), dataspace(nmonths,gridsize...), chunk=(nmonths,16,16), blosc=3)
        annualwind = create_dataset(group, "annualwind", datatype(Float32), dataspace(nyears,gridsize...), chunk=(nyears,16,16), blosc=3)

        for (y, year) in enumerate(years)
            print("$year: ")
            options = WindOptions(merge(windoptions(), Dict(:era_year => year)))
            windatlas, _, meanwind, windspeed = read_wind_datasets(options, 1:36000, 1:18000)
            monthdays = [Dates.daysinmonth(Date("$year-$m")) for m in 1:12]
            lasthour = cumsum(24*monthdays)
            firsthour = [1; lasthour[1:end-1] .+ 1]
            for m = 1:12
                monthlywind[12*(y-1) + m,:,:] =
                    mean(windspeed[firsthour[m]:lasthour[m], :, :], dims=1)
            end
            annualwind[y,:,:] = meanwind
        end
    end
    nothing
end
