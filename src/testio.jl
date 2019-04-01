export testread, testread2, testread3, testread4, testh5write, testh5read, testread6, testread7, testread5b, testread6b

function testread()
    filename = "D:/datasets/era5/era5windCF2016.h5"
    dataset = "windCF"
    yearlength = 8760 + 24*leapyear(2016)
    n = 100
    updateprogress = Progress(yearlength÷n, 1)
    for t=1:n:yearlength
        windCF = h5read(filename, dataset, (t:min(yearlength, t+n-1), 1:641, 1:1280))
        next!(updateprogress)
    end
end

function testread2(var,n)
    filename = "D:/datasets/era5/era5wind$(var)2016.h5"
    dataset = "wind$(var)"
    yearlength = 8760 + 24*leapyear(2016)
    for t=1:n:yearlength
        @time h5read(filename, dataset, (t:min(yearlength, t+n-1), 1:641, 1:1280))
    end
end

function testread3(var,n)
    filename = "D:/datasets/era5/era5wind$(var)2016.h5"
    dataname = "wind$(var)"
    yearlength = 8760 + 24*leapyear(2016)
    h5open(filename, "r") do file
        group = file["/"]
        dataset = group[dataname]
        for t=1:n:yearlength
            @time dataset[t:min(yearlength, t+n-1), 1:641, 1:1280]
        end       
    end
end

function testread4(n)
    filename = "D:/testwind.jld"
    yearlength = 8760 + 24*leapyear(2016)
    jldopen(filename, "r") do file
        for t=1:n:yearlength
            @time file["wind"][t:min(yearlength, t+n-1), 1:641, 1:1280]
        end       
    end
end

function testread5(n)
    filename = "D:/testwind.jld"
    yearlength = 8760 + 24*leapyear(2016)
    jldopen(filename, "r") do file
        for t=1:n:yearlength
            @time file["wind"][t:min(yearlength, t+n-1), 1:641, 1:1280]
        end       
    end
end

function testread5b(n)
    filename = "D:/testwind.jld"
    yearlength = 8760 + 24*leapyear(2016)
    jldopen(filename, "r") do file
        for lon=1:n:1280
            @time file["wind"][1:yearlength, 1:641, lon:min(1280, lon+n-1)]
        end       
    end
end

function testh5write(n)
    # r = repeat(rand(n,n,n), inner=(10,10,10))
    r = JLD.load("D:/testwind.jld", "wind");
    h5open("D:/rrr.h5", "w") do file
        group = file["/"]
        @time group["r", "compress", 3] = r
    end
end

testh5read() = h5read("D:/rrr.h5","r")

function testread6(n)
    filename = "D:/rrr.h5"
    dataset = "r"
    yearlength = 8760 + 24*leapyear(2016)
    for t=1:n:yearlength
        @time h5read(filename, dataset, (t:min(yearlength, t+n-1), 1:641, 1:1280))
    end
end

function testread6b(n)
    filename = "D:/rrr.h5"
    dataset = "r"
    yearlength = 8760 + 24*leapyear(2016)
    for lon=1:n:1280
        @time h5read(filename, dataset, (1:yearlength, 1:641, lon:min(1280, lon+n-1)))
    end
end

function testread7(n)
    filename = "D:/rrr.h5"
    dataset = "r"
    yearlength = 8760 + 24*leapyear(2016)
    h5open(filename, "r") do file
        group = file["/"]
        dataset = group[dataname]
        for t=1:n:yearlength
            @time dataset[t:min(yearlength, t+n-1), 1:641, 1:1280]
        end       
    end
end
