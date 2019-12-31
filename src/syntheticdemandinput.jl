function readsspdemand()
    return # Dict((region,model,year) => sspdemand)
end

# Can call this with ssp5 or ssp32 (global constants)
function make_sspregionlookup(ssp)
    _, _, regionlist, _, _ = loadregions("Global_GADM0")
    sspregions = fill("", length(regionlist))
    for (sspregionname, gadmregionnames) in ssp
        for gname in gadmregionnames
            index = findfirst(isequal(gname), regionlist)
            if index == nothing
                error("Region $gname is not a GADM level 0 region.")
            end
            sspregions[index] = sspregionname
        end
    end
    return sspregions # Vector{String}(length numcountries, indexed by gadm country code)
end

# load population dataset and sum globally by country
# Vector{Float64}(length numcountries, indexed by gadm country code)
function nationalpopulation(scenarioyear)
    pop = JLD.load("population_$scenarioyear.jld", "population")
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions("Global_GADM0")
    natpop = zeros(length(regionlist))
    for j = 1:latrange
        for i = 1:lonrange
            reg = regions[i,j]
            if reg > 0
                natpop[reg] += pop[i,j]
            end
        end
    end
    return natpop
end

function ieademandpercapita(scenarioyear)
    ieademand = JLD.load("ieademand2017.jld", "ieademand")  # TWh/year
    natpop = nationalpopulation(scenarioyear)               # million people
    return ieademand./natpop./1e6                           # MWh/year/capita   
end

function calcdemandmultipliers()
    sspdemand = readsspdemand()
    demandmult2050 = Dict{String, Float64}()
    for sspreg in keys(ssp5)
        demandmult2050image = sspdemand[sspreg,"image",2050] / sspdemand[sspreg,"image",2020]
        demandmult2050message = sspdemand[sspreg,"message",2050] / sspdemand[sspreg,"message",2020]
        demandmult2050[sspreg] = ((demandmult2050image + demandmult2050message)/2) ^ ((2050-2017)/(2050-2020))
    end
    return demandmult2050
end

function makesyntheticdemandinput(; optionlist...)
    println("Add countries (matrix of codes) and countryname vector to makeregions()!")
    println("Beware of readdlm() - there are 81 subregions called Nan, these become the non-number NaN. Fix in ogr2ogr.")
    options = WindOptions(merge(windoptions(), optionlist))
    regions, _, regionlist, _, pop, _, _, _, lonrange, latrange = read_datasets(options)

    @unpack scenarioyear, res = options
    gdp = JLD.load("gdp_$(scenarioyear).jld", "gdp")[lonrange,latrange]
    country = JLD.load("countries.jld", "countrycode")[lonrange,latrange]

    cellarea = rastercellarea.(latrange, res)
    popdens = pop ./ cellarea'

    # from IEA online stats (field Domestic Supply): 
    # https://www.iea.org/statistics/?country=MOROCCO&year=2016&category=Electricity&indicator=ElecGenByFuel&mode=table&dataTable=ELECTRICITYANDHEAT
    demandpercapita = ieademandpercapita(scenarioyear)  # vector of current national demand per capita (2017), MWh/year/capita
    ssp5region = make_sspregionlookup(ssp5)             # vector of SSP 5-region names for each country
    demandmult2050 = calcdemandmultipliers()            # Dict: SSP region name => multiplier 

    demand = zeros(length(regionlist))

    nlons, nlats = size(regions)
    for j = 1:nlats
        for i = 1:nlons
            reg = regions[i,j]
            countrycode = country[i,j]
            if reg > 0 && countrycode > 0
                sspreg = ssp5region[countrycode]
                isempty(sspreg) && error("Oops, no SSP region assigned to region $(regionlist[countrycode]).")
                demand[reg] += demandpercapita[countrycode] * pop[i,j]/1e6 * demandmult2050[sspreg]     # TWh/year
            end
        end
    end

    return regions, regionlist, demand, pop, popdens, gdp
    #save('tillville_eu10.mat', 'regions', 'regionlist', 'demand', 'pop', 'popdens', 'gdp');

end


    # ieademand2017 = [144.3+133.2+87.7+35.6, 326.8+2.4, 514.6, 598.6+120.1+91.7+8.5, 356.9+29.7,
    #                 63.2+38.9+15.3+18.4+14.0+37.4+3.4+5.5+7.7+7.7, 10.1+7.5+12.5, 168.6, 282.4+55.2+0.3,
    #                 75.5+67.1+72.3+44.6+29.7+60.1]
    # bpdemand2017 = [163.9+148.7+67.9+30.9, 295.4, 554.1, 654.2+86.4+116.6+2.2, 335.9+30.8,
    #                 55.1+45.8+16.3+12.0+5.6+0.7*84.5, 13.4+7.5+4.2, 170.3, 275.4+60.0,
    #                 70.1+63.2+87.0+32.9+27.6+63.6]

    # Final energy - electricity 2050 relative to 2020 in SSP2-34
    # [oecd] 
    # demandmult2050image = 36.240/32.322
    # demandmult2050message = 47.006/34.666

    # sspdemand = readsspdemand()
    # demandmult2050image = sspdemand["oecd","image",2050] / sspdemand["oecd","image",2020]
    # demandmult2050message = sspdemand["oecd","message",2050] / sspdemand["oecd","message",2020]
    # demandmult2050 = (demandmult2050image + demandmult2050message)/2

    # demand = ieademand2017 * demandmult2050^((2050-2017)/(2050-2020))




# load regions_Eurasia38_12.mat;
# regions38 = regions;
# regionlist38 = regionlist;

# load popdens.mat;
# popdens = pop;
# load popcount.mat;
# pop = uint32(pop);
# load gdp.mat;
# load regions_Eurasia21_12.mat;

# clear R;

# % TWh/year
# bpdemand2017 = [163.9+148.7+67.9+30.9  554.1  654.2+86.4+116.6+2.2  335.9+30.8  ...
#              55.1+45.8+295.4+16.3+12.0+5.6+0.7*84.5  170.3+13.4+7.5+4.2  275.4+60.0 ...
#              70.1+63.2+87.0+32.9+27.6+63.6  34.3+157.1  295.5+24.3+0.5*46.6  103.0 ...
#              22.8+60.1+0.5*46.6  1091.2 0 0 6495.1 0 0 0 0 0];

# % Final energy - electricity 2050 relative to 2020 in SSP2-34
# % [asia oecd ref] 
# demandmult2050image = [61.056/34.066  36.240/32.322  5.446/4.590];
# demandmult2050message = [62.145/30.766  47.006/34.666  6.767/4.511];
# demandmult2050 = (demandmult2050image + demandmult2050message)/2;

# demand = bpdemand2017;
# demand([1:8 10]) = demand([1:8 10]) * demandmult2050(2);
# demand([9 11 12 13]) = demand([9 11 12 13]) * demandmult2050(3);
# demand(16) = demand(16) * demandmult2050(1);

# demandrussia = demand(13);
# poprussia = sum(sum(pop(regions38>=24 & regions38<=30)));
# for i=13:15
#     demand(i) = demandrussia / poprussia * sum(sum(pop(regions==i)));
# end

# demandchina = demand(16);
# popchina = sum(sum(pop(regions38>=31 & regions38<=36)));
# for i=16:21
#     poprussiareg(i) = sum(sum(pop(regions==i)));
#     demand(i) = demandchina / popchina * sum(sum(pop(regions==i)));
# end

# save('tillville.mat', 'regions', 'regionlist', 'demand', 'pop', 'popdens', 'gdp');