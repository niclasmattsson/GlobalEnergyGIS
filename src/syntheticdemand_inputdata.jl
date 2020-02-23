export makesyntheticdemandinput

# Can call this with ssp5 or ssp32 (global constants)
# Since not all GADM level 0 regions are assigned to SSP regions, some regions are missed
# (mostly small island nations and the like). This captures 99.85% of global population.
function make_sspregionlookup(ssp)
    println("Lookup SSP R5.2 region names for each GADM level 0 region...")
    _, _, regionlist, _, _ = loadregions("Global_GADM0")
    sspregions = fill("", length(regionlist))
    for (sspregionname, gadmregionnames) in ssp
        for gname in gadmregionnames
            index = findfirst(isequal(Symbol(gname)), regionlist)
            if index == nothing
                error("Region $gname is not a GADM level 0 region.")
            end
            sspregions[index] = sspregionname
        end
    end
    return sspregions # Vector{String}(length numcountries, indexed by gadm country code)
end

function getnationalpopulation(scenarioyear)
    datafolder = getconfig("datafolder")
    filename = joinpath(datafolder, "nationalpopulation_$scenarioyear.jld")
    natpop = isfile(filename) ? JLD.load(filename, "natpop") : savenationalpopulation(scenarioyear)
    return natpop
end

# load population dataset and sum globally by country
# Vector{Float64}(length numcountries, indexed by gadm country code)
function savenationalpopulation(scenarioyear)
    println("Calculating population in all GADM level 0 regions...")
    datafolder = getconfig("datafolder")
    pop = JLD.load(joinpath(datafolder, "population_$scenarioyear.jld"), "population")   # persons per pixel
    regions, _, regionlist, lonrange, latrange = loadregions("Global_GADM0")
    natpop = zeros(length(regionlist))
    for j in latrange
        for i in lonrange
            reg = regions[i,j]
            if reg > 0
                natpop[reg] += pop[i,j]
            end
        end
    end
    JLD.save(joinpath(datafolder,"nationalpopulation_$scenarioyear.jld"), "natpop", natpop, compress=true)
    return natpop
end

# Get current national electricity demand from IEA online statistics (actually "domestic supply").
# https://www.iea.org/statistics/?country=MOROCCO&year=2016&category=Electricity&indicator=ElecGenByFuel&mode=table&dataTable=ELECTRICITYANDHEAT
function ieademand()
    println("Get current national electricity demand from IEA statistics...")
    datafolder = getconfig("datafolder")
    iea = CSV.read(joinpath(datafolder, "ieademand_2016.csv"))      # GWh/year
    _, _, regionlist, _, _ = loadregions("Global_GADM0")
    nationaldemand = zeros(length(regionlist))
    for row in eachrow(iea)
        country, demand = row
        country == uppercase(country) && continue   # skip IEA aggregated regions in uppercase, e.g. "ASIA"
        f = findfirst(regionlist .== Symbol(country))
        f == nothing && error("IEA country $country missing in GADM regions.")
        nationaldemand[f] = demand/1000
    end
    return nationaldemand        # TWh/year 
end

# Lookup electricity demand from a SSP database dataframe. Return a single dataframe row.
ssplookup(ssp, model, scen, region) =
    ssp[(ssp[:,:MODEL] .== model) .& (ssp[:,:SCENARIO] .== scen) .& (ssp[:,:REGION] .== region), :][1,:]

# Calculate regional multiplers for demand from 2016 to target year
function calcdemandmultipliers(scenarioyear)
    println("Calculate demand multipliers...")
    # First read electricity demand from the SSP database into a dataframe.
    datafolder = getconfig("datafolder")
    ssp = CSV.read(joinpath(datafolder, "SSP v2 Final Energy - Electricity.csv"))

    # SSP scenarios also include radiative forcing targets, e.g. SSP2-34
    # (only for IAM energy system results, not underlying population & GDP scenario)
    # We'll hardcode the 3.4 W/m2 scenario variant for now and make it a configurable option later. 
    scen = uppercase("$(scenarioyear[1:4])-34")
    year = parse(Int, scenarioyear[6:9])

    # We'll take the average result from IMAGE and MESSAGE models.
    demandmult2050 = Dict{String, Float64}()
    for sspreg in keys(ssp5)
        reg = "R5.2$sspreg"
        image = ssplookup(ssp, "IMAGE", scen, reg)
        message = ssplookup(ssp, "MESSAGE-GLOBIOM", scen, reg)
        # Take the year 2020 to represent 2016, then adjust for this difference.
        image_mult = image[Symbol(year)] / image[Symbol(2020)]
        message_mult = message[Symbol(year)] / message[Symbol(2020)]
        demandmult2050[sspreg] = ((image_mult + message_mult)/2) ^ ((year-2016)/(year-2020))
    end
    return demandmult2050
end

# function getsyntheticdemandregions()
#   datafolder = getconfig("datafolder")
#   filename = joinpath(datafolder, "regions_syntheticdemandregions.jld")
#   !isfile(filename) && saveregions("syntheticdemandregions", syntheticdemandregions)
#   demandregions, _, demandregionlist, _, _ = loadregions("syntheticdemandregions")
#     return demandregions, demandregionlist
# end

function makesyntheticdemandinput(; optionlist...)
    options = WindOptions(merge(windoptions(), optionlist))
    regions, _, regionlist, _, pop, _, _, _, lonrange, latrange = read_datasets(options)            # pop unit: people/grid cell

    @unpack scenarioyear, res, era_year, gisregion = options
    datafolder = getconfig("datafolder")
    gdp = JLD.load(joinpath(datafolder, "gdp_$(scenarioyear).jld"))["gdp"][lonrange,latrange]       # unit: USD(2010)/grid cell, PPP

    gadm0regions, _, gadm0regionlist, _, _ = loadregions("Global_GADM0")
    countrycodes = gadm0regions[lonrange,latrange]

    lats = (90-res/2:-res:-90+res/2)[latrange]
    cellarea = rastercellarea.(lats, res)
    popdens = pop ./ cellarea'

    natpop = getnationalpopulation(scenarioyear)        # people/GADMregion (vector, length numGADMcountries)
    demandpercapita = ieademand()./natpop * 1e6         # MWh/year/capita (vector, length numGADMcountries)
    ssp5region = make_sspregionlookup(ssp5)             # SSP 5-region names (vector, length numGADMcountries)
    demandmult = calcdemandmultipliers(scenarioyear)    # Dict: SSP region name => multiplier 

    numreg = length(regionlist)
    regionaldemand = zeros(numreg)
    regionalpop = zeros(numreg)
    regionalgdp = zeros(numreg)

    println("Calculating annual electricity demand for model regions...")
    nlons, nlats = size(regions)
    updateprogress = Progress(nlats, 1)
    for j = 1:nlats
        for i = 1:nlons
            reg = regions[i,j]
            countrycode = countrycodes[i,j]
            if reg > 0 && reg != NOREGION && countrycode > 0
                sspreg = ssp5region[countrycode]
                if isempty(sspreg)
                    # error("Oops, no SSP region assigned to region $(gadm0regionlist[countrycode]).")
                    continue
                end
                regionaldemand[reg] += demandpercapita[countrycode] * pop[i,j]/1e6 * demandmult[sspreg]     # TWh/year
                regionalpop[reg] += pop[i,j]            # unit: people
                regionalgdp[reg] += gdp[i,j]        # unit: USD(2010) 
            end
        end
        next!(updateprogress)
    end

    datafolder = getconfig("datafolder")
    syntheticdemanddata = joinpath(datafolder, "syntheticdemand", "data", "julia")
    mkpath(syntheticdemanddata)

    filename = joinpath(syntheticdemanddata, "regiondata_$(era_year)_$gisregion.h5")
    h5open(filename, "w") do file
        write(file, "regionlist", string.(regionlist))
        write(file, "regionaldemand", regionaldemand)
        write(file, "regionalpop", regionalpop)
        write(file, "regionalgdp", regionalgdp)
    end

    return regionlist, regionaldemand, regionalpop, regionalgdp
end
