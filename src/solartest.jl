
# download ERA5 TISR variable (top of atmosphere incident solar radiation)
# GlobalEnergyGIS.request_era5_vars("D:/testera5/tisrtest.nc", ["TOA_incident_solar_radiation"], "2018-06-01", "2018-06-15")


using NCDatasets, Dates, StatsBase, Plots, GlobalEnergyGIS
# https://www.itacanet.org/the-sun-as-a-source-of-energy/part-2-solar-energy-reaching-the-earths-surface/#2.1.-The-Solar-Constant
solarconstant(dt::DateTime) = 1361 * (1 + 0.034*cos(2*π*(dt - DateTime(year(dt))).value/1000/3600/24/365.25))

# solar constant info
# best historical record:
# https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/2016GL071866
# variability on different time scales:
# https://www.swsc-journal.org/articles/swsc/abs/2016/01/swsc160010/swsc160010.html



# Check if TISR is close to calculated insolation using solar constant projected to latitude/longitude
# location using solar position function. Note that while the solar position calculations are instantaneous
# positions, ERA5 radiation variables represent accumulated radiation over the hour *ending* at the
# indicated time. Therefore, solar positions must be shifted 30 minutes BACK to correspond to the midpoint
# time of the ERA5 accumulations.
# ERA5 "accumulations are over the hour ending at the forecast step":
# https://confluence.ecmwf.int//display/CKB/ERA5+data+documentation#ERA5datadocumentation-Meanratesandaccumulations
res = 0.28125
lons = (-180+res/2:res:180-res/2)
lats = (90-res/2:-res:-90+res/2)
println("Loading tisrtest.nc...")
ncdataset = Dataset("D:/testera5/tisrtest.nc")
dt = ncdataset["time"][1] - Minute(30)
tisr = replace(ncdataset["tisr"][:,:,:], missing => 0.0)
toa = tisr[:,:,1] ./ 3600
zen = acosd.(toa ./ solarconstant(dt))
zen2 = [min(90, solarposition(dt,lat,lon)[1]) for lon in lons, lat in lats]
toa2 = solarconstant(dt) .* cosd.(zen2)
plotly()
display(heatmap(toa',yflip=true))
display(heatmap(toa2',yflip=true))
f = zen.<85
mean(abs.(zen[f].-zen2[f]))  # absolute zenith error only 0.1 degrees off on average