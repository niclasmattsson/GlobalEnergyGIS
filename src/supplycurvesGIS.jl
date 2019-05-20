using GlobalEnergyGIS, Plots, StatsBase

CRF(r,T) = r / (1 - 1/(1+r)^T)
meandrop(x; dims=dims) = dropdims(mean(x, dims=dims), dims=dims)
sumdrop(x; dims=dims) = dropdims(sum(x, dims=dims), dims=dims)

function mean_skipNaN_dim12(xx)
	hours, regs, classes = size(xx)
	out = Vector{Float64}(undef,classes)
	for c = 1:classes
		hourtot = 0.0
		hourcount = 0
		for h = 1:hours
			regtot = 0.0
			regcount = 0
			for r = 1:regs
				x = xx[h,r,c]
				if !isnan(x)
					regtot += x
					regcount += 1
				end
			end
			if regcount != 0
				hourtot += regtot / regcount
				hourcount += 1
			end
		end
		out[c] = hourtot / hourcount
	end
	return out
end

# for plant_area in [0.03, 0.12], persons_per_km2 in [75, 150]

# end

function supplycurves_pv(reg, plant_area, persons_per_km2, minclasses, step)
	cfr,cfpva,cfpvb,cfcsa,cfcsb,pvr,pva,pvb,csa,csb =
		GISsolar(savetodisk=false, gisregion=reg, pvclasses_min=minclasses, pvclasses_max=(minclasses.+step),
			plant_area=plant_area, plant_persons_per_km2=persons_per_km2)

	# results=loadresults("regionset=eurasia21, carboncap=0",resultsfile="results_prev.jld2");
	# sum(results.params[:demand][1:8,:]), sum(results.params[:demand][16:21,:])
	regdemand = Dict("Europe8" => 4.4046006293422505e6, "China6" => 1.2380343090170678e7)		# GWh/year
	demand = regdemand[reg]
	meanCF_a = mean_skipNaN_dim12(cfpva)
	meanCF_b = mean_skipNaN_dim12(cfpvb)
	meanCF_r = mean_skipNaN_dim12(cfr)
	capacity_a = sumdrop(pva, dims=1)
	capacity_b = sumdrop(pvb, dims=1)
	capacity_r = sumdrop(pvr, dims=1)
	# investcost = Dict(:pv => 800, :wind => 1200)
	# fixedcost = Dict(:pv => 16, :wind => 36)
	# lifetime = Dict(:pv => 25, :wind => 25)
	totalcost_a = capacity_a .* (800 * CRF(0.05, 25) + 16)
	totalcost_b = capacity_b .* (1.1*800 * CRF(0.05, 25) + 16)
	totalcost_r = capacity_r .* (1200 * CRF(0.05, 25) + 24)
	annualenergy_a = capacity_a .* meanCF_a .* 8760
	annualenergy_b = capacity_b .* meanCF_b .* 8760
	annualenergy_r = capacity_r .* meanCF_r .* 8760
	lcoe_a = totalcost_a ./ annualenergy_a .* 1000
	lcoe_b = totalcost_b ./ annualenergy_b .* 1000
	lcoe_r = totalcost_r ./ annualenergy_r .* 1000
	demandshare_a = annualenergy_a / demand
	demandshare_b = annualenergy_b / demand
	demandshare_r = annualenergy_r / demand
	return demandshare_a, demandshare_b, demandshare_r, lcoe_a, lcoe_b, lcoe_r
end

function supplycurves_wind(reg, area_onshore, persons_per_km2, minclasses, step)
	cfa,cfb,cfoff,wa,wb,woff =
		GISwind(savetodisk=false, gisregion=reg, onshoreclasses_min=minclasses, onshoreclasses_max=(minclasses.+step),
			offshoreclasses_min=minclasses, offshoreclasses_max=(minclasses.+step),
			area_onshore=area_onshore, persons_per_km2=persons_per_km2)

	# results=loadresults("regionset=eurasia21, carboncap=0",resultsfile="results_prev.jld2");
	# sum(results.params[:demand][1:8,:]), sum(results.params[:demand][16:21,:])
	regdemand = Dict("Europe8" => 4.4046006293422505e6, "China6" => 1.2380343090170678e7)		# GWh/year
	demand = regdemand[reg]
	meanCF_a = mean_skipNaN_dim12(cfa)
	meanCF_b = mean_skipNaN_dim12(cfb)
	meanCF_off = mean_skipNaN_dim12(cfoff)
	capacity_a = sumdrop(wa, dims=1)
	capacity_b = sumdrop(wb, dims=1)
	capacity_off = sumdrop(woff, dims=1)
	totalcost_a = capacity_a .* (1200 * CRF(0.05, 25) + 36)
	totalcost_b = capacity_b .* (1.1*1200 * CRF(0.05, 25) + 36)
	totalcost_off = capacity_off .* (2000 * CRF(0.05, 25) + 60)
	annualenergy_a = capacity_a .* meanCF_a .* 8760
	annualenergy_b = capacity_b .* meanCF_b .* 8760
	annualenergy_off = capacity_off .* meanCF_off .* 8760
	lcoe_a = totalcost_a ./ annualenergy_a .* 1000
	lcoe_b = totalcost_b ./ annualenergy_b .* 1000
	lcoe_off = totalcost_off ./ annualenergy_off .* 1000
	demandshare_a = annualenergy_a / demand
	demandshare_b = annualenergy_b / demand
	demandshare_off = annualenergy_off / demand
	return demandshare_a, demandshare_b, demandshare_off, lcoe_a, lcoe_b, lcoe_off
end

function plotsupply(capac, lcoe; hold=false, options...)
	s = sortperm(lcoe)
	x = cumsum(capac[s])
	y = lcoe[s]
	if hold
		plot!(x,y; options...)
	else
		plot(x,y; options...)
	end
end

function plottriplepv(cap, lc, area, pop)
	demandshare_a, demandshare_b, demandshare_r = cap[:pv,area,pop]
	lcoe_a, lcoe_b, lcoe_r = lc[:pv,area,pop]
	plotsupply(demandshare_a, lcoe_a, label="PV class A", linewidth=3)
	plotsupply(demandshare_b, lcoe_b, label="PV class B", linewidth=3, hold=true)
	display(plotsupply(demandshare_r, lcoe_r, label="PV rooftop", linewidth=3, hold=true,
					size=(1000,500), xlim=(0,1.2), ylim=(0,100), title="Solar PV",
					xlabel="Potential supply (relative to annual electricity demand in 2050)", ylabel="€/MWh",
					tickfont=16, legendfont=16, guidefont=16, titlefont=20))
end

function plottriplewind(cap, lc, area, pop)
	demandshare_wa, demandshare_wb, demandshare_off = cap[:wind,area,pop]
	lcoe_wa, lcoe_wb, lcoe_off = lc[:wind,area,pop]
	plotsupply(demandshare_wa, lcoe_wa, label="Onshore class A", linewidth=3)
	plotsupply(demandshare_wb, lcoe_wb, label="Onshore class B", linewidth=3, hold=true)
	display(plotsupply(demandshare_off, lcoe_off, label="Offshore", linewidth=3, hold=true,
					size=(1000,500), xlim=(0,1.2), ylim=(0,100), title="Wind",
					xlabel="Potential supply (relative to annual electricity demand in 2050)", ylabel="€/MWh",
					tickfont=16, legendfont=16, guidefont=16, titlefont=20))
end

function mergevars(cap, lc)
	mergecap, mergelc = Dict(), Dict()
	for tech = [:pv, :wind], area = [:low, :high], pop = [75]
		cc,ll = cap[tech,area,pop], lc[tech,area,pop]
		mergecap[tech,area,pop] = [cc[1]; cc[2]; cc[3]]
		mergelc[tech,area,pop] = [ll[1]; ll[2]; ll[3]]	
	end
	return mergecap, mergelc
end

function plotmonopv(cap, lc)
	mcap, mlc = mergevars(cap, lc)
	plotsupply(mcap[:pv,:low,75], mlc[:pv,:low,75], label="area=3%", linewidth=3)
	display(plotsupply(mcap[:pv,:high,75], mlc[:pv,:high,75], label="area=12%", linewidth=3, hold=true,
					size=(1000,500), xlim=(0,1.2), ylim=(0,100), title="Solar PV (class A+B+rooftop)",
					xlabel="Potential supply (relative to annual electricity demand in 2050)", ylabel="€/MWh",
					tickfont=16, legendfont=16, guidefont=16, titlefont=20))
end

function plotmonowind(cap, lc)
	mcap, mlc = mergevars(cap, lc)
	plotsupply(mcap[:wind,:low,75], mlc[:wind,:low,75], label="area=5%", linewidth=3)
	display(plotsupply(mcap[:wind,:high,75], mlc[:wind,:high,75], label="area=20%", linewidth=3, hold=true,
					size=(1000,500), xlim=(0,1.2), ylim=(0,100), title="Wind (class A+B+offshore)",
					xlabel="Potential supply (relative to annual electricity demand in 2050)", ylabel="€/MWh",
					tickfont=16, legendfont=16, guidefont=16, titlefont=20))
end

function plotmono(cap, lc)
	mcap, mlc = mergevars(cap, lc)
	plotsupply(mcap[:pv,:low,75], mlc[:pv,:low,75], label="PV area=3%", line=(3, :solid, RGB([240,220,0]/255...)))
	plotsupply(mcap[:wind,:low,75], mlc[:wind,:low,75], label="wind area=5%", line=(3, :solid, RGB([149,179,215]/255...)), hold=true)
	plotsupply(mcap[:pv,:high,75], mlc[:pv,:high,75], label="PV area=12%", line=(3, :dash, RGB([240,220,0]/255...)), hold=true)
	display(plotsupply(mcap[:wind,:high,75], mlc[:wind,:high,75], label="wind area=20%", line=(3, :dash, RGB([149,179,215]/255...)), hold=true,
					size=(1000,500), xlim=(0,1.2), ylim=(0,100), title="Solar PV & wind (class A+B+rooftop+offshore)",
					xlabel="Potential supply (relative to annual electricity demand in 2050)", ylabel="€/MWh",
					tickfont=16, legendfont=16, guidefont=16, titlefont=20))
end

function calcvars(reg)
	cap, lc = Dict(), Dict()
	area = :low
	pop = 75

	step = 0.00005
	plantarea = area == :low ? 0.03 : 0.12
	demandshare_a, demandshare_b, demandshare_r, lcoe_a, lcoe_b, lcoe_r = supplycurves_pv(reg, plantarea, pop, [0:step:0.311-step;], step)
	cap[:pv,area,pop] = (demandshare_a, demandshare_b, demandshare_r)
	lc[:pv,area,pop] = (lcoe_a, lcoe_b, lcoe_r)

	onshorearea = area == :low ? 0.05 : 0.20
	step = 0.002
	demandshare_wa, demandshare_wb, demandshare_off, lcoe_wa, lcoe_wb, lcoe_off = supplycurves_wind(reg, onshorearea, pop, [0:step:25-step;], step)
	cap[:wind,area,pop] = (demandshare_wa, demandshare_wb, demandshare_off)
	lc[:wind,area,pop] = (lcoe_wa, lcoe_wb, lcoe_off)

	cap[:pv,:high,pop] = cap[:pv,:low,pop] .* 4
	lc[:pv,:high,pop] = lc[:pv,:low,pop]
	cap[:wind,:high,pop] = cap[:wind,:low,pop] .* 4
	lc[:wind,:high,pop] = lc[:wind,:low,pop]

	return cap, lc
end

function plotwindturbinecurve()
	speeds = powercurves[:,1]
	cfturbine = powercurves[:,2]
	cfpark = powercurves[:,3]
	display(plot(speeds, cfturbine, line=(3, :dash, RGB([149,179,215]/255...)), label="single turbine"))
	display(plot!(speeds, cfpark, line=(3, :solid, RGB([149,179,215]/255...)), label="park output",
					size=(1000,500), xlim=(0,32), ylim=(0,1.1), title="Wind turbine and park output",
					xlabel="Wind speed (m/s)", ylabel="Instantaneous capacity factor",
					tickfont=18, legendfont=18, guidefont=18, titlefont=20))
end


plotly()

cap, lc = calcvars("China6")

plottriplepv(cap, lc, :low, 75)
plottriplewind(cap, lc, :low, 75)
plotmono(cap, lc)


const powercurves = [
         0         0         0
    0.1000         0    0.0000
    0.2000         0    0.0000
    0.3000         0    0.0000
    0.4000         0    0.0000
    0.5000         0    0.0001
    0.6000         0    0.0001
    0.7000         0    0.0001
    0.8000         0    0.0002
    0.9000         0    0.0003
    1.0000         0    0.0003
    1.1000         0    0.0005
    1.2000         0    0.0006
    1.3000         0    0.0008
    1.4000         0    0.0010
    1.5000         0    0.0013
    1.6000         0    0.0016
    1.7000         0    0.0020
    1.8000         0    0.0025
    1.9000         0    0.0030
    2.0000         0    0.0036
    2.1000         0    0.0044
    2.2000         0    0.0052
    2.3000         0    0.0062
    2.4000         0    0.0073
    2.5000         0    0.0085
    2.6000         0    0.0098
    2.7000         0    0.0114
    2.8000         0    0.0130
    2.9000         0    0.0148
    3.0000         0    0.0168
    3.1000    0.0106    0.0190
    3.2000    0.0147    0.0213
    3.3000    0.0187    0.0237
    3.4000    0.0229    0.0263
    3.5000    0.0270    0.0291
    3.6000    0.0312    0.0321
    3.7000    0.0354    0.0352
    3.8000    0.0397    0.0385
    3.9000    0.0440    0.0419
    4.0000    0.0484    0.0455
    4.1000    0.0529    0.0492
    4.2000    0.0574    0.0531
    4.3000    0.0621    0.0572
    4.4000    0.0668    0.0615
    4.5000    0.0716    0.0659
    4.6000    0.0766    0.0705
    4.7000    0.0816    0.0753
    4.8000    0.0868    0.0802
    4.9000    0.0921    0.0854
    5.0000    0.0976    0.0908
    5.1000    0.1032    0.0964
    5.2000    0.1090    0.1023
    5.3000    0.1150    0.1084
    5.4000    0.1214    0.1147
    5.5000    0.1281    0.1213
    5.6000    0.1353    0.1282
    5.7000    0.1429    0.1354
    5.8000    0.1510    0.1428
    5.9000    0.1597    0.1506
    6.0000    0.1691    0.1587
    6.1000    0.1791    0.1670
    6.2000    0.1898    0.1757
    6.3000    0.2011    0.1847
    6.4000    0.2129    0.1940
    6.5000    0.2253    0.2037
    6.6000    0.2380    0.2136
    6.7000    0.2512    0.2238
    6.8000    0.2647    0.2344
    6.9000    0.2786    0.2453
    7.0000    0.2927    0.2564
    7.1000    0.3070    0.2679
    7.2000    0.3216    0.2796
    7.3000    0.3363    0.2917
    7.4000    0.3514    0.3040
    7.5000    0.3666    0.3166
    7.6000    0.3822    0.3294
    7.7000    0.3980    0.3425
    7.8000    0.4141    0.3558
    7.9000    0.4305    0.3694
    8.0000    0.4473    0.3832
    8.1000    0.4643    0.3972
    8.2000    0.4817    0.4114
    8.3000    0.4994    0.4257
    8.4000    0.5175    0.4402
    8.5000    0.5360    0.4549
    8.6000    0.5548    0.4697
    8.7000    0.5740    0.4845
    8.8000    0.5937    0.4994
    8.9000    0.6137    0.5144
    9.0000    0.6341    0.5293
    9.1000    0.6550    0.5442
    9.2000    0.6762    0.5590
    9.3000    0.6975    0.5738
    9.4000    0.7188    0.5884
    9.5000    0.7400    0.6028
    9.6000    0.7609    0.6170
    9.7000    0.7814    0.6309
    9.8000    0.8014    0.6446
    9.9000    0.8206    0.6579
   10.0000    0.8390    0.6709
   10.1000    0.8565    0.6835
   10.2000    0.8729    0.6958
   10.3000    0.8884    0.7076
   10.4000    0.9028    0.7189
   10.5000    0.9161    0.7299
   10.6000    0.9283    0.7403
   10.7000    0.9394    0.7503
   10.8000    0.9494    0.7599
   10.9000    0.9582    0.7690
   11.0000    0.9659    0.7776
   11.1000    0.9723    0.7858
   11.2000    0.9777    0.7935
   11.3000    0.9821    0.8009
   11.4000    0.9856    0.8078
   11.5000    0.9884    0.8143
   11.6000    0.9906    0.8205
   11.7000    0.9922    0.8264
   11.8000    0.9934    0.8319
   11.9000    0.9944    0.8371
   12.0000    0.9951    0.8420
   12.1000    0.9958    0.8467
   12.2000    0.9964    0.8512
   12.3000    0.9970    0.8554
   12.4000    0.9975    0.8595
   12.5000    0.9980    0.8633
   12.6000    0.9984    0.8670
   12.7000    0.9988    0.8706
   12.8000    0.9992    0.8740
   12.9000    0.9996    0.8774
   13.0000    1.0000    0.8806
   13.1000    1.0000    0.8837
   13.2000    1.0000    0.8867
   13.3000    1.0000    0.8896
   13.4000    1.0000    0.8925
   13.5000    1.0000    0.8953
   13.6000    1.0000    0.8980
   13.7000    1.0000    0.9006
   13.8000    1.0000    0.9031
   13.9000    1.0000    0.9056
   14.0000    1.0000    0.9080
   14.1000    1.0000    0.9103
   14.2000    1.0000    0.9126
   14.3000    1.0000    0.9147
   14.4000    1.0000    0.9167
   14.5000    1.0000    0.9187
   14.6000    1.0000    0.9205
   14.7000    1.0000    0.9223
   14.8000    1.0000    0.9239
   14.9000    1.0000    0.9254
   15.0000    1.0000    0.9268
   15.1000    1.0000    0.9281
   15.2000    1.0000    0.9293
   15.3000    1.0000    0.9304
   15.4000    1.0000    0.9314
   15.5000    1.0000    0.9323
   15.6000    1.0000    0.9330
   15.7000    1.0000    0.9337
   15.8000    1.0000    0.9343
   15.9000    1.0000    0.9349
   16.0000    1.0000    0.9353
   16.1000    1.0000    0.9357
   16.2000    1.0000    0.9361
   16.3000    1.0000    0.9363
   16.4000    1.0000    0.9366
   16.5000    1.0000    0.9368
   16.6000    1.0000    0.9369
   16.7000    1.0000    0.9371
   16.8000    1.0000    0.9372
   16.9000    1.0000    0.9373
   17.0000    1.0000    0.9373
   17.1000    1.0000    0.9374
   17.2000    1.0000    0.9374
   17.3000    1.0000    0.9374
   17.4000    1.0000    0.9375
   17.5000    1.0000    0.9375
   17.6000    1.0000    0.9375
   17.7000    1.0000    0.9375
   17.8000    1.0000    0.9375
   17.9000    1.0000    0.9375
   18.0000    1.0000    0.9375
   18.1000    1.0000    0.9375
   18.2000    1.0000    0.9375
   18.3000    1.0000    0.9375
   18.4000    1.0000    0.9375
   18.5000    1.0000    0.9375
   18.6000    1.0000    0.9375
   18.7000    1.0000    0.9375
   18.8000    1.0000    0.9375
   18.9000    1.0000    0.9375
   19.0000    1.0000    0.9375
   19.1000    1.0000    0.9375
   19.2000    1.0000    0.9375
   19.3000    1.0000    0.9375
   19.4000    1.0000    0.9375
   19.5000    1.0000    0.9375
   19.6000    1.0000    0.9375
   19.7000    1.0000    0.9375
   19.8000    1.0000    0.9375
   19.9000    1.0000    0.9375
   20.0000    1.0000    0.9375
   20.1000    1.0000    0.9375
   20.2000    1.0000    0.9375
   20.3000    1.0000    0.9375
   20.4000    1.0000    0.9375
   20.5000    1.0000    0.9375
   20.6000    1.0000    0.9375
   20.7000    1.0000    0.9375
   20.8000    1.0000    0.9375
   20.9000    1.0000    0.9375
   21.0000    1.0000    0.9375
   21.1000    1.0000    0.9375
   21.2000    1.0000    0.9375
   21.3000    1.0000    0.9375
   21.4000    1.0000    0.9375
   21.5000    1.0000    0.9375
   21.6000    1.0000    0.9375
   21.7000    1.0000    0.9375
   21.8000    1.0000    0.9375
   21.9000    1.0000    0.9375
   22.0000    1.0000    0.9375
   22.1000    1.0000    0.9375
   22.2000    1.0000    0.9374
   22.3000    1.0000    0.9372
   22.4000    1.0000    0.9370
   22.5000    1.0000    0.9366
   22.6000    1.0000    0.9361
   22.7000    1.0000    0.9353
   22.8000    1.0000    0.9343
   22.9000    1.0000    0.9329
   23.0000    1.0000    0.9311
   23.1000    1.0000    0.9289
   23.2000    1.0000    0.9260
   23.3000    1.0000    0.9226
   23.4000    1.0000    0.9183
   23.5000    1.0000    0.9131
   23.6000    1.0000    0.9068
   23.7000    1.0000    0.8994
   23.8000    1.0000    0.8906
   23.9000    1.0000    0.8803
   24.0000    1.0000    0.8684
   24.1000    1.0000    0.8547
   24.2000    1.0000    0.8392
   24.3000    1.0000    0.8216
   24.4000    1.0000    0.8020
   24.5000    1.0000    0.7804
   24.6000    1.0000    0.7566
   24.7000    1.0000    0.7308
   24.8000    1.0000    0.7030
   24.9000    1.0000    0.6734
   25.0000    1.0000    0.6421
   25.1000    0.9000    0.6093
   25.2000    0.8000    0.5753
   25.3000    0.7000    0.5403
   25.4000    0.6000    0.5047
   25.5000    0.5000    0.4688
   25.6000    0.4000    0.4328
   25.7000    0.3000    0.3972
   25.8000    0.2000    0.3622
   25.9000    0.1000    0.3282
   26.0000         0    0.2954
   26.1000         0    0.2641
   26.2000         0    0.2345
   26.3000         0    0.2067
   26.4000         0    0.1809
   26.5000         0    0.1571
   26.6000         0    0.1355
   26.7000         0    0.1159
   26.8000         0    0.0983
   26.9000         0    0.0828
   27.0000         0    0.0691
   27.1000         0    0.0572
   27.2000         0    0.0469
   27.3000         0    0.0381
   27.4000         0    0.0307
   27.5000         0    0.0244
   27.6000         0    0.0192
   27.7000         0    0.0149
   27.8000         0    0.0115
   27.9000         0    0.0086
   28.0000         0    0.0064
   28.1000         0    0.0046
   28.2000         0    0.0032
   28.3000         0    0.0022
   28.4000         0    0.0014
   28.5000         0    0.0009
   28.6000         0    0.0005
   28.7000         0    0.0003
   28.8000         0    0.0001
   28.9000         0    0.0000
   29.0000         0         0
   29.1000         0         0
   29.2000         0         0
   29.3000         0         0
   29.4000         0         0
   29.5000         0         0
   29.6000         0         0
   29.7000         0         0
   29.8000         0         0
   29.9000         0         0
   30.0000         0         0
   30.1000         0         0
   30.2000         0         0
   30.3000         0         0
   30.4000         0         0
   30.5000         0         0
   30.6000         0         0
   30.7000         0         0
   30.8000         0         0
   30.9000         0         0
   31.0000         0         0
   31.1000         0         0
   31.2000         0         0
   31.3000         0         0
   31.4000         0         0
   31.5000         0         0
   31.6000         0         0
   31.7000         0         0
   31.8000         0         0
   31.9000         0         0
   32.0000         0         0
   32.1000         0         0
   32.2000         0         0
   32.3000         0         0
   32.4000         0         0
   32.5000         0         0
   32.6000         0         0
   32.7000         0         0
   32.8000         0         0
   32.9000         0         0
   33.0000         0         0
   33.1000         0         0
   33.2000         0         0
   33.3000         0         0
   33.4000         0         0
   33.5000         0         0
   33.6000         0         0
   33.7000         0         0
   33.8000         0         0
   33.9000         0         0
   34.0000         0         0
   34.1000         0         0
   34.2000         0         0
   34.3000         0         0
   34.4000         0         0
   34.5000         0         0
   34.6000         0         0
   34.7000         0         0
   34.8000         0         0
   34.9000         0         0
   35.0000         0         0
]

nothing