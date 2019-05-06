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
		GISsolar(gisregion=reg, pvclasses_min=minclasses, pvclasses_max=(minclasses.+step),
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
		GISwind(gisregion=reg, onshoreclasses_min=minclasses, onshoreclasses_max=(minclasses.+step),
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


plotly()

cap, lc = calcvars("China6")

plottriplepv(cap, lc, :low, 75)
plottriplewind(cap, lc, :low, 75)
plotmono(cap, lc)

nothing
