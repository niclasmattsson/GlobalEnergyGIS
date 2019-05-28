function makedistances(gisregion; scenarioyear="ssp2_2050", res=0.01)
    regions, offshoreregions, regionlist, lonrange, latrange = loadregions(gisregion)
    numreg = length(regionlist)
    cellarea = rastercellarea.(latrange, res)

    pop = JLD.load("population_$scenarioyear.jld", "population")[lonrange,latrange]
    popdens = pop ./ cellarea'

    geocenters, popcenters = regioncenters(regions, numreg, popdens, lonrange, latrange)
    connected, connectedoffshore = connectedregions(regions, offshoreregions, numreg)
    distances = [greatcircledistance(popcenters[r1], popcenters[r2]) for r1 = 1:numreg, r2 = 1:numreg]
    regionpop = [sum(pop.==r) for r = 1:numreg]

    matopen("distances_$gisregion.mat", "w") do file
        write(file, "distances", distances)
        write(file, "connected", connected)
        write(file, "connectedoffshore", connectedoffshore)
        # write(file, "regionlist", regionlist)
        # write(file, "population", pop)
        # write(file, "demand", demand)
    end
end

# returns great circle distance in km between points given as (lat,lon) tuples (in degrees).
greatcircledistance(point1::Tuple, point2::Tuple) = haversine(point1, point2, 6371.0)

function regioncenters(regions, numreg, popdens, lonrange, latrange)
    geocenters = Vector{Tuple{Float64,Float64}}(undef, numreg)
    popcenters = Vector{Tuple{Float64,Float64}}(undef, numreg)
    for r = 1:numreg
        regmask = (regions.==r)
        regpopdens = popdens.*regmask
        geocenters[r] = sum(latrange.*regmask)/sum(regmask), sum(lonrange.*regmask)/sum(regmask)
        popcenters[r] = sum(latrange.*regpopdens)/sum(regpopdens), sum(lonrange.*regpopdens)/sum(regpopdens)
    end
    return geocenters, popcenters
end

function getneighbors(regions, r)
    regmask = (regions.==r)
    regmask_expanded = (imfilter(regmask, diskfilterkernel(1)) .> 0)    # expand mask with 1 pixel
    neighbormask = regmask_expanded .& .!regmask
    neighbors = unique(mergeregions[neighbormask])
    return neighbors[neighbors.>0 .& neighbors.<=numreg .& neighbors.!=r]
end

function connectedregions(regions, offshoreregions, numreg)
    connected, connectedoffshore = zeros(Bool, numreg, numreg), zeros(Bool, numreg, numreg)
    mergeregions = [(r > 0) ? r : offshoreregions[i] for (i,r) in enumerate(regions)]
    for r = 1:numreg     
        neighbors = getneighbors(regions, r)
        connected[r, neighbors] .= true

        neighborsoffshore = getneighbors(mergeregions, r)
        connectedoffshore[r, neighborsoffshore] .= .!connected[r, neighborsoffshore] 
    end
    return connected, connectedoffshore
end

# function showmap(; transmission=false, regionnames=true, superregions=false)
# end


# % adjust badly placed regional centers
# if numreg==21
#     latcenter(1) = 60; loncenter(1) = 13;     % NOR
#     latcenter(4) = 52.5; loncenter(4) = -2;   % UK
#     latcenter(5) = 44; loncenter(5) = 19.5;   % MED
#     latcenter(6) = 53; loncenter(6) = 20;     % BAL
#     latcenter(20) = 31; loncenter(20) = 95;   % CH_SW
#     lat(5) = 43; lon(5) = 17.5;   % MED
# elseif numreg==10
#     latcenter(1) = 60.3; loncenter(1) = 15;     % NOR
#     latcenter(5) = 52.5; loncenter(5) = -2;   % UK
#     latcenter(7) = 56;     % BAL
#     latcenter(2) = 42.5; loncenter(2) = 13;     % IT
# end

# if showsuperregions
#     % find borders between superregions
#     superregions = mergeregions;
#     superregions(superregions>0 & superregions<=8) = 1;
#     superregions(superregions>8 & superregions<=15) = 2;
#     superregions(superregions>15 & superregions<=21) = 3;
#     superregions(superregions>3) = 0;

#     superregionborders = zeros(size(regions), 'uint8');
#     for r = 1:3
#         [pixeldistance,ndx] = bwdist(superregions==r);
#         superregionborders(pixeldistance>0 & pixeldistance<=2 & superregions>0) = 1;
#         %disk = fspecial('disk', 1);
#         %expanded = filter2(disk, superregions==r) > 0.1;
#         %superregionborders(expanded & superregions>0 & superregions~=r) = 1;
#     end
#     superregionborders(655:700,:) = 0;  
#     regions(superregionborders > 0) = numreg + 2;
# end

# %newmap; mapcolors(numreg+1); meshm(mergeregions,R);
# %newmap; mapcolors(numreg+1); meshm(regions,R);
# %geoshow(lat, lon, 'DisplayType', 'point', 'MarkerEdgeColor', 'black')
# %distances,connected,connectedoffshore,regionlist'


# f1 = figure;
# if numreg==21
#     axesm('eqdconicstd','MapParallels',[50 25],'MapLatLimit',[12 78],'MapLonLimit',[-18 153]);
# elseif numreg==10
#     axesm('eqdconicstd','MapParallels',[60 45],'MapLatLimit',[25 80],'MapLonLimit',[-20 40]);
# else
#     axesm('eqdconicstd');
# end
# axis('off'); set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
# gridm; mlabel('south'); plabel;
# br = brewermap(8, 'Set2'); br(8,:) = []; br = [br; br; br];
# if numreg==10
#     br([8 10 2 6],:) = br([10 8 6 2],:);
# end
#     cmap = [0.3 0.3 0.45; br(1:numreg,:); 0.5 0.5 0.5];
# if showsuperregions
#     cmap = [cmap; 0 0 0]
# end
# colormap(cmap);
# meshm(regions,R);


# %f2 = newmap; mapcolors(numreg+1); meshm(regions,R);
# %axis([-0.219 1.871 0.265 1.31]);
# f2 = figure;
# if numreg==21
#     axesm('mollweid','MapLatLimit',[-15 90],'MapLonLimit',[-60 150]);
# elseif numreg==10
#     axesm('mollweid','MapLatLimit',[25 80],'MapLonLimit',[-20 40]);
# else
#     axesm('mollweid');
# end
# axis('off'); axis tight; set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]); gridm;
# %cmap = [0 0 0.8; distinguishable_colors(numreg+1, [0.7 0.7 0.7; 0 0 0; 0.7 0 0; 0 0.7 0; 0 0 0.7; 0.2 0 0; 0 0.2 0; 0 0 0.2; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5; 0.2 0 0.2; 0 0.2 0.6; 1 1 1; 1 1 0; 1 0.5 0])];
# %cmap(1,:) = [0.6 0.6 0.7]; cmap(end,:) = [0.7 0.73 0.7]; colormap(cmap);
# %br = brewermap(12, 'Paired'); br([1 9],:) = []; br = [br; br; br];
# %br = brewermap(9, 'Set1'); br(9,:) = []; br = [br; br; br];
# br = brewermap(8, 'Set2'); br(8,:) = []; br = [br; br; br];
# if numreg==10
#     br([8 10 2 6],:) = br([10 8 6 2],:);
# end
#     cmap = [0.3 0.3 0.45; br(1:numreg,:); 0.5 0.5 0.5];
# if showsuperregions
#     cmap = [cmap; 0 0 0]
# end
# colormap(cmap);

# meshm(regions,R);
# axis([-0.89 1.30 0.22 1.35]);

# for r1 = 1:numreg
#     for r2 = 1:numreg
#         distances(r1,r2) = 2*pi*6371/360*distance(lat(r1),lon(r1),lat(r2),lon(r2));
#         [lats,lons] = track2('gc',lat(r1),lon(r1),lat(r2),lon(r2));
#         if showtransmission
#             if connected(r1,r2)
#                 set(0,'CurrentFigure',f1); plotm(lats,lons,'k');
#                 set(0,'CurrentFigure',f2); plotm(lats,lons,'k');
#             elseif connectedoffshore(r1,r2)
#                 set(0,'CurrentFigure',f1); plotm(lats,lons,'w');
#                 set(0,'CurrentFigure',f2); plotm(lats,lons,'w');
#             end
#         end
#     end
#     if showregionnames
#         set(0,'CurrentFigure',f1);  %'Interpreter','none' gets rid of underscores
#         textm(latcenter(r1), loncenter(r1), regionlist{r1},'HorizontalAlignment','center','Interpreter','none');
#         set(0,'CurrentFigure',f2);
#         textm(latcenter(r1), loncenter(r1), regionlist{r1},'HorizontalAlignment','center','Interpreter','none');
#     end
# end

# %save("distances" + filesuffix, 'distances', 'connected', 'connectedoffshore', 'regionlist', 'population', 'demand')
