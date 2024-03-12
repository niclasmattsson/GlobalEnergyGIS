# Run this file with include() after adding packages below to the environment.
# E.g. includet("./GlobalEnergyGIS/src/external_scripts.jl")

using GMT

G = GlobalEnergyGIS

function geotest8()
    pol, nodes = sphtriangulate("@gshhs_c.txt", voronoi=:v, skip=true, nodes=true)

    # Compute distances in km
    Gtt = sphdistance(pol, region=:global360, inc=1, voronoi=true, nodes=nodes, dist_unit=:k)
    t_cpt = makecpt(cmap=:hot, range=(0,3500))

    # Make a basic image plot and overlay contours, Voronoi polygons and coastlines
    grdimage(Gtt, proj=(name=:ortho, center=(-140,30)), figsize=16, xshift=2, yshift=5)
    grdcontour!(Gtt, cont=500, annot=(int=1000, labels=(font=(10,:Helvetica,:white),)),
                range=500, labels=(Line=[0 90 203 -10; 175 60 170 -30; -50 30 220 -5],),
                pen=((annot=true, pen=(0.75,:white)), (contour=true, pen=(0.25,:white))) )

    GMT.plot!(pol, pen=(0.25, :green, :dotted))
    coast!(shore=1, land=:steelblue, area=(0,1,1),
        frame=(annot=30, grid=30, title="Distances from GSHHG crude coastlines"), show=true)
end

function geotest9()
    table_5 = gmtread("@Table_5_11.txt")    # The data used in this example
    T = gmtinfo(table_5, nearest_multiple=(dz=25, col=2))
    makecpt(color=:jet, range=T.text[1][3:end])  # Make it also the current cmap

    subplot(grid=(2,2), limits=(0,6.5,-0.2,6.5), col_axes=(bott=true,), row_axes=(left=true,),
            figsize=8, margins=0.1, panel_size=(8,0), tite="Delaunay Triangulation")
        # First draw network and label the nodes
        net_xy = triangulate(table_5, M=true)
        GMT.plot(net_xy, lw=:thinner)
        GMT.plot(table_5, marker=:circle, ms=0.3, fill=:white, MarkerLine=:thinnest)
        GMT.text(table_5, font=6, rec_number=0)

        # Then draw network and print the node values
        GMT.plot(net_xy, lw=:thinner, panel=(1,2))
        GMT.plot(table_5, marker=:circle, ms=0.08, fill=:black)
        GMT.text(table_5, zvalues=true, font=6, justify=:LM, fill=:white, pen="", clearance="1p", offset=("6p",0), noclip=true)

        # Finally color the topography
        GMT.contour(table_5, pen=:thin, mesh=(:thinnest,:dashed), labels=(dist=2.5,), panel=(2,1))
        GMT.contour(table_5, colorize=true, panel=(2,2))
    subplot("show")
end

function ehub500()
    buses = G.CSV.read(G.in_datafolder("Bus_data_EHUB500 - buses_original.csv"), G.DataFrame)
    unique!(buses, ["x-coordinate", "y-coordinate"])
    xy = Matrix(buses[:, ["x-coordinate", "y-coordinate"]])

    # bbox = collect(Iterators.flatten(extrema(uxy, dims=1)))
    bbox = [4, 32, 54.5, 71.5]
    ds = GMTdataset(xy)
    ds.bbox = bbox
    ds.ds_bbox = bbox
    
    pol, nodes = sphtriangulate(ds, voronoi=:v, skip=true, nodes=true)

    t_cpt = makecpt(cmap=:categorical, range=(0,size(xy,1),1), wrap=:w)
    GMT.plot(pol, proj=:moll, region=bbox, close=true, cmap=t_cpt, alpha=65, wrap=:w)
    GMT.scatter!(nodes, fill=:red, markersize="2p")
    coast!(land=nothing, DCW=(country="SE,NO,DK,FI", pen=(0.25,:black)), frame=(annot=:auto, ticks=:auto, grid=:auto), show=true)

    gmtwrite("ehub500_id.shp", pol)
    # gmtwrite("ehub500.geojson", pol)
    return pol
end

function ehub500_test()
    buses = G.CSV.read(G.in_datafolder("Bus_data_EHUB500 - buses_original.csv"), G.DataFrame)
    unique!(buses, ["x-coordinate", "y-coordinate"])
    xy = Matrix(buses[:, ["x-coordinate", "y-coordinate", "bus_id"]])

    # bbox = collect(Iterators.flatten(extrema(uxy, dims=1)))
    bbox = [4, 32, 54.5, 71.5, 5500, 90000]
    # ds = GMTdataset(xy)
    # ds.bbox = bbox
    # ds.ds_bbox = bbox
    ds = mat2ds(xy)
    
    pol, nodes = sphtriangulate(ds, voronoi=:v, skip=true, nodes="nodes.txt", verbose=true)

    gmtwrite("ehub500_id.shp", pol)
    # gmtwrite("test.geojson", pol)
    pol
end

function readshape()
    df = G.GDF.read("ehub500_id.shp")
end
