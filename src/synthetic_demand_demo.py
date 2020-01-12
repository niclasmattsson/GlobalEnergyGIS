# Example script to 
# (i) train machine learning model (gradient boosting tree regression)
# (ii) read in a file with regions and
# (iii) output synthetic time series for regions (based on (i))

from scipy.io import loadmat
import pandas as pd
from pylab import rcParams
rcParams['figure.figsize'] = 30,15
rcParams['axes.titlesize'] = 20
from sklearn.ensemble import GradientBoostingRegressor
import numpy as np
import geopandas as gpd
import sys
from scipy.spatial import ConvexHull
from osgeo import osr
from osgeo import ogr
from shapely.geometry import Polygon
#import matplotlib.pylab as plt

REGION = sys.argv[1]
# REGION = 'eu10' 
# REGION = 'MENA'

# inputs
INPUT_DIR = 'data/'
OUTPUT_DIRECTORY = 'output/'
INPUT_DIR_TS = INPUT_DIR+'temperature_ranked_regions/'
INPUT_FILE_REGIONS = INPUT_DIR+REGION+'.mat'
INPUT_FILE_COUNTRY_CODES_AND_NAMES = INPUT_DIR+'codes_and_country_names_MOD.xlsx'
INPUT_FILE_GDP = INPUT_DIR+'BNP.xlsx'
INPUT_FILE_POPULATION_AGGREGATE = INPUT_DIR+'POP.xlsx'
INPUT_FILE_ML_DATA = INPUT_DIR+"df_model_features.csv"
INPUT_FILE_GPW_GLOBAL_POPULATION_WORLD = INPUT_DIR + 'gpw_v4_population_count_rev10_2010_30_min.shp'
# workflow file names
FILENAME_MERRA2_GRID = OUTPUT_DIRECTORY+'merra2.shp'
FILENAME_REGIONS_SHP = OUTPUT_DIRECTORY+'regions_'+REGION+'.shp'

datadict = dict()
loadmat(INPUT_FILE_REGIONS,datadict);

regions_regions = datadict['regions']
regions_pop = datadict['pop']
regions_popdens = datadict['popdens']
regions_gdp = datadict['gdp']
regions_regionlist = datadict['regionlist']
regions_demand = datadict['demand']

df_codes_and_country_names = pd.read_excel(INPUT_FILE_COUNTRY_CODES_AND_NAMES)
dict_WB_country_name_to_country_code = df_codes_and_country_names.set_index('Countryname').to_dict()['ISO']
dict_WB_country_code_to_country_name = df_codes_and_country_names.set_index('ISO').to_dict()['Countryname']

gdp_adjusted = pd.read_excel(INPUT_FILE_GDP)
gdp_adjusted.set_index('Region',inplace=True)
gdp_adjusted['2015'] = (gdp_adjusted[2010] + gdp_adjusted[2020])/2
dict_country_code_to_gdp = gdp_adjusted.to_dict()['2015']

pop_aggregate = pd.read_excel(INPUT_FILE_POPULATION_AGGREGATE, skipfooter=1)
pop_aggregate = pop_aggregate[pop_aggregate.Model == 'IIASA GDP']
pop_aggregate.set_index('Region',inplace=True)
pop_aggregate['2015'] = (pop_aggregate[2010] + pop_aggregate[2020])/2
dict_country_code_to_pop_aggregate = pop_aggregate.to_dict()['2015']

df_all_features = pd.read_csv(INPUT_FILE_ML_DATA, index_col='time', parse_dates=True)
included_countries = sorted(df_all_features.country.unique())
for country_name in included_countries:
    if not country_name in dict_WB_country_name_to_country_code.keys():
        print(country_name, "missing!")

df_all_features['ppp2015'] = df_all_features.country.apply(lambda x: dict_country_code_to_gdp[dict_WB_country_name_to_country_code[x]]/dict_country_code_to_pop_aggregate[dict_WB_country_name_to_country_code[x]])

# normalize features
gi = 0
setup_df_all_features_normalized = []
for country in included_countries:
    print(country)
    gi += 1
    country_data = df_all_features[df_all_features.country==country].copy(deep=True)
    country_data['normalized_demand_per_capita'] = country_data.demand_per_capita/country_data.demand_per_capita.mean()
    setup_df_all_features_normalized.append(country_data)
df_all_features = pd.concat(setup_df_all_features_normalized)

column_names = ['hour', 'ranked_month', 'weekend01','temperature_top3_mean',
                'ppp2015','temperature1_qlow', 'temperature1_qhigh','temperature1_mean']

full_model = GradientBoostingRegressor(n_estimators=100, subsample=0.75, verbose=1, n_iter_no_change=3, validation_fraction=0.2)
X = df_all_features[column_names]
y = df_all_features['normalized_demand_per_capita']
full_model.fit(X,y)
print("Fitted GradientBoostingRegressor model")

# Setup MERRA-2 grid
merra_grid = np.zeros((361,576),dtype='2float32')
for i in range(361):
    for j in range(576):
        merra_grid[i,j] = [-90 + (1.0/2)*i,-180+(5.0/8)*j]
# Flip coordinates so that matrix starts from north-west
merra_grid_normal = np.flip(merra_grid,axis=0)
merra_grid_normal[0,0]
# Shape file
target = osr.SpatialReference()
target.ImportFromEPSG(4326)
merra_fn = FILENAME_MERRA2_GRID
merra_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(merra_fn)
sr = target
shp_lyr = merra_ds.CreateLayer('merra2',sr,ogr.wkbPoint)
shp_lyr.CreateField(ogr.FieldDefn('coord',ogr.OFTString))
for i in range(361):
    for j in range(576):
        (lat,lon) = merra_grid_normal[i,j]
        
        p = ogr.Geometry(ogr.wkbPoint)
        p.AddPoint(float(lon),float(lat))
        shp_row = ogr.Feature(shp_lyr.GetLayerDefn())
        shp_row.SetGeometry(p)
        shp_row.SetField('coord',"lat:"+str(lat)+" lon:"+str(lon))
        shp_lyr.CreateFeature(shp_row)
del merra_ds
gdf_merra2 = gpd.read_file(FILENAME_MERRA2_GRID)

# Regions as collections of points, cast to MERRA-2 GRID
target = osr.SpatialReference()
target.ImportFromEPSG(4326)
regions_fn = FILENAME_REGIONS_SHP
regions_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(regions_fn)
sr = target
shp_lyr = regions_ds.CreateLayer('cell',sr,ogr.wkbPoint)
shp_lyr.CreateField(ogr.FieldDefn('region',ogr.OFTString))
y_delta = 180/2160
x_delta = 360/4320
region = 0
for yi in np.arange(0,1080):
    y_start = 90-yi*y_delta
    y_end = y_start-y_delta
    print(".",end="")
    #print("[",y_start,",",y_end,"]\t")    
    for xi in np.arange(0,regions_regions.shape[1]):
        x_start = -180+xi*x_delta
        x_end = x_start + x_delta
        p = ogr.Geometry(ogr.wkbPoint)
        p.AddPoint(float(x_start+x_delta/2),float(y_start+y_delta/2))
        shp_row = ogr.Feature(shp_lyr.GetLayerDefn())
        shp_row.SetGeometry(p)
        region = regions_regions[yi,xi]
        shp_row.SetField('region',str(region))
        shp_lyr.CreateFeature(shp_row)
del regions_ds

population_global = gpd.read_file(INPUT_FILE_GPW_GLOBAL_POPULATION_WORLD)
gpd_regions = gpd.read_file(FILENAME_REGIONS_SHP)
gpd_regions['region'] = gpd_regions['region'].apply(int)

gdf_merra2 = gpd.read_file(FILENAME_MERRA2_GRID)
for region_code in sorted(gpd_regions.region.unique()):
    print(region_code)
    # assume: region with first code should be skipped (encoded northern hemisphere)
    # assume: region with highest code should be skipped (encoded northern hemisphere)
    if region_code>0 and region_code<len(gpd_regions.region.unique())-1:
        print("Analyzing region", region_code)
        region_name = regions_regionlist[region_code-1][0][0]
        region = gpd_regions.query('region==@region_code')
        #ax = region.plot()

        points = np.array(list(zip(region.geometry.x, region.geometry.y)))
        hull = ConvexHull(points)
        hullpoints = []
        for vertex in hull.vertices:
            point = points[vertex]
            hullpoints.append(point)
        polyhull = Polygon((np.flip(hullpoints,axis=0)))

        gs_polyhull = gpd.GeoSeries(polyhull) 
        gdf_polyhull = gpd.GeoDataFrame(gs_polyhull)
        gdf_polyhull = gdf_polyhull.rename(columns={0:'geometry'}).set_geometry('geometry')
        #gdf_polyhull.boundary.plot(ax=ax,color='green')

        poly = gdf_polyhull.head(1).geometry.values[0]
        within = gdf_merra2.within(poly)
        subset = gdf_merra2[within]
        #subset.plot(ax=ax,color='red',alpha=0.3)

        bounds = gdf_polyhull.geometry.bounds
        minx = bounds.minx.values[0]
        miny = bounds.miny.values[0]
        maxx = bounds.maxx.values[0]
        maxy = bounds.maxy.values[0]

        population_overlap = population_global[population_global.within(poly)]
        # if no population entry is entirely within polygon, try boundary overlap
        if len(population_overlap) == 0:
            population_overlap = population_global[population_global.overlaps(poly)]            
        population_overlap_candidates = population_overlap.sort_values('DN',ascending=False).head(20)
        population_overlap_selected_indices = []
        d = dict()
        for k in np.arange(1,20):
            polyk = population_overlap_candidates.head(k).tail(1).geometry.values[0]
            within = region.within(polyk)
            if np.sum(within)>0:
                population_overlap_selected_indices.append(k-1) # iloc starts at 0, but head at 1               
            if len(population_overlap_selected_indices)==5:
                break
        # If we did not find the required number of population entries, pad with largest population
        # (Need to evaluate how this works for very small regions)
        if len(population_overlap_candidates)<5:
            diff = 5 - len(population_overlap_candidates)
            for di in np.arange(diff):
                idx = diff - 1 + di
                print(idx)
                population_overlap_selected_indices[idx] = population_overlap_selected_indices[0]
                        
        population_overlap_selected = population_overlap_candidates.iloc[population_overlap_selected_indices]
        #population_overlap_selected.plot(column='DN',ax=ax,cmap='Wistia')
        #plt.show()

        region_id = str(region_code) + "_" + region_name
        for k in np.arange(5):
            ag = population_overlap_selected.head(k+1).tail(1)
            lat = ag.geometry.centroid.y.values[0]
            lon = ag.geometry.centroid.x.values[0]
            coordinates = str(round(lat,6))+","+str(round(lon,6))
            d[region_id+"_"+str(k)] = coordinates                        
        #fn_population_ranked_time_series = OUTPUT_DIRECTORY+'/regions/population_ranked_time_series_'+region_id+".npy"
        #np.save(fn_population_ranked_time_series, d)
        
        region_where = np.where(regions_regions==region_code)
        region_gdp = regions_gdp[region_where]
        region_gdp_tot = np.sum(region_gdp)
        region_pop = regions_pop[region_where]
        region_pop_tot = np.sum(region_pop)
        region_gdp_per_cap = region_gdp_tot / region_pop_tot
        region_gdp_per_cap = region_gdp_per_cap / 1000 # Thousands of US $
        region_gdp_per_cap = round(region_gdp_per_cap,2)
                
        k = 1
        fn_population_ranked_time_series = INPUT_DIR_TS + region_id + '_' + str(k) + '.csv'
        print(fn_population_ranked_time_series)
        df_temp = pd.read_csv(fn_population_ranked_time_series, parse_dates=True)
        df_temp.set_index(pd.DatetimeIndex(df_temp['time']),inplace=True)
        df_temp.drop('time',axis=1,inplace=True)
        df_temp['year'] = df_temp.index.year
        df_temp = df_temp[df_temp.year==2015]                            # Use year 2015
        df_temp_top1 = df_temp
        
        k = 2
        fn_population_ranked_time_series = INPUT_DIR_TS + region_id + '_' + str(k) + '.csv'
        print(fn_population_ranked_time_series)
        df_temp = pd.read_csv(fn_population_ranked_time_series, parse_dates=True)
        df_temp.set_index(pd.DatetimeIndex(df_temp['time']),inplace=True)
        df_temp.drop('time',axis=1,inplace=True)
        df_temp['year'] = df_temp.index.year
        df_temp = df_temp[df_temp.year==2015]                            # Use year 2015
        df_temp_top2 = df_temp

        k = 3
        fn_population_ranked_time_series = INPUT_DIR_TS + region_id + '_' + str(k) + '.csv'
        print(fn_population_ranked_time_series)
        df_temp = pd.read_csv(fn_population_ranked_time_series, parse_dates=True)
        df_temp.set_index(pd.DatetimeIndex(df_temp['time']),inplace=True)
        df_temp.drop('time',axis=1,inplace=True)
        df_temp['year'] = df_temp.index.year
        df_temp = df_temp[df_temp.year==2015]                            # Use year 2015
        df_temp_top3 = df_temp

        df_country_features = df_temp_top1.copy(deep=True)
        df_country_features.rename({'temperature':'temperature1'},axis=1, inplace=True)
        df_country_features['month'] = df_country_features.index.month
        group_country_months = df_country_features.groupby('month')
        months_temperature_means = group_country_months.aggregate({'temperature1':np.mean})
        months_temperature_order = np.argsort(months_temperature_means.values.T).T.flatten()+1
        ranks = np.arange(13,24+1)
        dict_month_to_rank = dict(zip(months_temperature_order,ranks))
        dict_month_to_temperature1_mean = months_temperature_means.to_dict()['temperature1']
        df_country_features['ranked_month'] = df_country_features['month'].apply(lambda x: dict_month_to_rank[x])
        
        df_country_features['temperature'] = (df_temp_top1.temperature + df_temp_top2.temperature + df_temp_top3.temperature)/3
        df_country_features.rename({'temperature':'temperature_top3_mean'},axis=1, inplace=True)
        df_country_features['hour'] = df_country_features.index.hour
        df_country_features['weekend'] = df_country_features.index.weekday>4
        df_country_features['weekend01'] = df_country_features['weekend'].apply(lambda x: 1 if x else 0)
        df_country_features['temperature1_qlow'] = df_temp_top1.temperature.quantile(0.05)
        df_country_features['temperature1_qhigh'] = df_temp_top1.temperature.quantile(0.95)
        df_country_features['temperature1_mean'] = df_temp_top1.temperature.mean()
        #df_country_features['region_gdp_per_cap'] = region_gdp_per_cap
        df_country_features['ppp2015'] = region_gdp_per_cap
        
        #plt.title(region_id  + " " + str(region_gdp_per_cap))
        fn_region = OUTPUT_DIRECTORY+'region_'+region_id+'.png'
        #plt.savefig(fn_region)
        #plt.show()
        
        X = df_country_features[column_names].values
        predicted = full_model.predict(X)
        #plt.plot(predicted)
        #plt.title(region_id)
        #plt.savefig(OUTPUT_DIRECTORY+'prediction_'+region_id+'.png')

        # Then scale back
        predicted_region_demand_twh = predicted*regions_demand[0][region_code-1]/(365*24)

        # Normalize if slight deviation from 0
        predicted_region_demand_twh = predicted_region_demand_twh*(regions_demand[0][region_code-1]/np.sum(predicted_region_demand_twh))

        predicted_region_demand_wh_per_cap = 10**12*predicted_region_demand_twh/region_pop_tot
        setup = np.array([predicted_region_demand_twh,predicted_region_demand_wh_per_cap]).T
        setup_columns = ['predicted_region_demand_twh','predicted_region_demand_wh_per_cap']
        df_prediction = pd.DataFrame(setup, columns=setup_columns)
        timeindex = pd.date_range('2050-01-01', periods=8760, freq='H')
        df_prediction.set_index(timeindex,inplace=True)
        fn_prediction = OUTPUT_DIRECTORY+'synthetic2050_region'+region_id+'.csv'
        df_prediction.to_csv(fn_prediction)
