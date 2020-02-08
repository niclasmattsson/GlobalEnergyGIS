# Example script to 
# (i) train machine learning model (gradient boosting tree regression)
# (ii) read in a file with regions and
# (iii) output synthetic time series for regions (based on (i))

import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
import numpy as np
import sys
import h5py

REGION = sys.argv[1]
# REGION = 'eu10' 
# REGION = 'MENA'

# inputs
INPUT_DIR = 'data/'
OUTPUT_DIRECTORY = 'output/'
INPUT_DIR_TS = INPUT_DIR+'temperature_ranked_regions/'
INPUT_FILE_REGIONS = INPUT_DIR+'julia/regiondata_'+REGION+'.h5'
INPUT_FILE_COUNTRY_CODES_AND_NAMES = INPUT_DIR+'codes_and_country_names_MOD.xlsx'
INPUT_FILE_GDP = INPUT_DIR+'BNP.xlsx'
INPUT_FILE_POPULATION_AGGREGATE = INPUT_DIR+'POP.xlsx'
INPUT_FILE_ML_DATA = INPUT_DIR+"df_model_features.csv"
INPUT_FILE_GPW_GLOBAL_POPULATION_WORLD = INPUT_DIR + 'gpw_v4_population_count_rev10_2010_30_min.shp'

h5data = h5py.File(INPUT_FILE_REGIONS, 'r')

regions_regionlist = h5data['regionlist']
regions_demand = h5data['regionaldemand']
regions_pop = h5data['regionalpop']
regions_gdp = h5data['regionalgdp']

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

column_names = ['hour', 'ranked_month', 'weekend01', 'temperature_top3_mean', 'ppp2015', 'temperature1_qlow', 'temperature1_qhigh', 'temperature1_mean']
# column_names = ['hour', 'temperature1_month_mean', 'weekend01', 'temperature1', 'ppp2015', 'temperature1_qlow', 'temperature1_qhigh', 'temperature1_mean']

full_model = GradientBoostingRegressor(n_estimators=100, subsample=0.75, verbose=1, n_iter_no_change=3, validation_fraction=0.4)
X = df_all_features[column_names]
y = df_all_features['normalized_demand_per_capita']
full_model.fit(X,y)
print("Fitted GradientBoostingRegressor model")



for region_code in range(len(regions_regionlist)):
    print(region_code)
    region_name = regions_regionlist[region_code]

    region_id = str(region_code+1) + "_" + region_name
    region_gdp_per_cap = round(regions_gdp[region_code] / regions_pop[region_code] / 1000, 2)  # Thousands of US$ per capita
            
    fn_population_ranked_time_series = INPUT_DIR + 'julia/tempdata_' + REGION + '_' + region_id + '.csv'
    print(fn_population_ranked_time_series)
    df_temp = pd.read_csv(fn_population_ranked_time_series, parse_dates=True)
    df_temp.set_index(pd.DatetimeIndex(df_temp['time']),inplace=True)
    df_temp.drop('time',axis=1,inplace=True)
    df_temp['year'] = df_temp.index.year

    df_country_features = df_temp.copy(deep=True)
    df_country_features['month'] = df_country_features.index.month
    group_country_months = df_country_features.groupby('month')
    months_temperature_means = group_country_months.aggregate({'temp1':np.mean})
    months_temperature_order = np.argsort(months_temperature_means.values.T).T.flatten()+1
    ranks = np.arange(13,24+1)
    dict_month_to_rank = dict(zip(months_temperature_order,ranks))
    dict_month_to_temperature1_mean = months_temperature_means.to_dict()['temp1']
    df_country_features['ranked_month'] = df_country_features['month'].apply(lambda x: dict_month_to_rank[x])
    df_country_features['temperature1_month_mean'] = df_country_features['month'].apply(lambda x: dict_month_to_temperature1_mean[x])

    df_country_features['temperature1'] = df_temp.temp1
    df_country_features['temperature_top3_mean'] = (df_temp.temp1 + df_temp.temp2 + df_temp.temp3)/3
    # df_country_features.rename({'temperature':'temperature_top3_mean'},axis=1, inplace=True)
    df_country_features['hour'] = df_country_features.index.hour
    df_country_features['weekend'] = df_country_features.index.weekday>4
    df_country_features['weekend01'] = df_country_features['weekend'].apply(lambda x: 1 if x else 0)
    df_country_features['temperature1_qlow'] = df_temp.temp1.quantile(0.05)
    df_country_features['temperature1_qhigh'] = df_temp.temp1.quantile(0.95)
    df_country_features['temperature1_mean'] = df_temp.temp1.mean()
    #df_country_features['region_gdp_per_cap'] = region_gdp_per_cap
    df_country_features['ppp2015'] = region_gdp_per_cap

    #plt.title(region_id  + " " + str(region_gdp_per_cap))
    # fn_region = OUTPUT_DIRECTORY+'region_'+region_id+'.png'
    #plt.savefig(fn_region)
    #plt.show()
    
    X = df_country_features[column_names].values
    predicted = full_model.predict(X)
    #plt.plot(predicted)
    #plt.title(region_id)
    #plt.savefig(OUTPUT_DIRECTORY+'prediction_'+region_id+'.png')

    # Then scale back
    predicted_region_demand_twh = predicted*regions_demand[region_code]/(365*24)

    # Normalize if slight deviation from 0
    predicted_region_demand_twh = predicted_region_demand_twh*(regions_demand[region_code]/np.sum(predicted_region_demand_twh))

    predicted_region_demand_wh_per_cap = 10**12*predicted_region_demand_twh/regions_pop[region_code]
    setup = np.array([predicted_region_demand_twh,predicted_region_demand_wh_per_cap]).T
    setup_columns = ['predicted_region_demand_twh','predicted_region_demand_wh_per_cap']
    df_prediction = pd.DataFrame(setup, columns=setup_columns)
    timeindex = pd.date_range('2050-01-01', periods=8760, freq='H')
    df_prediction.set_index(timeindex,inplace=True)
    fn_prediction = OUTPUT_DIRECTORY+'synthetic2050_' + REGION + '_'+region_id+'.csv'
    df_prediction.to_csv(fn_prediction)
