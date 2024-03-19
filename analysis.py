import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
import seaborn as sns
import os

# LM DATA FORMATTING NOTEBOOK
def read_data(rate_df, methane_value, filename):
    
    ''' 
    Reads in hourly measurement data in csv format and makes minor formatting and unit changes. Original units are micrograms m-3.
    
    Args:
        rate_df: dataframe of rate constants for each OH-species reaction and the molecular mass for each measured species
        methane_value (int): the constant hourly calculated methane concentration (unit: micrograms m-3)
        filename (str): the name of the CSV file, all of which are contained in a folder labelled DEFRA
        
    Returns:
        a pandas dataframe with the hourly measured data that has been cleaned. All units converted from micrograms m-3 to molecules cm-3
    
    '''

    # specifying the file path: contained as a .feather within the DEFRA folder
    filepath = 'DEFRA/' + str(filename) + '.csv'
    # reading in the data and replacing with np.nan
    df = (pd.read_csv(filepath, skiprows=10)).replace(['No data', 'No Data'],np.nan)
        
    # COLUMN FORMATTING
    # removing all columns labelled Status
    df = df[df.columns.drop(list(df.filter(regex='Status')))]
  
    # keeping only universal columns contained within the rate_df
    gas_list = list(rate_df.index)
    df2 = df[[col for col in df.columns if col in gas_list]]
    
    # adding the methane concentration column
    df2['methane'] = [methane_value]*len(df2)
    
    # DATE TIME FORMATTING
    # replacing the 24:00:00 with date time recognised 00:00:00
    df = df.replace('24:00:00','00:00:00')

    # combining Date and Time columns into a Date_Time column with correct datetime format
    df2['Date_Time'] = pd.to_datetime(df['Date']+' '+df['Time'], format='%Y-%m-%d %H:%M:%S')
    # setting index to Date Time
    df2 = df2.set_index('Date_Time')
    #replacing all remaining columns to relevant numeric data type
    df2 = df2.apply(pd.to_numeric)
    
    # DATA CLEANING
    # ensuring no negative values are present as these would be outliers
    df2 = df2[df2.min(axis=1) >= 0]
    # ensuring no NaN values are present in the Ozone column as this would prevent OH calculation
    df2 = df2[df2['Ozone'].notna()]
    
    # UNIT CONVERSION: molecules cm-3
    # generating a list of all gas species
    list_gas = list(df2.columns)
 
    for i in range(len(list_gas)): 
        gas = list_gas[i]
        # accessing molecular mass of each species
        mol_mass = rate_df.loc[gas,'Molecular mass']
        # diving each value by their molecular mass
        df2[gas] = df2[gas].divide(mol_mass) # mol cm-3
    
    # multiplying by Avogadros Constant
    df2 = df2.multiply(6.02214076e23) # molecules cm-3
    
    return df2

def calculate_OH(df, rate_df, temp_dependant_rates, rate_type, year):
    '''
    Calculates the OH radical concentration by assuming the steady state approximation.
    
    Args:
        df: dataframe containing hourly measured species concentration data (units: molecules cm-3)
        rate_df: dataframe of rate constants for each OH-species reaction and the molecular mass for each measured species
        temp_dependant_rates: dataframe containing the calculated rate constants for the OH-NO reaction over 200-300K
        rate_type (str): which selects which calculated OH-NO rate constant to use
            - literature k (JPL literature)
            - all k (all data)
            - float f k (Fcent unfixed)
            - float m k (m unfixed)
            - 298 - 258 k (data only over the 298-258 range)
        year (int): the year the DEFRA data was obtained from
    
    returns:
        Creates a folder for each year and rate constant type in the notebook working directory, format:
            'rate_type-year'
        Writes out the final dataframe as a .feather file. File naming is done as:
            '/London,Marylebone-year-temperature-rate_type.feather'
    
    '''
    
    # Raising an error if incorrect inputs used (rate_type and year)
    if rate_type not in ['literature k','all k','float f k', 'float m k','298-258 k']:
        raise Exception('Incorrect input for rate_type used, please input one of: \n - literature k, all k, float f k, float m k, 298-258')
        
    if year not in [2005, 2010, 2020]:
        raise Exception('Incorrect input for year, please use one of: \n - 2005, 2010, 2020')
    
    # REMOVING MISSING GASES
    # list of all gases in the rate_df that aren't present in the DEFRA data        
    missing_gases = list(set(list(rate_df.index)).difference(df.columns))
    # removing the missing gases from the rate_df
    updated_rate_df = (rate_df[rate_df.index.isin(missing_gases)==False]).reset_index()
    
    # list of all gases
    list_gas = list(updated_rate_df['Gas'])
    # specifying the file path for the chosen rate constant
    file_path = str(rate_type)+'-'+str(year)
    #creating directory of filepath and checking if it already exists
    if os.path.isdir(file_path) == True:
        pass
    if os.path.isdir(file_path) ==False:
        os.makedirs(file_path)
    
    # CALCULATING OH
    # looping over the temperature to calculate OH for each rate constant
    for i in range(len(temp_dependant_rates)):
        # empty list to append the OH radical concentration
        OH_conc = []
        # creating a copy of the original dataframe to work on
        df_instance = df
        # accessing current loop temperature
        temperature = temp_dependant_rates['Temperature (K)'][i]
        # selecting OH-NO rate constant indicated in function initialisation for loop temperature
        chosen_rate_value = temp_dependant_rates.loc[i,rate_type]

        # replacing the old NO rate constant with the selected one
        updated_rate_df.loc[updated_rate_df.Gas == 'Nitric oxide', 'Rate constant'] = chosen_rate_value

        # calculating rate x concentration for each gas
        for j in range(len(list_gas)):
            # accessing each gas
            gas = list_gas[j]
            # retrieving the related rate constant from rate_df
            kgas = updated_rate_df['Rate constant'][j] # cm3 molecule-1 s-1

            # retrieving the gas concentration for whole time period
            gas_conc = df_instance[gas].tolist() #molecules cm-3

            # empty list to append the concentration*rate
            list_kgas = []

            # appending concentration*k
            for z in range(len(gas_conc)):
                gas_mult = gas_conc[z]*kgas
                list_kgas.append(gas_mult) # s-1

            # adding a new column for each k*gas species
            column_header = 'k_'+str(gas)
            df_instance[column_header] = list_kgas

        # SUMMATION OF ALL SINK GASES
        # accessing only the k_sink columns
        dfsinks = (df_instance[df_instance.columns[pd.Series(df_instance.columns).str.startswith('k_')]]).drop('k_Ozone',axis=1)
        # calculating the sum of the sinks for each time point
        ksinks = list(dfsinks.sum(axis='columns'))
        # adding as a separate column in the dataframe
        df_instance['sum_destructive'] = ksinks

        # CALCULATING THE OH CONCENTRATION
        # appending ozone and sink data to a list
        kozone = list(df_instance['k_Ozone'])
        kdestructive = list(df_instance['sum_destructive'])

        # dividing the ozone concentration by the total sink concentration
        for x in range(len(kozone)):
            OH = ((kozone[x]))/kdestructive[x]
            OH_conc.append(OH)
        # adding the OH concentration as a separate column
        df_instance['OH concentration'] = OH_conc
        
        # WRITING OUT
        # resetting index for writing to .feather
        df_instance = df_instance.reset_index()
        
        # outputting the dataframe to a feather format
        df_instance.to_feather(str(file_path) + '/London,Marylebone' + '-'+ str(year) + '-' + str(temperature) + '-' + str(rate_type) + '.feather')

# ANALYSIS NOTEBOOK
def read_data_analysis(rate_type, year, temp_dependant_rates):
    '''
    Format the data into a dataframe more suitable for analysis. The total non-methane VOC concentration, NOx and NMVOC/NOx ratio was calculated for each hour. Any rows where the NMVOCs equalled 0 were dropped.
    
    Args:
        rate_type (str): chosen rate type
            - literature k (JPL literature)
            - all k (all data)
            - float f k (Fcent unfixed)
            - float m k (m unfixed)
            - 298 - 258 k (data only over the 298-258 range)
        year (int): the year the DEFRA data was obtained from
        temp_dependant_rates: dataframe containing the calculated rate constants for the OH-NO reaction over 200-300K
    
    returns:
        A dataframe with the columns: 
                'Date_Time','Nitric oxide', 'OH concentration','NMVOC/NOx ratio','MVOC sum','NMVOC sum', 'Nitrogen dioxide',
                'Ozone','methane','NOx','Carbon monoxide','Temperature'
 
    '''
 
    # Raising an error if incorrect inputs used (rate_type and year)
    if rate_type not in ['literature k','all k','float f k', 'float m k','298-258 k']:
        raise ValueError('Incorrect input for rate_type used, please input one of: \n - literature k, all k, float f k, float m k, 298-258 k')
        
    if year not in [2005, 2010, 2020]:
        raise ValueError('Incorrect input for year, please use one of: \n - 2005, 2010, 2020')
    
    # READING IN ALL FILES FOR EACH RATE CONSTANT
    # defining the folder filepath for accessing the data
    filepath = str(rate_type)+'-'+str(year)+'/'
    # list to append each file to (w/ .feather)
    file_feather_list = []
    # list to append each file to (w/out .feather)
    file_list = []
    # list to append the temperature to
    temp_list = []
    # list to append each temperature df to in order to create a single combined df
    df_list = []
    
    # looping through every temperature to access all files of a rate constant
    for i in range(len(temp_dependant_rates)):
        # accessing temperature
        temperature = temp_dependant_rates.loc[i,'Temperature (K)']
        # filename w/ .feather
        filename_feather = 'London,Marylebone-'+str(year)+'-'+str(temperature)+'-'+str(rate_type)+'.feather'
        # filename w/out .feather
        filename = 'London,Marylebone-'+str(year)+'-'+str(temperature)+'-'+str(rate_type)
        # appending to relevant list
        file_feather_list.append(filename_feather)
        file_list.append(filename)
        temp_list.append(temperature)
        
    # creating a dictionary of dataframes with the key as the file name
    dataframes = (pd.read_feather(str(filepath)+file_feather) for file_feather in file_feather_list)
    dictionary = dict(zip(file_list, dataframes))

    for j in range(len(file_list)):
        # accessing key (filename)
        name = file_list[j]
        # accessing temperature of particular dataframe
        temperature = int(name[23:26])
        # accessing relevent dataframe from dictionary
        df = dictionary[name]
        # adding a temperature column
        df.loc[:,'Temperature'] = temperature

        # CALCULATING A VOC/NOx ratio for analysis
        # selecting VOC columns
        NMVOC_col_list = ['ethane', 'ethene', 'ethyne', 'propane', 'propene', 'iso-butane',
                           'n-butane', '1-butene', 'trans-2-butene', 'cis-2-butene', 'iso-pentane',
                           'n-pentane', '1,3-butadiene', 'trans-2-pentene', '1-pentene',
                           '2-methylpentane', 'isoprene', 'n-hexane', 'n-heptane', 'iso-octane',
                           'n-octane', 'benzene', 'toluene', 'ethylbenzene', 'm+p-xylene',
                           'o-xylene', '1,2,3-trimethylbenzene', '1,2,4-trimethylbenzene',
                           '1,3,5-trimethylbenzene']
        MVOC_col_list = NMVOC_col_list + ['methane']
        
        # summating VOC columns
        NMVOC_sum = df[NMVOC_col_list].sum(axis=1)
        MVOC_sum = df[MVOC_col_list].sum(axis=1)

        # accessing NO values
        NO_list = list(df['Nitric oxide'])
        NO2_list = list(df['Nitrogen dioxide'])

        NOx_list = []

        for i in range(len(df)):
            NOx_list.append(NO_list[i]+NO2_list[i])

        df['NOx'] = NOx_list

        # adding ratio to the dataframe
        df['NMVOC/NOx ratio'] = [NMVOC_sum[i]/NOx_list[i] for i in range(len(NOx_list))]
        df['NMVOC sum'] = NMVOC_sum
        df['MVOC sum'] = MVOC_sum
        
        # keeping only relevent columns
        df_w_summation = df[['Date_Time','Nitric oxide', 'OH concentration','NMVOC/NOx ratio','MVOC sum','NMVOC sum', 'Nitrogen dioxide',
                             'Ozone','methane','NOx','Carbon monoxide','Temperature']]
        # creating a copy to prevent warning being produced
        df_w_summation = df_w_summation.copy()
        
        # adding a column containing the rate type and year
        df_w_summation.loc[:,'rate type'] = rate_type
        df_w_summation.loc[:,'year'] = year
        
        # appending formatted dataframe to final list
        df_list.append(df_w_summation)
        
        # DATAFRAME COMBINATION
        final_df = pd.concat(df_list)
        # removing any rows were no non-methane VOCs were measured
        final_df = final_df[final_df['NMVOC sum'] != 0]
        
        # resetting the index
        final_df = final_df.set_index('Date_Time')
        final_df_reset = final_df.reset_index()
    
    return final_df_reset


def average(df, time_period, format_type):
    '''
    Calculates the daily, monthly and yearly average of all species
    
    Args:
        df: pandas dataframe containing all hourly data with OH concentration
            - output from the read_data_analysis function
        time_period (str):
            - 'year','month','day','hour'
            - if format type = continuous then only available for month and day
        format_type (str):
            'grouped': averages grouped by the rate type, time period and temperature. The daily and hourly mean is also grouped by month.
            'continous': maintains the datetime index and finds the hourly or monthly average.
            
    Returns:
        An averaged dataframe depending on the desired averaging technique.

    '''
    # Raising an error if incorrect inputs used
    if format_type not in ['grouped','continuous']:
        raise ValueError('Incorrect input for format_type used, please input either: \n - grouped, continuous')
    if (format_type == 'grouped') & (time_period not in ['year','month','day','hour']):
        raise ValueError('Incorrect input for time_period used, please input one of: \n - year, month, day, hour')
    if (format_type == 'continuous') & (time_period not in ['month','day']):
        raise ValueError('Incorrect input for time_period used, please input one of: \n - month, day')
        
    # ensuring that the date_time column is in datetime format
    df['Date_Time'] = pd.to_datetime(df['Date_Time'], format='%Y-%m-%d %H:%M:%S')
    # setting as index
    df = df.set_index('Date_Time')
    
    # obtaining a list of unique rate constant types
    rate_type_set = set(df['rate type'])
    rate_type_list = list(rate_type_set)
    
    # GROUPED
    if format_type == 'grouped':
        # yearly average
        if time_period == 'year':
            # adding a column containing the year 
            df['year'] = df.index.year
            # calculating the mean grouped by the rate type, year and temperature
                # as_index=False so that each grouping has a separate column
            df = df.groupby(['rate type','year','Temperature'], as_index=False).mean()
        
        # monthly average
        if time_period == 'month':
            # adding a column containing the month
            df['month'] = df.index.month
            # calculating the mean grouped by the rate type, month and temperature
                # as_index=False so that each grouping has a separate column
            df = df.groupby(['rate type','month','Temperature'],as_index=False).mean()
        
        # daily average
        if time_period == 'day':
            # adding columns for the month and day
            df['day'] = df.index.day
            df['month'] = df.index.month
            # calculating the mean grouped by the rate type, day, month and temperature
                # as_index=False so that each grouping has a separate column
            df = df.groupby(['rate type','day','Temperature','month'],as_index=False).mean()
            
        # hour average by month
        if time_period == 'hour':
            # adding month and hour column
            df['month'] = df.index.month
            df['hour'] = df.index.hour
            # calculating the mean grouped by the rate type, hour, month and temperature
                # as_index=False so that each grouping has a separate column
            df = df.groupby(['rate type','hour','Temperature','month'], as_index=False).mean()

    # CONTINOUS
    if format_type == 'continuous':
        if time_period == 'month':
            monthly_avg_list = []

            for i in range(len(rate_type_list)):
                rate_type = rate_type_list[i]
                rate_df = df[df['rate type'] == rate_type]
#                 rate_df = rate_df.set_index('Date_Time')
                monthly_avg = rate_df.resample('M').mean(numeric_only=True)
                monthly_avg['rate type'] = [rate_type]*len(monthly_avg)

                monthly_avg_list.append(monthly_avg)

            df = pd.concat(monthly_avg_list)

        if time_period == 'day':
            daily_avg_list = []

            for i in range(len(rate_type_list)):
                rate_type = rate_type_list[i]
                rate_df = df[df['rate type'] == rate_type]
                daily_avg = rate_df.resample('d').mean(numeric_only=True)
                daily_avg['rate type'] = [rate_type]*len(daily_avg)

                daily_avg_list.append(daily_avg)

            df = pd.concat(daily_avg_list)
    
    df = df.reset_index()

    return df

def OH_percentage_change(df, df_year, time_period, temp_dependant_rates):
    '''
    Calculates the percentage change in OH concentration compared to the literature value for monthly and daily averaged data.

    Args:
        df: averaged dataframe (daily, monthly)
        df_year: dataframe containing full data for that year (non-averaged)
        time_period (str): the period over which the data was averaged over 
            -'day', 'month'
        temp_dependant_rates: dataframe containing the calculated rate constants for the OH-NO reaction over 200-300K

    returned:
        A dataframe containing an additional 'OH percentage change' column is returned. Only the newly determined C.J Pervical rate constants are included.
    
    '''
    
    # Raising an error if incorrect inputs used
    if time_period not in ['day','month']:
        raise ValueError('Incorrect input for time_period used, please input either: \n - day or month')
        
    # empty list to contain updated dataframes
    df_list = []
    
    # creating a dataframe which only contains the literature rate constant data
    literature_df = df[df['rate type'] == 'literature k']
    
    # obtaining a list of unique rate constant types
    rate_type_set = set(df_year['rate type'])
    rate_type_list = list(rate_type_set)

    # looping over each rate type
    for i in range(len(rate_type_list)):
        # empty list to contain all percentage change values for each rate type
        percentage_list = []
        # accessing current loop rate type
        rate_type = rate_type_list[i]
        # filtering dataframe to current loop rate type
        rate_df = df[df['rate type'] == rate_type]
        
        # monthly averaged data
        if time_period == 'month':
            # looping over each month from 1 to 12
            for month in set(df['month']):
                # filtering rate df to current loop month
                month_df = rate_df[rate_df['month'] == month]
                # filtering literature df to current loop month
                lit_month_df = literature_df[literature_df['month'] == month]

                # looping over each temperature (200-300)
                for z in range(len(temp_dependant_rates)):
                    # accessing current temperaure
                    temp = temp_dependant_rates['Temperature (K)'][z]
                    # filtering the month,rate df to contain current loop temperature
                    temp_df = month_df[month_df['Temperature'] == temp]
                    # filtering the month, literature df to contain current loop temperature
                    lit_temp_df = lit_month_df[lit_month_df['Temperature'] == temp]

                    # accessing the literature OH value
                    lit_OH = float(lit_temp_df['OH concentration'])
                    # accessing the rate OH value
                    df_OH = float(temp_df['OH concentration'])

                    # calculating the percentage change from literature
                    percentage_change = ((df_OH-lit_OH)/lit_OH)*100
                    # appending to list
                    percentage_list.append(percentage_change)

            # adding a percentage change column to each rate_df
            # creating a copy to avoid error
            rate_df = rate_df.copy()
            rate_df['OH percentage change'] = percentage_list
            # appending to list of dataframes
            df_list.append(rate_df)
        
        # daily averaged data
        if time_period == 'day':
            # looping over each month from 1 to 12
            for month in set(df['month']):
                # filtering rate df to current loop month
                month_df = rate_df[rate_df['month'] == month]
                # filtering literature df to current loop month
                lit_month_df = literature_df[literature_df['month'] == month]

                # finding unique number of days in that month
                number_days = set(month_df['day'])
                number_days_list = list(number_days)
                
                # looping over each day in that month
                for x in number_days_list:
                    # filtering rate df to current loop day
                    day_df = month_df[month_df['day'] == x]
                    # filtering literature df to current loop day
                    lit_day_df = lit_month_df[lit_month_df['day'] == x]

                    # looping over each temperature (200-300)
                    for z in range(len(temp_dependant_rates)):
                        # accessing current temperaure
                        temp = temp_dependant_rates['Temperature (K)'][z]
                        # filtering the month,rate df to contain current loop temperature
                        temp_df = day_df[day_df['Temperature'] == temp]
                        # filtering the month, literature df to contain current loop temperature
                        lit_temp_df = lit_day_df[lit_day_df['Temperature'] == temp]

                        # accessing the literature OH value
                        lit_OH = float(lit_temp_df['OH concentration'])
                        # accessing the rate OH value
                        df_OH = float(temp_df['OH concentration'])

                        # calculating the percentage change from literature
                        percentage_change = ((df_OH-lit_OH)/lit_OH)*100
                        # appending to a list
                        percentage_list.append(percentage_change)
                        
            # adding a percentage change column to each rate_df
            # creating a copy to avoid error
            rate_df = rate_df.copy()
            rate_df['OH percentage change'] = percentage_list
            # appending to list of dataframes
            df_list.append(rate_df)
 
    # combining into 1 dataframe
    df_combined = pd.concat(df_list)
    # dropping the literature value as percentage change will always be 1
    df_combined = df_combined[df_combined['rate type'] != 'literature k']
    
    return df_combined

def calculate_OH_percentage_mean(df):
    '''
    Calculates the mean OH concentration and percentage change in OH from the literature value. Grouped by rate constant type and temperature.
    
    Args:
        df: dataframe containing species measurements with calculated OH concentration and percentage change
        
    Returns:
        Dataframe containing the mean OH concentration and percentage change by rate constant type and temperature.
    '''
    # empty lists to contain values
    OH_mean_list = []
    OH_perc_mean_list = []
    temp_list = []
    df_rate_list = []
    
    # list of rate constants
    rate_list = ['298-258 k','all k','float f k','float m k']

    # looping over each rate constant type
    for rate in rate_list:
        # filtering to contain only one rate type
        rate_df = df[df['rate type'] == rate]
        # looping over temperatures
        for i in range(200,310,10):
            # filtering to contain one temperature
            temp_df = rate_df[rate_df['Temperature'] == i]
            # calculating mean OH at temperature
            OH_mean = temp_df['OH concentration'].mean()
            #calculating mean % change from literature at temperature
            OH_perc_mean = temp_df['OH percentage change'].mean()
            
            # appending final values to corresponding list
            OH_mean_list.append(OH_mean)
            OH_perc_mean_list.append(OH_perc_mean)
            temp_list.append(i)
            df_rate_list.append(rate)
    
    #creating final dataframe
    df_perc_mean = pd.DataFrame({'rate type':df_rate_list, 'temperature': temp_list, '% change in OH mean':OH_perc_mean_list, 'mean OH': OH_mean_list})
    
    return df_perc_mean

def check_collapse(df, significance, temp_dependant_rates):
    '''
    Checks if a collapse in the OH radicals has occured, based on an hourly reduction in their concentration to the chosen significance.
    
    Args:
        df: dataframe of measurements including the OH concentration
        significance (float): the degree of significance to test to
        temp_dependant_rates: dataframe containing the calculated rate constants for the OH-NO reaction over 200-300K
        
    Returns:
        final_df: dataframe containing data from hour of collapse
        pre_final_df: dataframe containing data from hour prior to collapse
        post_final_df: dataframe containing data from hour after collapse
            - if hour before/after not available then skipped
    '''
    # temperature list
    temp_list = []
    # list to contain dataframes
    df_list = []
    # list to contain dataframe with the NO values which caused collapse
    pre_df_list = []
    post_df_list = []
    # list to hold index values
    index = []
    # list to hold the date at which collapse occured 
    collapse_dates_post = []
    collapse_dates_pre = []
    
    # looping through each temperature
    for i in range(len(temp_dependant_rates)):
        # list to append the difference between current and previous time
        difference = []
        
        # accessing the current temperature and filtering dataframe to the temperature
        temp = temp_dependant_rates['Temperature (K)'][i]
        temp_list.append(temp)
        temp_df = df[df['Temperature'] == temp]
        
        # accessing OH values
        OH = list(temp_df['OH concentration'])

        # calculating the difference between current and prior OH
        for i in range(len(temp_df)):
            # current OH
            OH1 = OH[i]
            # prior OH
            OH2 = OH[i-1]
            # appending the difference
            difference.append(OH2-OH1)

        # adding the difference as a column to the current temp dataframe
        temp_df = temp_df.copy()
        temp_df['difference'] = difference
        temp_df = temp_df.reset_index()

        # filtering to only include differences greater than or equal to chosen significance
        collapsed_df = temp_df[temp_df['difference'] >= significance]
        collapsed_df = collapsed_df.copy()
        # checking in datetime
        collapsed_df['Date_Time'] = pd.to_datetime(collapsed_df['Date_Time'])
        
        # finding the hour before and after collapse
        pre_coll_df_date = list(collapsed_df['Date_Time'] - pd.Timedelta(hours=1))
        post_coll_df_date = list(collapsed_df['Date_Time'] + pd.Timedelta(hours=1))
        # filtering to before and after date times
        pre_coll_df = temp_df[temp_df['Date_Time'].isin(pre_coll_df_date)]
        post_coll_df = temp_df[temp_df['Date_Time'].isin(post_coll_df_date)]
        
        pre_df_list.append(pre_coll_df)
        post_df_list.append(post_coll_df)
        
        collapsed_df = collapsed_df.set_index('Date_Time')
        
        # adding a month, day and hour column
        collapsed_df['month'] = collapsed_df.index.month
        collapsed_df['day'] = collapsed_df.index.day
        collapsed_df['hour'] = collapsed_df.index.hour
        
        collapsed_df = collapsed_df.reset_index()
        df_list.append(collapsed_df)
    
    # outputting the final dataframes
    final_df = pd.concat(df_list)
    pre_final_df = pd.concat(pre_df_list)
    post_final_df = pd.concat(post_df_list)
    
    pre_final_df = pre_final_df.reset_index()
    post_final_df = post_final_df.reset_index()
        
    return final_df, pre_final_df, post_final_df

def number_of_collapses(allk,k298,floatf,floatm,lit):
    '''
    Counts the number of collapses occuring for each of the rate constant types.
    
    Args:
        allk: dataframe using the all data OH-NO rate constant
        k298: dataframe using the 298-258 K OH-NO rate constant
        floatf: dataframe using the OH-NO rate constant where Fcent is unfixed
        floatm: dataframe using the OH-NO rate constant where m is unfixed
        lit: dataframe using the JPL literature OH-NO rate constant
        
    Returns:
        dataframe containing the number of collapses and percentage change from the literature number. Grouped by temperature.
        
    '''
    
    # counting the number of collapse points grouped by temperature
    no_allk = allk.groupby('Temperature').count()
    no_298 = k298.groupby('Temperature').count()
    no_floatf = floatf.groupby('Temperature').count()
    no_floatm = floatm.groupby('Temperature').count()
    no_lit = lit.groupby('Temperature').count()
    
    # generating a new dataframe with the temperature and OH concentration of each rate constant
    df = pd.DataFrame({'Temperature':list(no_allk.index),
                      'all k': no_allk['OH concentration'],
                      '298': no_298['OH concentration'],
                        'float f': no_floatf['OH concentration'],
                       'float m': no_floatm['OH concentration'],
                       'lit': no_lit['OH concentration']})
    # dropping the second temperature column generated
    df = df.drop('Temperature', axis=1)
    # resetting the index
    df = df.reset_index()
    #
    literature_no = list(df['lit'])
    # list of unique rate constant types
    rate_list = ['all k','298','float f','float m']

    # looping through each rate constant type
    for i in rate_list:
        
        rate_number = list(df[i])
        # calculating the difference between the new rate constant number and the literature number
        diff = np.subtract(rate_number,literature_no)
        # appending the rate constant difference as a column
        df[(i+' difference')] = diff
        # calculating the percentage change
        percentage_change = diff/literature_no
        # appending the rate constant percentage change as a column
        df[(i+'%')] = percentage_change
        
    return df

def collapse_means(df_lit, df_collapse, df_pre_collapse, df_post_collapse):
    '''
    Calculate the mean value of NO, NO2, NMVOC and O3 before, during and after a collapse, in addition to annual stats.
    
    Args:
        - df_lit: dataframe containing the DEFRA measurements using the JPL literature rate constants
        - df_collapse: dataframe containing the data associated with points of collapse
        - df_pre_collapse: dataframe containing the data associated with points before collapse
        - df_post_collapse: dataframe containing the data associated with points after collapse
        
    Returns:
        - dataframe containing the mean NO, NO2, NMVOC sum and O3 for before, during, after a collapse and annually.
    '''
    
    # selecting the required columns
    pre_collapse = df_pre_collapse[['Nitric oxide','NMVOC sum','Nitrogen dioxide','Ozone']]
    post_collapse = df_post_collapse[['Nitric oxide','NMVOC sum','Nitrogen dioxide','Ozone']]
    collapse = df_collapse[['Nitric oxide','NMVOC sum','Nitrogen dioxide','Ozone']]
    df_lit = df_lit[['Nitric oxide','NMVOC sum','Nitrogen dioxide','Ozone']]
    
    # calculating stats
    pre_stats = pre_collapse.describe()
    post_stats = post_collapse.describe()
    col_stats = collapse.describe()
    lit_stats = df_lit.describe()
    # accessing mean value
    df_pre_stats = (pd.DataFrame(pre_stats.iloc[1])).reset_index()
    df_post_stats = (pd.DataFrame(post_stats.iloc[1])).reset_index()
    df_col_stats = (pd.DataFrame(col_stats.iloc[1])).reset_index()
    df_lit_stats = (pd.DataFrame(lit_stats.iloc[1])).reset_index()
    # adding a descriptor
    df_pre_stats['Statistic Type'] = 'Before'
    df_post_stats['Statistic Type'] = 'After'
    df_col_stats['Statistic Type'] = 'Collapse'
    df_lit_stats['Statistic Type'] = 'Annual'
    # creating a final dataframe
    complete_stats = pd.concat([df_pre_stats, df_col_stats, df_post_stats, df_lit_stats])
    
    return complete_stats

def AQ_limits(collapsed_df, df_year):
    '''
    Calculates the 8 hour mean of ozone and carbon monoxide and the hourly mean of nitrogen dioxide for comparison with DEFRA air quality limits.
    
    Args:
        collapsed_df: dataframe containing the collapse point data for chosen year
        df_year: dataframe containing complete measurement data for year chosen
        
    Returns:
        dataframe with the average 8 hour mean for ozone and carbon monoxide and the hour mean of nitrogen dioxide in molecules cm-3
    '''
    
    # making sure that the index is reset
    collapsed_df = collapsed_df.reset_index()
    # creating an empty dataframe
    total_8_hour = pd.DataFrame()

    # calculating the 8hr averages
    hour8 = []

    for z in range(len(collapsed_df)):
        # accessing each row
        row_df = collapsed_df.loc[z]

        # accessing rows 8 hours ahead of collapse time
        for i in range(1,8,1):
            hour_added = i
            hour8.append(row_df['Date_Time'] + pd.Timedelta(hours=hour_added))

        # filtering to include only rows within 8 hour period
        lit_8hour = df_year[df_year.Date_Time.isin(hour8)]
        # working out 8hr mean
        lit_8hour_mean = (pd.DataFrame(lit_8hour.mean(numeric_only=True))).T
        # filtering to relevant species only
        lit_8hour_mean = lit_8hour_mean[['Ozone','Carbon monoxide']]
        # adding to one complete dataframe
        total_8_hour = pd.concat([lit_8hour_mean, total_8_hour], axis = 0)
        # resetting list
        hour8 = []
        
    # filtering to NO2 only for 1 hour mean
    collapsed_df_NO2 = collapsed_df[['Nitrogen dioxide']]
    
    # Overall averages
    mean_8hr = pd.DataFrame(total_8_hour.mean(numeric_only=True))
    mean_1hr = pd.DataFrame(collapsed_df_NO2.mean(numeric_only=True))
    
    # adding to single dataframe
    total_df = pd.concat([mean_8hr, mean_1hr], axis = 0)
    total_df = total_df.rename(columns={0:'mean value / molecules cm-3'})
        
    return total_df