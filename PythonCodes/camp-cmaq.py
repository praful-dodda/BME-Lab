import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PseudoNetCDF as pnc
import os

from scipy import interpolate
import seaborn as sns
sns.set()

cmaq_data_folder = "DATA/CMAQ/2019/"

cmaq_file_path = os.path.join(cmaq_data_folder, "HR2DAY_LST_ACONC_EQUATES_v532_12US1_2019.nc")
# read cmaq data into pnc 
cmaq_data = pnc.pncopen(cmaq_file_path , format='ioapi')

aqf = cmaq_data.subsetVariables(['PM25_AVG', 'PM10_AVG'])

aqs_path = os.path.join("DATA", "AQS","daily_88101_2019.csv")

aqsdata = pd.read_csv(aqs_path, parse_dates = ["Date Local"])

def calc_stats(df, obs_col, model_col, model_name = None):
    """
    This function calculates statistics like R2, RMSE, MAE, correlation, and bias
    and returns as a formatted dataframe
    """
    stats = pd.DataFrame()
    # stats['R2'] = [df[model_col].corr(df[obs_col])]
    
    stats['RMSE'] = [np.sqrt(((df[model_col] - df[obs_col])**2).mean())]
    # variance of error
    stats["SE"] = (df[model_col] - df[obs_col]).std()
    stats["ME"] = (df[model_col] - df[obs_col]).mean()
    # stats['MSE'] = [((df[model_col] - df[obs_col])**2).mean()]
    # stats['MAE'] = [np.abs(df[model_col] - df[obs_col]).mean()]
    stats['Correlation'] = [df[model_col].corr(df[obs_col], method='pearson')]
    # no. of observations
    stats['NO'] = [len(df)]

    # mean of observations
    stats['MO'] = [df[obs_col].mean()]
    # sd of observations
    stats['SO'] = [df[obs_col].std()]

    # same for model
    stats['MZ'] = [df[model_col].mean()]
    stats['SZ'] = [df[model_col].std()]

    if model_name:
        stats.index = [model_name]

    return stats


# Create a column with 8 digit Alphanumeric code called "Site ID" similar to "AQS12345" for each unique combination of Longitude and Latitude
aqsdata["Site ID"] = aqsdata[["Latitude", "Longitude"]].apply(lambda x: "AQS" + str(int(x[0]*10000)), axis=1)

starttime = aqf.getTimes()[0].replace(tzinfo=None)

time = aqsdata['Date Local'].dt.to_pydatetime()

dt = np.array([(t - starttime).total_seconds() // (24*3600) for t in time], 'i')

xx = np.array([(t - starttime).days for t in time], 'i')
yy = np.array([(t - starttime).total_seconds() // (24*3600) for t in time], 'i')

lat = aqsdata.Latitude.values
lon = aqsdata.Longitude.values

# Convert lon/lat to indices and add a surface index
aqsdata['I'], aqsdata['J'] = aqf.ll2ij(lon, lat)
aqsdata['K'] = aqsdata['I'] * 0
aqsdata['T'] = dt

# Only include data in time/space domain
nx = aqf.NCOLS
ny = aqf.NROWS
nt = len(aqf.dimensions['TSTEP'])
aqssubset = aqsdata.query(
    f'I > 0 and I < {nx} and J > 0 and J < {ny} and T > 0 and T < {nt}'
).copy()

# Index model ozone and add to obs DataFrame
i = aqssubset['I'].values
j = aqssubset['J'].values
k = aqssubset['K'].values
t = aqssubset['T'].values

aqssubset['PM25_MOD'] = aqf.variables['PM25_AVG'][t, k, j, i]

renamed_df = aqssubset.rename(columns={"PM25_MOD": "CMAQ_PM2.5", "Arithmetic Mean": "AQS_PM2.5"}).copy()

camp_df = renamed_df[['CMAQ_PM2.5', 'AQS_PM2.5']].copy()

def CMAQ_CAMP(camp_df: pd.DataFrame,log_transform = False,
              plot = "raw") -> pd.DataFrame:
    """
    This function applies CAMP-Correction to the CMAQ data and plots the results.
    The inputs are:
    1. camp_df: A dataframe with the CMAQ and AQS data
    2. obs_col: The name of the column containing the AQS data
    3. model_col: The name of the column containing the CMAQ data
    4. log_transform: If True, the data is log transformed before applying CAMP
    5. plot: If True, the results are plotted
    6. plot_name: The name of the plot file
    """
    df = camp_df.copy()
    n_cut = 10
    # Stratify the data into 10 percentile bins of increasing predicted values
    df['bin'] = pd.qcut(df['CMAQ_PM2.5'], n_cut, labels=False)

    # Calculate the mean and variance of observed values for each bin
    obs_mean = df.groupby('bin')['AQS_PM2.5'].mean() # equivalent to lambda_1 in the paper
    obs_var = df.groupby('bin')['AQS_PM2.5'].var()
    pred_mean = df.groupby('bin')['CMAQ_PM2.5'].mean()
    pred_var = df.groupby('bin')['CMAQ_PM2.5'].var()
    l_1 = np.array([])
    l_2 = np.array([])
    for each_bin in range(n_cut):
        # calculate the mean and variance of the observed PM2.5 in each bin
        obs_mean[each_bin] = df.query(f'bin == {each_bin}')['AQS_PM2.5'].mean()
        obs_var[each_bin] = df.query(f'bin == {each_bin}')['AQS_PM2.5'].var()
        values_in_each_bin = df.query(f'bin == {each_bin}')['AQS_PM2.5'].count()

        lambda_1 = (1/values_in_each_bin) * df.query(f'bin == {each_bin}')['AQS_PM2.5'].sum()
        lambda_2 = (1/(values_in_each_bin - 1)) * np.sum((df.query(f'bin == {each_bin}')['AQS_PM2.5']  - lambda_1)**2)
        
        l_1 = np.append(l_1, lambda_1)
        l_2 = np.append(l_2, lambda_2)

    if plot == "raw":
        raw_stats = calc_stats(df, 'CMAQ_PM2.5', 'AQS_PM2.5', 'AQS-CMAQ')
        plt.figure(figsize=(14, 9))
        # Create a scatter plot of the observation-prediction pairs
        plt.scatter(df['CMAQ_PM2.5'], df['AQS_PM2.5'], label='Observation-Prediction Pairs',
                    alpha=0.02, s=0.2, color='k')

        plt.xlim(0, 40)
        plt.ylim(0, 40)

        # Add bin lines to the plot
        for i in range(10):
            plt.axvline(df.groupby('bin')['CMAQ_PM2.5'].min()[i], ls='--', c='gray')
            plt.text(df.groupby('bin')['CMAQ_PM2.5'].median()[i], plt.ylim()[1], f'bin {i+1}', rotation=90, va='top')

        # Add 1:1 line
        plt.plot([0, 40], [0, 40], color='r', label='1:1 Line', linestyle='dotted')
        plt.plot(df.groupby('bin')['CMAQ_PM2.5'].median().values, 
                l_1, 'o-', color = "b", label = "lambda_1")
        # add points

        plt.plot(df.groupby('bin')['CMAQ_PM2.5'].median().values, 
                l_2, 'o-', color = "g", label = "lambda_2")
        # Add axis labels and a legend
        plt.xlabel('Modeled PM2.5 Concentration')
        plt.ylabel('Observed PM2.5 Concentration')
        plt.legend(loc = "upper right")

        plt.annotate(f"R$^2$ = {raw_stats['Correlation'].values[0]:.2f}", (0.75, 0.55), xycoords='axes fraction', color='k')
        plt.annotate(f"RMSE = {raw_stats['RMSE'].values[0]:.2f}", (0.75, 0.50), xycoords='axes fraction', color='k')
        plt.annotate(f"ME = {raw_stats['ME'].values[0]:.2f}", (0.75, 0.45), xycoords='axes fraction', color='k')
        plt.annotate(f"SE = {raw_stats['SE'].values[0]:.2f}", (0.75, 0.40), xycoords='axes fraction', color='k')
        plt.annotate(f"NO = {raw_stats['NO'].values[0]:.0f}", (0.75, 0.35), xycoords='axes fraction', color='k')
        plt.annotate(f"MO = {raw_stats['MO'].values[0]:.2f}", (0.75, 0.30), xycoords='axes fraction', color='k')
        plt.annotate(f"SO = {raw_stats['SO'].values[0]:.2f}", (0.75, 0.25), xycoords='axes fraction', color='k')
        plt.annotate(f"MZ = {raw_stats['MZ'].values[0]:.2f}", (0.75, 0.20), xycoords='axes fraction', color='k')
        plt.annotate(f"SZ = {raw_stats['SZ'].values[0]:.2f}", (0.75, 0.15), xycoords='axes fraction', color='k')

        # Show the plot
        plt.show()
    else:
        df['custom_camp'] = interpolate.interp1d(pred_mean, l_1, fill_value="extrapolate")(df['CMAQ_PM2.5'])
        camp_corr_stats = calc_stats(df, 'custom_camp', 'AQS_PM2.5', 'CAMP Corrected')
        plt.figure(figsize=(14, 9))
        # Create a scatter plot of the observation-prediction pairs
        plt.scatter(df['CMAQ_PM2.5'], df['AQS_PM2.5'], label='Observation-Prediction Pairs',
                    alpha=0.02, s=0.2, color='k')

        plt.scatter(df['custom_camp'], df['AQS_PM2.5'], label='CAMP Corrected PM2.5',
                    alpha=0.02, s=0.2, color='r')

        plt.xlim(0, 40)
        plt.ylim(0, 40)

        # Add bin lines to the plot
        for i in range(10):
            plt.axvline(df.groupby('bin')['CMAQ_PM2.5'].min()[i], ls='--', c='gray')
            plt.text(df.groupby('bin')['CMAQ_PM2.5'].median()[i], plt.ylim()[1], f'bin {i+1}', rotation=90, va='top')

        # Add 1:1 line
        plt.plot([0, 40], [0, 40], color='r', label='1:1 Line', linestyle='dotted')
        plt.plot(df.groupby('bin')['CMAQ_PM2.5'].median().values, 
                l_1, 'o-', color = "b", label = "lambda_1")
        # add points

        plt.plot(df.groupby('bin')['CMAQ_PM2.5'].median().values, 
                l_2, 'o-', color = "g", label = "lambda_2")
        # Add axis labels and a legend
        plt.xlabel('Predicted PM2.5 Concentration')
        plt.ylabel('Observed PM2.5 Concentration')
        plt.legend(loc = "upper right")

        # annotate with all statistics calculated and saved in camp_corr_stats in the bottom right corner
        plt.annotate(f"R$^2$ = {camp_corr_stats['Correlation'].values[0]:.2f}", (0.75, 0.55), xycoords='axes fraction', color='r')
        plt.annotate(f"RMSE = {camp_corr_stats['RMSE'].values[0]:.2f}", (0.75, 0.50), xycoords='axes fraction', color='r')
        plt.annotate(f"ME = {camp_corr_stats['ME'].values[0]:.2f}", (0.75, 0.45), xycoords='axes fraction', color='r')
        plt.annotate(f"SE = {camp_corr_stats['SE'].values[0]:.2f}", (0.75, 0.40), xycoords='axes fraction', color='r')
        plt.annotate(f"NO = {camp_corr_stats['NO'].values[0]:.0f}", (0.75, 0.35), xycoords='axes fraction', color='r')
        plt.annotate(f"MO = {camp_corr_stats['MO'].values[0]:.2f}", (0.75, 0.30), xycoords='axes fraction', color='r')
        plt.annotate(f"SO = {camp_corr_stats['SO'].values[0]:.2f}", (0.75, 0.25), xycoords='axes fraction', color='r')
        plt.annotate(f"MZ = {camp_corr_stats['MZ'].values[0]:.2f}", (0.75, 0.20), xycoords='axes fraction', color='r')
        plt.annotate(f"SZ = {camp_corr_stats['SZ'].values[0]:.2f}", (0.75, 0.15), xycoords='axes fraction', color='r')


        # Show the plot
        plt.show()