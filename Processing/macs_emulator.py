import pandas as pd
from scipy.stats import poisson
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

#%%

# def perform_peak_calling(df,window_size = 500):
#     df = df.sort_values(by=['chromosome','start'])
#     df['local_background'] = df.groupby('chromosome')['count'].transform(lambda x: x.rolling(window_size,center=True,min_periods=1).mean())
#     df['p_value'] = poisson.sf(df['count'] - 1, df['local_background'])
#     significance_level = 0.05 /len(df) #bonferonni correction
#     df['significant'] = df['p_value'] < significance_level
#     return df

#%%

def calculate_excluded_lambdas(df, bin_size=50, local_bp=1000):
    # Define window sizes in terms of number of bins
    window_size_1k = local_bp // bin_size
    window_size_5k = 5000 // bin_size
    window_size_10k = 10000 // bin_size

    # Global lambda: mean across the data
    df['lambda_BG'] = df['count'].mean()

    # Adjusting the window size for accurate exclusion of the current bin
    def adjusted_rolling(series, window_size):
        shifted_up = series.shift(1)
        shifted_down = series.shift(-1)
        mean_up = shifted_up.rolling(window=window_size, min_periods=1, center=True).mean()
        mean_down = shifted_down.rolling(window=window_size, min_periods=1, center=True).mean()
        return (mean_up + mean_down) / 2  # Average the shifted rolling means

    # Applying the adjusted rolling calculation for each scale
    df['lambda_1k'] = df.groupby('chromosome')['count'].transform(lambda x: adjusted_rolling(x, window_size_1k))
    df['lambda_5k'] = df.groupby('chromosome')['count'].transform(lambda x: adjusted_rolling(x, window_size_5k))
    df['lambda_10k'] = df.groupby('chromosome')['count'].transform(lambda x: adjusted_rolling(x, window_size_10k))

    # Using max of global lambda and each local lambda as the final lambda for the bin
    df['lambda_local'] = df[['lambda_BG', 'lambda_1k', 'lambda_5k', 'lambda_10k']].max(axis=1)

    return df


def perform_peak_calling(df):
    # Vectorized calculation of p-values using Poisson survival function
    df['p_value'] = poisson.sf(df['count'] - 1, df['lambda_local'])

    adjusted_p_vals = multipletests(df['p_value'], method='fdr_bh')[1]
    df['adjusted_p_value'] = adjusted_p_vals
    df['significant'] = df['adjusted_p_value'] < 0.05

    # This defines peaks as continuous regions of significance
    df['peak'] = (df['significant'].astype(int).diff().ne(0).cumsum())
    significant_df = df[df['significant']]

    peak_groups = significant_df.groupby('peak')['count'].agg(['mean', 'max', 'size'])
    peak_groups = peak_groups[peak_groups['size'] > 1]
    return significant_df, peak_groups

#%%
control_df = pd.read_csv("../ENCODE Data/eClip_control_and_target/control_binned.csv")

#%%
lambda_df = calculate_excluded_lambdas(control_df)
#%%
significant_df, peak_groups = perform_peak_calling(lambda_df)


#%% plot p and q values
import matplotlib.pyplot as plt
import numpy as np
df = lambda_df
df['-log10(p_value)'] = -np.log10(df['p_value'])
# plt.figure(figsize=(12, 6))
# plt.scatter(df['start'], df['-log10(p_value)'], c=df['significant'], cmap='viridis', alpha=0.5)
# plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')  # Threshold line for significance
# plt.title('Manhattan Plot of Peaks')
# plt.xlabel('Genomic Position')
# plt.ylabel('-log10(p-value)')
# plt.show()