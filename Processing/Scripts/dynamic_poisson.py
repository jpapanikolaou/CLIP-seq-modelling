import pandas as pd
import numpy as np
from scipy.stats import nbinom
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
#%%
control_data = pd.read_csv("ENCODE Data/eClip_control_and_target/Control/control_results.csv")
target_data = pd.read_csv("ENCODE Data/eClip_control_and_target/Target/target_results.csv")
control_data_stats = control_data.copy()
target_data_stats = target_data.copy()
#%%

def calculate_local_lambda(data, window_size):
    """
    Calculate the local lambda for each bin using a rolling mean, excluding the bin itself.
    """
    # Ensure rolling mean doesn't include the current row in its calculation
    rolled_means = data['counts'].rolling(window=window_size, min_periods=1, center=True).mean()
    adjusted_lambda = rolled_means.shift(1).combine_first(rolled_means.shift(-1))
    return adjusted_lambda

def peak_detection(df, window_size=10000, alpha=0.05): #window size of 10000 matches MACS
    """
    Identify peaks based on local lambda calculations and adjust for multiple testing using FDR.
    """
    data = df.copy()  # Work on a copy to ensure no side effects or warnings
    data['local_lambda'] = calculate_local_lambda(data, window_size)
    data['p_value'] = data.apply(
        lambda row: nbinom.sf(row['counts'] - 1, row['local_lambda'], 0.5) if row['local_lambda'] > 0 else 1,
        axis=1
    )

    # Apply FDR control using Benjamini-Hochberg
    _, adjusted_p_values, _, _ = multipletests(data['p_value'], alpha=alpha, method='fdr_bh')
    data['adjusted_p_value'] = adjusted_p_values
    data['is_significant'] = data['adjusted_p_value'] < alpha

    data['FDR'] = np.nan
    data['is_significant'] = data['adjusted_p_value'] < alpha

    return data


#%%
chr_1_data = control_data_stats[control_data_stats['chromosome'] == 'chr1']
results = peak_detection(chr_1_data)
#%%
