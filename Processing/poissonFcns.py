from scipy.stats import poisson
import numpy as np
from statsmodels.stats.multitest import multipletests


import pandas as pd
#%% Global poisson functions
def apply_poisson_distribution(df):
    coverages = df['coverage']
    lamda_val = np.mean(coverages)
    df['global_poisson_prob'] = poisson.pmf(coverages,lamda_val)
    return df

'''
def add_poisson_p_value(df):
    lambda_val = np.mean(df['coverage'])
    df['global_poisson_p_value'] = df['coverage'].apply(lambda x: 1 - poisson.cdf(x, lambda_val))
    return df
'''

#%% Local poisson functions
def apply_dynamic_poisson(df, window_size=50):
    """Applies Poisson PMF with a dynamic lambda based on a moving average of coverages."""
    # Calculate moving average of coverages as lambda
    lambda_vals = df['coverage'].rolling(window=window_size, min_periods=1).mean()
    df['dynamic_poisson_prob'] = [poisson.pmf(k, lam) for k, lam in zip(df['coverage'], lambda_vals)]
    return df

#%% Get spread

def get_spread(df,col1,col2,spread_name):
    df[spread_name] = df[col1]-df[col2]
    return df

#%% apply multiple p-val tests

def apply_poisson_p_value_with_correction(df):
    lambda_val = np.mean(df['coverage'])
    # Calculate original p-values
    p_values = df['coverage'].apply(lambda x: 1 - poisson.cdf(x, lambda_val))
    # Store original p-values in the DataFrame
    df['global_poisson_p_value'] = p_values
    # Flag significant original p-values
    df['original_is_significant'] = p_values < 0.05
    # Apply Bonferroni correction
    corrected_p_values = multipletests(p_values, method='bonferroni')
    df['global_poisson_p_value_corrected'] = corrected_p_values[1]
    df['corrected_is_significant'] = corrected_p_values[0]
    return df#%% compare correction to original
#%%
def compare_original_p_vals_with_correction(df):
    original = df['global_poisson_p_value']
    corrected = df['global_poisson_p_value_corrected']
    return original, corrected
