#%%
import sys
import pandas as pd
import importlib
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/Scripts")
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/")
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Evaluation")
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/PlayingWithData")

#%%
import binning_V3
import macs_emulator
import macs_comparison
import compare_canonical
import RNA
import manual_gc
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)
importlib.reload(compare_canonical)
importlib.reload(compare_canonical)
importlib.reload(RNA)

#%%
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
               'chrX', 'chrY', 'chrM']

#%%

def run_entire_pipeline(experiment_file_path,macs_output_path):
    binned = binning_V3.bin_reads_to_dataframe(experiment_file_path,200)
    df = macs_emulator.calculate_excluded_lambdas(binned)
    df_analyzed,significant_df, peak_group = macs_emulator.perform_peak_calling(df)
    macs_data = macs_comparison.read_encode(macs_output_path)
    overlap_df = macs_comparison.find_overlaps_optimized(significant_df,macs_data)
    confusion_matrix = macs_comparison.build_confusion_matrix(significant_df,macs_data,overlap_df)
    print("Confusion matrix:")
    print("-----------------")
    print(confusion_matrix)
    return df_analyzed[df_analyzed['chromosome'].isin(chromosomes)]

def run_pipeline_gc_normalization(experiment_file_path,macs_output_path):
    # binned = binning_V3.bin_reads_to_dataframe(experiment_file_path,200)
    binned = manual_gc.merge_main(experiment_file_path)
    binned['count']  = binned['normalized_count']
    binned.rename(columns={'Chromosome':'chromosome','Start':'start','End':'end'},inplace=True)
    print("binned head")
    print(binned.head())
    df = macs_emulator.calculate_excluded_lambdas(binned)
    df_analyzed,significant_df, peak_group = macs_emulator.perform_peak_calling(df)
    macs_data = macs_comparison.read_encode(macs_output_path)
    overlap_df = macs_comparison.find_overlaps_optimized(significant_df,macs_data)
    confusion_matrix = macs_comparison.build_confusion_matrix(significant_df,macs_data,overlap_df)
    print("Confusion matrix:")
    print("-----------------")
    print(confusion_matrix)
    return df_analyzed[df_analyzed['chromosome'].isin(chromosomes)]

#%%
control_file_path = "ENCODE Data/eClip_control_and_target/Control/ControlFiltering/control_filtered_reads.bam"
macs_control_path = "ENCODE Data/eClip_control_and_target/Control/control_output_peaks.csv"
control_df = run_entire_pipeline(control_file_path,macs_control_path)
#%%
normalized_control_df = run_pipeline_gc_normalization(control_file_path,macs_control_path)

#%% spin up a q-q plot
import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt
import statsmodels.api as sm
chr1_df = control_df[control_df['chromosome']=='chr1']


chr1_df_sorted = chr1_df.sort_values(by='p_value',ascending=True)
p_values_sorted = chr1_df_sorted['p_value']
lambda_local_sorted = chr1_df_sorted['lambda_local']

n = len(p_values_sorted)

theoretical_quantiles = np.linspace(0,1,n)

observed_quantiles = np.interp(theoretical_quantiles,np.linspace(0,1,n),p_values_sorted)

#%%
plt.figure(figsize=(6, 6))
plt.plot(theoretical_quantiles, observed_quantiles, 'o', label='Observed Quantiles')
plt.plot([0, 1], [0, 1], 'r--', label='45-degree line')  # Ideal line for perfect uniform distribution
plt.title('Q-Q Plot against a Uniform Distribution')
plt.xlabel('Theoretical Quantiles (Uniform)')
plt.ylabel('Observed Quantiles')
plt.legend()
plt.grid(True)
plt.show()


#%% trying plot V2

chr1_df = chr1_df.copy()
chr1_df.loc[:, 'neg_log_p'] = -np.log10(chr1_df['p_value'])

n = len(chr1_df)
uniform_quantiles = np.linspace(1/(n+1), 1-1/(n+1), n)
expected_neg_log_p = -np.log10(uniform_quantiles)  # Use uniform quantiles for expected -log10 p-values

#%%
plt.figure(figsize=(8, 6))
plt.scatter(expected_neg_log_p, np.sort(chr1_df['neg_log_p']), color='blue', edgecolor='k', s=20)
plt.plot([0, max(expected_neg_log_p)], [0, max(expected_neg_log_p)], 'r--', label='Expected line')
plt.title('Q-Q Plot of -log10 p-values')
plt.xlabel('Expected -log10 p-values')
plt.ylabel('Observed -log10 p-values')
plt.grid(True)
plt.legend()
plt.show()