# read in data
import pandas as pd
import sys
import importlib
sys.path.append('/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/Scripts')
import model_funcs
import poisson_peaks
importlib.reload(model_funcs)
importlib.reload(poisson_peaks)
#%%
control_data = pd.read_csv("../ENCODE Data/eClip_control_and_target/Control/control_results.csv")
target_data = pd.read_csv("../ENCODE Data/eClip_control_and_target/Target/target_results.csv")
control_data_stats = control_data.copy()
target_data_stats = target_data.copy()
#%%
all_chr1_peaks,chr1_lambda_global = poisson_peaks.main_peak_calling(control_data_stats,'chr1')

#%%
poisson_peaks.create_poisson_qq_plot(control_data_stats, all_chr1_peaks, chr1_lambda_global)
#%%

chr_1_data = poisson_peaks.subset_chromosome(control_data_stats, 'chr1')
chr_1_data_no_zeroes = poisson_peaks.remove_zeroes(chr_1_data)
lambda_1 = poisson_peaks.calculate_lambda(chr_1_data_no_zeroes, set())
print('lambda 1:',lambda_1)
lambda_2 = poisson_peaks.calculate_lambda(chr1_data,set())
print('lambda 2:',lambda_2)

#%%

frequency_zeroes = chr_1_data['counts'].value_counts()
frequency_no_zeroes = chr_1_data_no_zeroes['counts'].value_counts()