import sys
import pandas as pd
import importlib
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/Scripts")
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/")
sys.path.append("/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Evaluation")

#%%
import binning_V3
import macs_emulator
import macs_comparison
import gene_peak_count
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)

#%% run binning_V3

control_file_path = "ENCODE Data/eClip_control_and_target/Control/ControlFiltering/control_filtered_reads.bam"
control_binned = binning_V3.bin_reads_to_dataframe(control_file_path,50)
control_binned.to_csv("ENCODE Data/eClip_control_and_target/control_binned.csv", index=False)

#%% run macs_emulator

control_df = pd.read_csv("ENCODE Data/eClip_control_and_target/control_binned.csv")
lambda_df = macs_emulator.calculate_excluded_lambdas(control_df)
regular_df,significant_df, peak_groups = macs_emulator.perform_peak_calling(lambda_df)
regular_df.to_csv("ENCODE Data/eClip_control_and_target/our_control_peaks.csv", index=False)

#%% run macs_comparison
# load data
macs_control_path = "ENCODE Data/eClip_control_and_target/Control/control_output_peaks.csv"
macs_control_data = macs_comparison.read_encode(macs_control_path)
our_data_path = "ENCODE Data/eClip_control_and_target/our_control_peaks.csv"
our_control_data = pd.read_csv(our_data_path)
significant_our_control_data = our_control_data[our_control_data['significant']==True]
# note - overlap_df can have more tuples than one of the participating dataframes,
# if there are multiple overlaps within a given tuple
overlap_df = macs_comparison.find_overlaps_optimized(significant_our_control_data,macs_control_data)
confusion_matrix = macs_comparison.build_confusion_matrix(significant_our_control_data,macs_control_data,overlap_df)

#%% see how we did relative to macs :)
print(confusion_matrix)

#%% convert peaks to gene counts
gene_count_path = "ENCODE Data/eClip_control_and_target/hgTablesCanonicalFilt.txt"
peak_path = "ENCODE Data/eClip_control_and_target/test2_peaks_csv.csv"
peak_gene_count = gene_peak_count(peak_path,gene_count_path)
