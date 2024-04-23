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
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)

#%% run binning_V3

control_file_path = "ENCODE Data/eClip_control_and_target/Control/ENCFF913CBD.bam"
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

# method to generate overlap_df will take a bit of time ~1-2 minutes
overlap_df = macs_comparison.find_overlaps_optimized(our_control_data,macs_control_data)
confusion_matrix = macs_comparison.build_confusion_matrix(our_control_data,macs_control_data,overlap_df)