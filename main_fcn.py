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
import compare_canonical
import RNA
import gc_normalization
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)
importlib.reload(compare_canonical)
importlib.reload(compare_canonical)
importlib.reload(RNA)

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
    return df_analyzed



#%%
control_file_path = "ENCODE Data/eClip_control_and_target/Control/ControlFiltering/control_filtered_reads.bam"
macs_control_path = "ENCODE Data/eClip_control_and_target/Control/control_output_peaks.csv"
control_df = run_entire_pipeline(control_file_path,macs_control_path)