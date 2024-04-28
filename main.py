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
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)
importlib.reload(compare_canonical)
importlib.reload(compare_canonical)
importlib.reload(RNA)

#%% run binning_V3 on control

control_file_path = "ENCODE Data/eClip_control_and_target/Control/ControlFiltering/control_filtered_reads.bam"
control_binned = binning_V3.bin_reads_to_dataframe(control_file_path,250)
control_binned_path = "ENCODE Data/eClip_control_and_target/control_binned.csv"
control_binned.to_csv(control_binned_path, index=False)

#%% running binning_V3 on target
target_file_path = "ENCODE Data/eClip_control_and_target/Target/TargetFiltering/target_filtered_reads.bam"
target_binned = binning_V3.bin_reads_to_dataframe(target_file_path,200)
target_binned_path = "ENCODE Data/eClip_control_and_target/target_binned.csv"
target_binned.to_csv(target_binned_path, index=False)

#%% run macs_emulator

control_df = pd.read_csv(control_binned_path)
control_lambda_df = macs_emulator.calculate_excluded_lambdas(control_df)
control_df_analyzed,control_significant_df, control_peak_group = macs_emulator.perform_peak_calling(control_lambda_df)

control_df.to_csv("ENCODE Data/eClip_control_and_target/our_control_peaks.csv", index=False)

target_df = pd.read_csv(target_binned_path)
target_lambda_df = macs_emulator.calculate_excluded_lambdas(target_df)
target_df_analyzed,target_significant_df, target_peak_groups = macs_emulator.perform_peak_calling(target_lambda_df)


# control_df_analyzed.to_csv("ENCODE Data/eClip_control_and_target/our_control_peaks.csv", index=False)
target_df_analyzed.to_csv("ENCODE Data/eClip_control_and_target/our_target_peaks.csv", index=False)

#%% run macs_comparison
# load data
macs_control_path = "ENCODE Data/eClip_control_and_target/Control/control_output_peaks.csv"
macs_control_data = macs_comparison.read_encode(macs_control_path)
macs_control_data.columns = ['chr','start','end','length','abs_summit','pileup','-log10(p_value)',
                             'fold_enrichment','-log10(q_value)','name']
#%%
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
gene_count_path = "GenomeBrowserData/hgTableFilt.txt"
target_peak_path = "ENCODE Data/eClip_control_and_target/our_target_peaks.csv"
control_peak_path = "ENCODE Data/eClip_control_and_target/our_control_peaks.csv"
target_overlap_df = compare_canonical.compare_canonical(gene_count_path,target_peak_path)
control_overlap_df = compare_canonical.compare_canonical(gene_count_path,control_peak_path)
#%% do RNA annotations

ncrna_path = "oRNAment/Homo_sapiens_ncRNA_oRNAment.csv"
int_to_string_path = "/oRNAment/Homo_sapiens_string_to_int_ID_conversion.csv"
merged_df = RNA.do_merging(target_peak_path,ncrna_path,int_to_string_path)
