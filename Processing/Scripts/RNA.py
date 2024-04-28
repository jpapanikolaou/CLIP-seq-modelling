from Bio import  motifs
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import pyranges as pr

#%%
#
# peak_path = "../../ENCODE Data/eClip_control_and_target/our_target_peaks.csv"
# peaks = pd.read_csv(peak_path)
# significant_peaks = peaks[peaks['significant']==True]
# #%%
# ncrna_path = "../../oRNAment/Homo_sapiens_ncRNA_oRNAment.csv"
# ncrna_df = pd.read_csv(ncrna_path)
# #%%
# ncrna_df_columns = ncrna_df.columns.tolist()
# ncrna_df_columns[0] = "align_here"
# ncrna_df_columns[-1] = 'End'
# ncrna_df_columns[-2] = 'Start'
# ncrna_df_columns[-4] = 'Chromosome'
# ncrna_df.columns = ncrna_df_columns
#
# #%%
# int_to_string_path = "../../oRNAment/Homo_sapiens_string_to_int_ID_conversion.csv"
# int_to_string_df = pd.read_csv(int_to_string_path)
# int_to_string_column_names = int_to_string_df.columns.tolist()
# int_to_string_column_names[-1] = "align_here"
# int_to_string_column_names[-3] = "gene_name"
# int_to_string_df.columns = int_to_string_column_names
#
# #%%
# significant_peaks_columns = significant_peaks.columns.tolist()
# significant_peaks_columns[0] = 'Chromosome'
# significant_peaks_columns[1] = "Start"
# significant_peaks_columns[2] = "End"
# significant_peaks.columns = significant_peaks_columns
#
# #%% data cleaning
#
# ncrna_df = ncrna_df[~ncrna_df['Start'].astype(str).str.contains(',')]
# ncrna_df = ncrna_df[~ncrna_df['End'].astype(str).str.contains(',')]
#
# # Now convert 'Start' and 'End' columns to integers
# ncrna_df['Start'] = ncrna_df['Start'].astype(int)
# ncrna_df['End'] = ncrna_df['End'].astype(int)
#
# # Now you can proceed with converting to PyRanges
# ncrna_ranges = pr.PyRanges(ncrna_df)
# #%%
# significant_ranges = pr.PyRanges(significant_peaks)
#
# #%%
# overlap = significant_ranges.join(ncrna_ranges,how='left')
#
#
# #%%
#
# merged_df = pd.merge(overlap.df, int_to_string_df, on='align_here', how='inner')
#
# # Print the merged DataFrame
# print(merged_df)
#
# #%%
#
#
#
def do_merging(peak_path,ncrna_path,int_to_string_path):
    peaks = pd.read_csv(peak_path)
    significant_peaks = peaks[peaks['significant'] == True]

    #process ncrna
    # ncrna_path = "../../oRNAment/Homo_sapiens_ncRNA_oRNAment.csv"
    ncrna_df = pd.read_csv(ncrna_path)
    ncrna_df_columns = ncrna_df.columns.tolist()
    ncrna_df_columns[0] = "align_here"
    ncrna_df_columns[-1] = 'End'
    ncrna_df_columns[-2] = 'Start'
    ncrna_df_columns[-4] = 'Chromosome'
    ncrna_df.columns = ncrna_df_columns

    #process int_to_string
    # int_to_string_path = "../../oRNAment/Homo_sapiens_string_to_int_ID_conversion.csv"
    int_to_string_df = pd.read_csv(int_to_string_path)
    int_to_string_column_names = int_to_string_df.columns.tolist()
    int_to_string_column_names[-1] = "align_here"
    int_to_string_column_names[-3] = "gene_name"
    int_to_string_df.columns = int_to_string_column_names

    #prepare for pyranges
    significant_peaks_columns = significant_peaks.columns.tolist()
    significant_peaks_columns[0] = 'Chromosome'
    significant_peaks_columns[1] = "Start"
    significant_peaks_columns[2] = "End"
    significant_peaks.columns = significant_peaks_columns

    #execute pyranges
    ncrna_df = ncrna_df[~ncrna_df['Start'].astype(str).str.contains(',')]
    ncrna_df = ncrna_df[~ncrna_df['End'].astype(str).str.contains(',')]

    # Now convert 'Start' and 'End' columns to integers
    ncrna_df['Start'] = ncrna_df['Start'].astype(int)
    ncrna_df['End'] = ncrna_df['End'].astype(int)

    # Now you can proceed with converting to PyRanges
    print('data preprocessing done, now making pyranges objects')
    ncrna_ranges = pr.PyRanges(ncrna_df)
    significant_ranges = pr.PyRanges(significant_peaks)
    print("pyranges objects created")

    overlap = significant_ranges.join(ncrna_ranges,how='left')

    print('merging df (final step in function)')
    merged_df = pd.merge(overlap.df, int_to_string_df, on='align_here', how='inner')

    return merged_df

#%%
# peak_path = "../../ENCODE Data/eClip_control_and_target/our_target_peaks.csv"
# merged_df = do_merging(peak_path)
#%%
