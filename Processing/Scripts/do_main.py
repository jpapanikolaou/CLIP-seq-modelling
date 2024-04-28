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
importlib.reload(binning_V3)
importlib.reload(macs_emulator)
importlib.reload(macs_comparison)
importlib.reload(compare_canonical)
#%%

def process_raw(path):
    #filter with samtools

    #bin
    binned = binning_V3.bin_reads_to_dataframe(path,50)
    binned_path = path+"_binned.csv"
    binned.to_csv(binned_path, index=False)
    binned.to_csv(binned_path, index=False)

    #run analysis
    df = pd.read_csv(binned_path)
    lambda_df = macs_emulator.calculate_excluded_lambdas(df)
    df_analyzed,significant_df, peak_group = macs_emulator.perform_peak_calling(lambda_df)
    df_analyzed.to_csv(path+"_our_peaks.csv", index=False)
    return df_analyzed