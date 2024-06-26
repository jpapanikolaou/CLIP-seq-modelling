import os
import sys
sys.path.append('/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/Scripts')
import binning_script
import importlib
importlib.reload(binning_script)
os.listdir()
#%%
#For eClip_control_and_target_CONTROL
# data_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Control/ENCFF913CBD.bam"
# output_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Control"

#%% for eClip_control_and_target_TARGET
data_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Target/ENCFF862SNV.bam"
output_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Target"

binning_script.analyze_bam_coverage(bam_path=data_loc, output_path=output_loc, output_file_name='target.csv')
#%%
