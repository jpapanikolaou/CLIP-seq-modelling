import sys
sys.path.append('/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing/Scripts')
import importlib
import early_model_script
importlib.reload(early_model_script)
#%% eClip_control

control_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Control/control.csv"
eClip_control_df = early_model_script.process_coverage_data(control_loc)

#%% eClip_target
target_loc = "/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/bioinformatics_28_23_3013_s20/data_found/ENCODE Data/eClip_control_and_target/Target/target.csv"
eClip_target_df = early_model_script.process_coverage_data(target_loc)
#%%