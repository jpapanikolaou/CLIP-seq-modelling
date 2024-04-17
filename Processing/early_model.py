import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/johnpapanikolaou/Columbia/spring-junior/Genomics/Project/Processing')
import importlib
import poissonFcns
importlib.reload(poissonFcns)
#%%
coverages = pd.read_csv('coverage.csv')

#%%
#convert coverages to a dict
coverages.set_index('chrom',inplace=True)
coverage_dict = {k: v[['index', 'coverage']] for k, v in coverages.groupby('chrom')}

#%% #apply global poisson
for key in coverage_dict.keys():
    coverage_dict[key] = poissonFcns.apply_poisson_distribution(coverage_dict[key])
# for key in coverage_dict.keys():
#     coverage_dict[key] = poissonFcns.add_poisson_p_value(coverage_dict[key])


#%% apply local poisson - observation - it 'smoothes out' values
for key in coverage_dict.keys():
    coverage_dict[key] = poissonFcns.apply_dynamic_poisson(coverage_dict[key])

#%% # get spread between global_poisson_prob and dynamic_poisson_probability
for key in coverage_dict.keys():
    coverage_dict[key] = poissonFcns.get_spread(coverage_dict[key],'global_poisson_prob',
                                    'dynamic_poisson_prob','global_local_poisson_spread'
                                    )


#%% apply multiple p value tests
for key in coverage_dict.keys():
    coverage_dict[key] = poissonFcns.apply_poisson_p_value_with_correction(coverage_dict[key])


#%% Plotting
plt.figure(figsize=(10,4))
plt.plot(coverage_dict['chr7']['index'],coverage_dict['chr7']['global_local_poisson_spread'])
#%%
