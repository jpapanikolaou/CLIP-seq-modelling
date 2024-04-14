import pandas as pd
import matplotlib.pyplot as plt
#%% #read in data
binning_df = pd.read_csv('coverage.csv')
coverage = binning_df.iloc[:,1]
#%% #make a plot
plt.figure(figsize=(10,4))
plt.plot(coverage)
plt.ion()
plt.show()
print("HEllo")
#%%