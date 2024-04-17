import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
#%%
coverage_df = pd.read_csv('../coverage.csv')
coverage_df.drop('Unnamed: 0',axis=1,inplace=True)
#%%
def plot_chrom_df(chrom_df,chrom_name):
    df = chrom_df[chrom_df['chrom']==chrom_name]
    plt.figure(figsize=(10,4))
    plt.title(chrom_name)
    plt.plot(df[['index','coverage']])

#%%
chr_names = [coverage_df['chrom'].unique()][0]

for chrom in chr_names:
    plot_chrom_df(coverage_df,chrom)

#%%
plots_dir = '../Coverage_Plots'
def plot_chrom_df(chrom_df, chrom_name, ax):
    df = chrom_df[chrom_df['chrom'] == chrom_name]
    ax.plot(df['index'], df['coverage'], label=chrom_name)
    ax.set_title(chrom_name)
    ax.set_xlabel('Bin Index')
    ax.set_ylabel('Coverage')

# Get the unique chromosome names from the dataframe
chr_names = coverage_df['chrom'].unique()

# Define the number of chromosomes per figure (up to 5)
chroms_per_fig = 5

# Calculate the number of figures needed
num_figs = len(chr_names) // chroms_per_fig + (len(chr_names) % chroms_per_fig > 0)

# Plot 5 chromosomes per figure
for fig_num in range(num_figs):
    plt.figure(figsize=(10, 20))  # Define a new figure

    # Determine the chromosomes for the current figure
    start_idx = fig_num * chroms_per_fig
    end_idx = min(start_idx + chroms_per_fig, len(chr_names))

    # Create a subplot for each chromosome in the current figure
    for i, chrom in enumerate(chr_names[start_idx:end_idx], start=1):
        ax = plt.subplot(chroms_per_fig, 1, i)
        plot_chrom_df(coverage_df, chrom, ax)

    plt.tight_layout()

    # Define the path for the current figure and save it
    fig_path = os.path.join(plots_dir, f'coverage_plot_{fig_num + 1}.png')
    plt.savefig(fig_path)
    plt.close()  # Close the figure after saving to free up memory