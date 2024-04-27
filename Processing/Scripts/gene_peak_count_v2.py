import pandas as pd

def count_peaks_in_genes(peaks_file_path, genes_file_path):
    # Load data into pandas DataFrames
    peaks_df = pd.read_csv(peaks_file_path)
    genes_df = pd.read_csv(genes_file_path)

    # Make sure data types are efficient and chromosome columns are strings
    peaks_df[['chromosome', 'start', 'end']] = peaks_df[['chromosome', 'start', 'end']].astype({'chromosome': 'str', 'start': 'int', 'end': 'int'})
    genes_df[['chromosome', 'start', 'end']] = genes_df[['chromosome', 'start', 'end']].astype({'chromosome': 'str', 'start': 'int', 'end': 'int'})

    # Initialize result dataframe
    gene_peak_counts = pd.DataFrame(columns=['Gene', 'Peak Count'])

    # Group by chromosome for both datasets to reduce comparison scope
    grouped_peaks = peaks_df.groupby('chromosome')
    grouped_genes = genes_df.groupby('chromosome')

    # Iterate over chromosomes that exist in both peak and gene data
    for chromosome in set(grouped_peaks.groups).intersection(set(grouped_genes.groups)):
        print("Beginning chromosome", chromosome)
        peaks = grouped_peaks.get_group(chromosome)
        genes = grouped_genes.get_group(chromosome)

        # Use an interval index for genes
        genes_interval_index = pd.IntervalIndex.from_arrays(genes['start'], genes['end'], closed='both')

        # Check each peak against the gene intervals
        for _, peak in peaks.iterrows():
            peak_interval = pd.Interval(peak['start'], peak['end'], closed='both')

            # Find overlaps
            overlapping_genes = genes[genes_interval_index.overlaps(peak_interval)]

            # Increment counts
            for gene in overlapping_genes['Gene'].unique():
                if gene in gene_peak_counts.index:
                    gene_peak_counts.loc[gene] += 1
                else:
                    gene_peak_counts.loc[gene] = 1

    # Convert the series to a DataFrame and reset the index
    gene_peak_counts = gene_peak_counts.reset_index()
    gene_peak_counts.columns = ['Gene', 'Peak Count']

    # Sort the DataFrame by peak count in descending order
    gene_peak_counts = gene_peak_counts.sort_values(by='Peak Count', ascending=False)

    return gene_peak_counts
