import pandas as pd

def count_peaks_in_genes(peaks_file_path, genes_file_path):
    # Load significant peaks and gene datasets into pandas DataFrames
    peaks_df = pd.read_csv(peaks_file_path)
    genes_df = pd.read_csv(genes_file_path, sep="\t")

    # Define a function to check if a peak overlaps with a gene
    def peak_overlaps_gene(peak_start, peak_end, gene_start, gene_end):
        return (peak_start <= gene_end) and (peak_end >= gene_start)

    # Initialize a dictionary to store the count of peaks associated with each gene
    gene_peak_counts = {}

    # Iterate through each peak
    for index, peak in peaks_df.iterrows():
        peak_chromosome = peak['chromosome']
        peak_start = int(peak['start'])
        peak_end = int(peak['end'])

        # Iterate through each gene
        for _, gene in genes_df.iterrows():
            gene_chromosome = gene[0]
            gene_start = int(gene[1])
            gene_end = int(gene[2])

            # Check if the peak overlaps with the gene
            if peak_chromosome == gene_chromosome and peak_overlaps_gene(peak_start, peak_end, gene_start, gene_end):
                gene_symbol = gene[3]

                # Increment the count of peaks associated with the gene
                gene_peak_counts[gene_symbol] = gene_peak_counts.get(gene_symbol, 0) + 1

    # Convert the dictionary to a DataFrame for easier analysis
    gene_peak_counts_df = pd.DataFrame(list(gene_peak_counts.items()), columns=['Gene Symbol', 'Peak Count'])

    # Sort the DataFrame by peak count in descending order
    gene_peak_counts_df = gene_peak_counts_df.sort_values(by='Peak Count', ascending=False)

    return gene_peak_counts_df

# # Example usage:
# peaks_file_path = "significant_our_control.csv"
# genes_file_path = "hgTablesAllGenesGeneSymbolFilt2.txt"
# result_df = count_peaks_in_genes(peaks_file_path, genes_file_path)
# print(result_df)
