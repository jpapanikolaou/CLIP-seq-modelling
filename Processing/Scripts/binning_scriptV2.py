import pysam
import numpy as np
import pandas as pd
def analyze_bam_coverage(bam_path, output_path, output_file_name):
    bin_size = 10000  # Set the bin size

    # Open BAM file
    bam = pysam.AlignmentFile(bam_path, 'rb')

    chromosome_coverage = {}
    total_count = 0

    # Process reads
    for read in bam.fetch():
        if not read.is_unmapped and not read.is_duplicate:  # Exclude unmapped and duplicate reads
            total_count += 1
            if total_count % 100000 == 0:
                print(f"Processed {total_count} reads")

            chrom = bam.get_reference_name(read.reference_id)
            chrom_length = bam.get_reference_length(chrom)
            num_bins = int(np.ceil(chrom_length / bin_size))

            # Initialize the chromosome array if not already
            if chrom not in chromosome_coverage:
                chromosome_coverage[chrom] = np.zeros(num_bins, dtype=int)

            # Calculate the bin index and increment the coverage
            bin_index_start = read.reference_start // bin_size
            bin_index_end = read.reference_end // bin_size
            if bin_index_end >= num_bins:
                bin_index_end = num_bins - 1  # Ensure the index is within the range
            chromosome_coverage[chrom][bin_index_start:bin_index_end + 1] += 1
    for chrom in chromosome_coverage:
        indices = range(len(chromosome_coverage[chrom]))
        counts = chromosome_coverage[chrom]
        chromosome_coverage[chrom] = pd.DataFrame({'index': indices, 'counts': counts})
    return chromosome_coverage

def save_dict_to_csv_filtered(chromosome_coverage, output_file): #WARNING: this is for humans
    dfs = []
    for chromosome, df in chromosome_coverage.items():
        df['chromosome'] = chromosome  # Add a column for the chromosome name
        dfs.append(df)

    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(dfs, ignore_index=True)

    # Reorder columns if needed to have 'chromosome' as the first column
    combined_df = combined_df[['chromosome', 'index', 'counts']]

    # filter combined_df
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                   'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                   'chrX', 'chrY', 'chrM']
    combined_df = combined_df[combined_df['chromosome'].isin(chromosomes)]

    # Write the DataFrame to a CSV file
    combined_df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")