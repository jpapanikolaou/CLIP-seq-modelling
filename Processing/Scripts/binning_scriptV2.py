import pysam
import numpy as np
import pandas as pd
def analyze_bam_coverage(bam_path, output_path, output_file_name):
    bin_size = 100000  # Set the bin size

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

