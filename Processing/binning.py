import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%

bam_path = '../bioinformatics_28_23_3013_s20/data_found/ENCODE Data/RIP-sq_targetRBP_PABPC1/bam_bai_data/ENCFF000VAJ.bam'
#%%
bin_size = 100000
genome_size = 2913022393
bam = pysam.AlignmentFile(bam_path, 'rb')
num_bins = int(np.ceil(genome_size / bin_size))
read_counts = np.zeros(num_bins)

#%%
coverage = np.zeros(num_bins,dtype=int)

#%% some legacy code for binning that doesn't take chromosomes into account
'''
total_count=0
for read in bam.fetch():
    total_count+=1
    if total_count % 100000 == 0:
        print(total_count)
    if not read.is_unmapped:
        bin_index_start = read.reference_start // bin_size
        bin_index_end = read.reference_end // bin_size
        coverage[bin_index_start:bin_index_end] += 1
'''

#%% take chromosomes into account


chromosome_coverage={}
total_count=0
# Process reads
for read in bam.fetch():
    if not read.is_unmapped:
        total_count+=1
        if total_count % 100000 == 0:
            print(total_count)
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

        # Increment coverage for all bins the read spans
        chromosome_coverage[chrom][bin_index_start:bin_index_end + 1] += 1

#%%
coverage_lengths = [len(chromosome_coverage[chrom]) for chrom in chromosome_coverage]

#%% Write everything into a file to save compute time
chrom_data = []
for chrom,coverages in chromosome_coverage.items():
    for index,coverage in enumerate(coverages):
        chrom_data.append((chrom,index,coverage))
chrom_data_df = pd.DataFrame(chrom_data,columns=['chrom','index','coverage'])
chrom_data_df.to_csv('coverage.csv')