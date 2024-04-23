import pysam
import pandas as pd
import numpy as np

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
               'chrX', 'chrY', 'chrM']
#%%

def bin_reads_to_dataframe(bam_file_path, bin_size):
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    all_bins = {}  # Initialize as dictionary

    total_count = 0
    removed_count = 0
    for read in bam.fetch():
        if (not read.is_unmapped) and (not read.is_duplicate):
            total_count += 1
            if total_count % 100000 == 0:
                print(f"Processed {total_count} reads")

            chrom = read.reference_name
            if chrom not in all_bins:
                chrom_length = bam.lengths[bam.get_tid(chrom)]
                num_bins = int(np.ceil(chrom_length / bin_size))
                all_bins[chrom] = np.zeros(num_bins, dtype=int)

            bin_index_start = read.reference_start // bin_size
            bin_index_end = (read.reference_end - 1) // bin_size
            if bin_index_end >= len(all_bins[chrom]):
                bin_index_end = len(all_bins[chrom]) - 1  # Ensure index is within range

            all_bins[chrom][bin_index_start] += 1
            '''all_bins[chrom][bin_index_start:bin_index_end + 1] += 1 culprit for multiple binnings
                                                                        of the same read'''
        else:
            removed_count +=1

    print("removed count",removed_count)
    bam.close()

    # Create a DataFrame from the results
    rows = []
    for chrom, counts in all_bins.items():
        for index, count in enumerate(counts):
            if count > 0:  # Only include bins that have counts
                start = index * bin_size
                end = start + bin_size
                rows.append({'chromosome': chrom, 'start': start, 'end': end, 'count': count})

    return pd.DataFrame(rows, columns=['chromosome', 'start', 'end', 'count'])


#%%


control_file_path = "ENCODE Data/eClip_control_and_target/Control/ENCFF913CBD.bam"
bam_file = pysam.AlignmentFile(control_file_path, "rb")
# chrom_lengths = extract_regions(bam_file)


#%%
control_binned = bin_reads_to_dataframe(control_file_path,50)

#%%
control_binned.to_csv("ENCODE Data/eClip_control_and_target/control_binned.csv", index=False)

