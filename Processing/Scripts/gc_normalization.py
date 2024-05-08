#%%
import Bio.SeqUtils
import pysam
from Bio.SeqUtils import GC_skew
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

#%%


def bin_reads_and_calculate_gc_skew(bam_file_path, bin_size):
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    all_bins = {}
    gc_skews = {}  # Store concatenated sequences for GC skew calculation

    total_count = 0
    for read in bam.fetch():
        if not read.is_unmapped and not read.is_duplicate and not read.is_secondary:
            total_count += 1
            if total_count % 100000 == 0:
                print(f"Processed {total_count} reads")
            chrom = read.reference_name
            if chrom not in all_bins:
                chrom_length = bam.lengths[bam.get_tid(chrom)]
                num_bins = int(np.ceil(chrom_length / bin_size))
                all_bins[chrom] = np.zeros(num_bins, dtype=int)
                gc_skews[chrom] = [''] * num_bins  # Initialize with empty strings

            bin_index = read.reference_start // bin_size
            if bin_index < len(all_bins[chrom]):
                all_bins[chrom][bin_index] += 1
                sequence = read.query_sequence
                gc_skews[chrom][bin_index] += sequence

    bam.close()

    # Calculate GC skew and create DataFrame
    rows = []
    for chrom, bins in all_bins.items():
        for index, count in enumerate(bins):
            if count > 0:  # Only include bins that have counts
                gc_sequences = gc_skews[chrom][index:index + bin_size]  # Extract reads within the window
                gc_contents = [Bio.SeqUtils.gc_fraction(seq) for seq in gc_sequences]  # Calculate GC content for each read
                average_gc_content = sum(gc_contents) / len(gc_contents)  # Calculate average GC content
                start = index * bin_size
                end = start + bin_size
                rows.append({'chromosome': chrom, 'start': start, 'end': end,
                             'count': count, 'gc_content': average_gc_content})

        return pd.DataFrame(rows, columns=['chromosome', 'start', 'end', 'count', 'gc_skew'])

#%%

def has_sufficient_data(start, end, df):
    # Check if there are neighboring bins with non-zero counts
    neighboring_bins = df[(df['start'] >= start - 1) & (df['end'] <= end + 1) & (df['count'] > 0)]
    return len(neighboring_bins) > 0

def perform_normalization(count, gc_skew):
    # Simple normalization based on the GC skew
    if gc_skew==0:
        return count
    normalized_count = count / gc_skew
    return normalized_count

def normalize_calculated_skew(df):
    normalized_counts = []
    for index,row in df.iterrows():
        start = row['start']
        end = row['end']
        gc_content = row['gc_skew']
        count = row['count']

        if has_sufficient_data(start, end, df):
            normalized_count = perform_normalization(count, gc_content)
            normalized_counts.append(normalized_count)
        else:
            normalized_counts.append(count)
    df['normalized_count']=normalized_counts

#%%

control_file_path = "../../ENCODE Data/eClip_control_and_target/Control/ControlFiltering/control_filtered_reads.bam"
calculated_skew = bin_reads_and_calculate_gc_skew(control_file_path, 200)
#%%
# normalized_data = normalize_dynamic_windows(calculated_skew, 10000)
#%%
