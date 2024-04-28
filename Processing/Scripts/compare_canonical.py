import pyranges as pr
import pandas as pd

#%%
def compare_canonical(genes_path,peaks_path,is_macs=False):
    genes_df = pd.read_csv(genes_path,sep='\t')
    peaks_df = pd.read_csv(peaks_path)
    if is_macs:
        peaks_df = peaks_df.rename(columns={"chr":"Chromosome","start":"Start","end":"End"})
    else:
        peaks_df = peaks_df[peaks_df['significant']==True]
        peaks_df = peaks_df.rename(columns={"chromosome":"Chromosome","start":"Start","end":"End"})
    genes_df = genes_df.rename(columns={"#hg38.knownCanonical.chrom":"Chromosome",
                                        "hg38.knownCanonical.chromStart":"Start",
                                        "hg38.knownCanonical.chromEnd":"End"})
    gr_genes = pr.PyRanges(genes_df)
    gr_peaks = pr.PyRanges(peaks_df)
    overlap = gr_genes.join(gr_peaks)
    overlap_df = overlap.df
    return overlap_df


# #%%
#
# genes_df = pd.read_csv("../GenomeBrowserData/hgTableFilt.txt", sep="\t")
# target_peaks_df = pd.read_csv("../ENCODE Data/eClip_control_and_target/our_target_peaks.csv")
# control_peaks_df = pd.read_csv("../ENCODE Data/eClip_control_and_target/our_control_peaks.csv")
#
# target_peaks_df = target_peaks_df[target_peaks_df['significant']==True]
# control_peaks_df = control_peaks_df[control_peaks_df['significant']==True]
#
# #%%
# genes_df = genes_df.rename(columns={"#hg38.knownCanonical.chrom":"Chromosome",
#                                     "hg38.knownCanonical.chromStart":"Start",
#                                     "hg38.knownCanonical.chromEnd":"End"},
#                            )
# control_peaks_df = control_peaks_df.rename(columns={"chromosome":"Chromosome","start":"Start","end":"End"})
# target_peaks_df = target_peaks_df.rename(columns={"chromosome":"Chromosome","start":"Start","end":"End"})
#
# #%%
#
# gr_genes = pr.PyRanges(genes_df)
# gr_target = pr.PyRanges(target_peaks_df)
# gr_control = pr.PyRanges(control_peaks_df)
#
#
# #%%
# control_overlap = gr_genes.join(gr_control)
# target_overlap = gr_genes.join(gr_target)
# #%%
# control_overlap_df = control_overlap.df
# target_overlap_df = target_overlap.df
# #%%
