import pysam
# idea is to only have nonduplicate reads

def filter_duplicates(input_bam, output_bam):
    input_bamfile = pysam.AlignmentFile(input_bam, "rb")

    output_bamfile = pysam.AlignmentFile(output_bam, "wb", template=input_bamfile)


    for read in input_bamfile:
        if not read.is_duplicate: # is_duplicate is part of the pysam.aligned segment object
            output_bamfile.write(read)

    input_bamfile.close()
    output_bamfile.close()