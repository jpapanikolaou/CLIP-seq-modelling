# ClIP-Seq Processor

## Data install details

1. This program is designed to work wth bam files. You can download the bam files from the ENCODE database.
Data that we used is present in ```datainfo.txt``` in the root directory of this repo.
2. Path organization format is as follows:
```ENCODE Data/<experiment_name>/<experiment_type (control or target)>/<bam_indexing_output.bam>```

## Run MACS

1. Run macs with these arguments: ```macs2 callpeak --nomodel -t <input.bam> -g hs -n <output_name>```
2. Convert the ```*_peaks.xls```  file to a CSV manually (this is a one-time thing, and lets one compare the model's output with MACS)

## How to run the pipeline

1. Install BAM files from ENCODE in your preferred data path
2. Install samtools, if you havene't already. Note that you will need to use ```brew install``` if 
you're on a Mac, as the samtools verision on conda is out of date and doesn't have ```markdup``` on the path
3. Do filtering

First, run:
```bash'''
samtools markdup <input file name.bam> <output file name.bam>
```
This allows binning_V3.py to effectively call read.isDuplicate() to filter out duplicates. Results recapitulate MACS


Second, run
```bash'''
samtools index <output file name.bam>
```

This indexes the files and allows bam_file.fetch() to be called, which greatly speeds up the binning process

4. Set up main.py with the appropriate filtered files.

macs_path in ```python main.py``` should be 
```macs_path = "ENCODE Data/<experiment_name>/<experiment_type>/<bam_indexing_output.bam>```, NOT the .bambai extension file

When it prints removed_count, this number should NOT be zero.

5. Run the rest of ```python main.py```, this bit should be fairly self-explanatory :)

## Converting peaks.csv to gene_count.csv

1. Navigate to https://genome.ucsc.edu/ and click on Genome Browser. Select "Go".

2. Now you should be at the UCSC Genome Browser for Humans page. Click Tools near the top then table browser.

3. In select dataset make group "Genes and Gene Predictions" select GENCODE V44 for track and make the table knownCanonical. For region of interest select genome and for output format select "selected fields from primary and related tables".
   
4. After selecting get output, select chrom, chromStart, chromEnd, geneSymbol, description, and hg38 kgXref.

5. Click get output and save.

6. Run ```grep -v -e "non-protein coding RNA" -e "antisense RNA" -e "pseudogene" -e "alt" -e "fix" -e "random"  <outputfile.txt > <outputfilefilt.txt``` in the same working directory as the file you saved

7. Run pipeline
