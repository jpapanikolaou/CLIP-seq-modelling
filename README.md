# ClIP-Seq Processor

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