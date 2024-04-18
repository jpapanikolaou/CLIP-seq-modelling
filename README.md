# ClIP-Seq Processor

## How to run the pipeline

1. Install BAM files from ENCODE in your preferred data path
2. Run filtering.py to filter the reads (remove duplicates)
3. Use do_binning.py to bin the reads. This will take ~5 min for 42M reads.
4. do_binning.py will write a CSV to a location that you specify with reads
    so you don't have to run things multiple times
5. Run model_calling.py to execute the model script, which will create
    two dataframes with the results of the model - once for control, once
    for target
6. The 'meat' of this project is in 'Project/Processing/Scripts'. Everything
    in the pipeline I described to you calls from there