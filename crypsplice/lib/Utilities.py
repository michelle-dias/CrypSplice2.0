# Utility.py 


import pandas as pd

# Utility functions that are necessary throughout the modules #


# Read GTF File
## reads in the gencode GTF file 
def read_gtf_file(file_path):
    col_names = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"] 
    dtype_dict = {
        "chrom": str,   # Chrom column as string
        "start": int,   # Start column as integer
        "end": int      # End column as integer
    }
    gtf = pd.read_csv(file_path, sep="\t", comment="#", header=None, names=col_names, dtype=dtype_dict)
    return gtf

