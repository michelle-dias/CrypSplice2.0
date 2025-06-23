# FilterJunctions.py
import numpy as np
import glob
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

## NOTE: general junction length and read cutoffs are in the ExtractJunctions module - for speed of that module


    

# Cryptic Load Filter
## sets 'juncCounts' to 0 if 'originCounts'  are below a cutoff 
def juncCut_filter(outdir, processors, junction_cutoff, origin_cutoff):
    junction_bedFiles = glob.glob(outdir+"*ExtractJunctions.bed")
    params = [(junction_bedFile_path, junction_cutoff, origin_cutoff) for junction_bedFile_path in junction_bedFiles]
    with ProcessPoolExecutor(max_workers=processors) as executor:
        filtered_bedFiles = list(executor.map(juncCut_filter_wrapper, params))
    if len(filtered_bedFiles)==len(junction_bedFiles):
        return(1)
    else:
        return(0)

# Cryptic Load Filter Wrapper
## wrapper to pass multiple params to parallelization
def juncCut_filter_wrapper(params):
    return juncCut_filter_parallelized(*params)

# Cryptic Load Filter Parallelization
## parallelization of load filter
def juncCut_filter_parallelized(junction_bedFile_path, junction_cutoff, origin_cutoff):
    junction_bedFile = pd.read_csv(junction_bedFile_path, sep="\t")
    # dropping rows where origin col is below origin_cutoff or junction col is below junction_cutoff
    junction_bedFile = junction_bedFile[
        (junction_bedFile[[col for col in junction_bedFile.columns if col.endswith('juncCounts')]].gt(junction_cutoff).all(axis=1)) &
        (junction_bedFile[[col for col in junction_bedFile.columns if col.endswith('originCounts')]].gt(origin_cutoff).all(axis=1))]
    junction_bedFile.to_csv(junction_bedFile_path, sep="\t", index=None) 
    return(junction_bedFile_path)


# PoverAM Filter (Cryptic Junctions)
## retaining rows with greater than 'M' reads in all condition specific samples 
## retaining rows with greater than 'A' reads in 'P' proportion of samples 
def PoverAM_filter(sjdf, nc, nt, P, A, M):
    control_df = sjdf.iloc[:, 4:4+nc]
    treated_df = sjdf.iloc[:, 4+nc:4+nc+nt]
    ## retaining rows with at least 'M' reads in all rows within controls and treated ##
    control_df = control_df[(control_df>= M).sum(axis=1) == nc]
    treated_df = treated_df[(treated_df>= M).sum(axis=1) == nt]
    ## retaining rows with greater than 'A' reads in 'P' proportion of samples greater than ##
    ## decimal P (1.5) rounded up automatically -> P=0.5 and there are 3 samples, 2 will have to be above cutoff
    control_df=control_df[((control_df >= A).sum(axis=1) >= nc*P)]
    treated_df=treated_df[((treated_df >= A).sum(axis=1) >= nt*P)]
    control_indices = control_df.index.tolist()
    treated_indices = treated_df.index.tolist()
    superset_indices = set(control_indices) | set(treated_indices)
    sjdf_filtered = sjdf.loc[list(superset_indices)]
    sjdf_filtered=sjdf_filtered.reset_index(drop=True)
    return(sjdf_filtered)


# PoverAM Filter (Cryptic Load)
## retaining rows with greater than 'A' reads in 'P' proportion of samples 
def PoverA_filter_CL(sjdf, nc, nt, P, A):
    control_df = sjdf.iloc[:, 1:1+nc]
    treated_df = sjdf.iloc[:, 1+nc:1+nc+nt]
    ## retaining rows with greater than 'A' reads in 'P' proportion of samples greater than ##
    ## decimal P (1.5) rounded up automatically -> P=0.5 and there are 3 samples, 2 will have to be above cutoff
    control_df=control_df[((control_df >= A).sum(axis=1) >= nc*P)]
    treated_df=treated_df[((treated_df >= A).sum(axis=1) >= nt*P)]
    control_indices = control_df.index.tolist()
    treated_indices = treated_df.index.tolist()
    superset_indices = set(control_indices) | set(treated_indices)
    sjdf_filtered = sjdf.loc[list(superset_indices)]
    sjdf_filtered=sjdf_filtered.reset_index(drop=True)
    return(sjdf_filtered)


# Gene Filter
## removing junctions spanning multiple genes or junctions not attributed to a gene
def gene_filter(sjdf):
    sjdf['gene_id'] = sjdf['gene_id'].astype(str)
    sjdf['gene_id'].replace('nan', np.nan, inplace=True)
    sjdf = sjdf.dropna(subset=['gene_id'])
    sjdf = sjdf[~sjdf['gene_id'].str.contains(',')]
    return(sjdf)



