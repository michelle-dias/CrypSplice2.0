# FilterJunctions.py
import numpy as np

## NOTE: general junction length and read cutoffs are in the ExtractJunctions module - for speed of that module



# Load Quality Filter
## sets 'juncCounts' to 0 if both 'originCounts' and 'juncCounts' are below a cutoff 
## want to maintain specificity for load clustering by not dropping everything
def load_quality_filter(junction_counts, read_cutoff):
    #The code sets 'juncCounts' to 0 if both 'originCounts' and 'juncCounts' are below a cutoff.
    lowQual_junctions = junction_counts.filter(like='originCounts').lt(read_cutoff).all(axis=1)
    junction_counts.loc[lowQual_junctions, junction_counts.columns[junction_counts.columns.str.endswith('juncCounts')]] = 0
    return(junction_counts)



# PoverA Filter
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


# Gene Filter
## removing junctions spanning multiple genes or junctions not attributed to a gene
def gene_filter(sjdf):
    sjdf['gene'] = sjdf['gene'].astype(str)
    sjdf['gene'].replace('nan', np.nan, inplace=True)
    sjdf = sjdf.dropna(subset=['gene'])
    sjdf = sjdf[~sjdf['gene'].str.contains(',')]
    return(sjdf)



