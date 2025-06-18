# ExtractJunctions.py

import sys, os
import pandas as pd
from collections import defaultdict
import concurrent.futures as cf
import subprocess
from concurrent.futures import ProcessPoolExecutor
import gc

sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities

# Extracting junction locations and counts using regtools #


# Get Junction Counts
## master function for obtaining junctions using pysam
def get_junction_counts(outdir, processors, control, treated, strand, gtf_path, length_cutoff):
    # adjust strand for regtools #
    if strand == 1:
        regtools_strand = "FR"
    elif strand == 2:
        regtools_strand = "RF"
    else:
        regtools_strand = "XS"
    samples=control+treated
    pl=[]
    for sample in samples:
        pl.append([sample,outdir,regtools_strand, length_cutoff])
    with cf.ProcessPoolExecutor(max_workers=int(processors)) as (executor):
        result = list(executor.map(run_regtools, pl))
    if 0 in result:
        return(0)
    else:
        params = [(sample, outdir) for sample in samples]
        with ProcessPoolExecutor(max_workers=processors) as executor:
            results = list(executor.map(adjust_juctions_wrapper, params))
        if len(results)==len(samples):
            return(1)
        else:
            return(0)



    
# Run Regtools
## using regtools to extract junctions from each bam file 
def run_regtools(file):
    input_bam=file[0]; output_bed=file[1]+file[0].split("/")[-1].replace(".bam",".regtools.juncExtract.bed")
    strand=str(file[2])
    length_cutoff=str(file[3])
    command = ['regtools', 'junctions', 'extract','-s', strand,'-m', length_cutoff,'-o', output_bed,input_bam]
    command_string = ' '.join(command)
    print(command_string)
    log = subprocess.run(command, stderr=subprocess.DEVNULL, shell=False)
    return(1)
 

# Adjust Junctions Wrapper
## wrapper to pass multiple parameters to parallel processor
def adjust_juctions_wrapper(params):
    return adjust_juctions(*params)


# Adjust Junctions
## adjusting junction coordinates based on regtools protocol
def adjust_juctions(sample, outdir):
    sample_name = os.path.basename(sample).replace(".bam", "")
    junction_bed_path=outdir+sample_name+".regtools.juncExtract.bed"
    adjusted_junction_bed_path=outdir+sample_name+".ExtractJunctions.bed"
    junction_bed= pd.read_csv(junction_bed_path, sep="\t", header=None)
    junction_bed[[12, 13]] = junction_bed[10].str.split(pat=',', n=1, expand=True)
    # Adding/subtracting block sizes
    junction_bed[1] = junction_bed[1].astype(int) + junction_bed[12].astype(int)
    junction_bed[2] = junction_bed[2].astype(int) - junction_bed[13].astype(int)
    junction_bed = junction_bed[[0, 1, 2, 5, 4]]
    junction_bed.columns = ["chrom", "start", "end", "strand", sample_name]
    junction_bed = junction_bed.groupby(["chrom", "start", "end", "strand"]).sum().reset_index()
    # finding origin counts
    junction_counts = get_origin_counts(junction_bed)
    # adding juncID label
    junction_counts['juncID'] = (junction_counts['chrom'].astype(str) + ':' + junction_counts['start'].astype(str) + '-' + junction_counts['end'].astype(str) + '(' + junction_counts['strand'].astype(str) + ')')
    junction_counts.to_csv(adjusted_junction_bed_path, sep="\t", index=None)
    # removing original regtools file
    os.system("rm "+junction_bed_path)
    return(adjusted_junction_bed_path)
    
    
# Getting Junction Origin Counts 
# getting junction origin counts for each file using pandas merge and subset functions
def get_origin_counts(juncCounts):
    sample_columns = [col for col in juncCounts.columns if col not in ['chrom', 'start', 'end', 'strand']]
    # summing origins [chrom, start] for + strand
    counts_by_start = juncCounts[juncCounts['strand'] == '+'].groupby(['chrom', 'start', 'strand'])[sample_columns].sum()
    counts_by_start.rename(columns={col: f'{col}_originCounts' for col in counts_by_start.columns if col not in ['chrom', 'start', 'end', 'strand']}, inplace=True)
    counts_by_start = counts_by_start.reset_index()
    # summing origins [chrom, end] for - strand
    counts_by_end = juncCounts[juncCounts['strand'] == '-'].groupby(['chrom', 'end','strand'])[sample_columns].sum()
    counts_by_end.rename(columns={col: f'{col}_originCounts' for col in counts_by_end.columns if col not in ['chrom', 'start', 'end', 'strand']}, inplace=True)
    counts_by_end = counts_by_end.reset_index()
    # changing ending of columns with jucntion counts to "juncCounts" 
    juncCounts.rename(columns={col: f'{col}_juncCounts' for col in juncCounts.columns if col not in ['chrom', 'start', 'end', 'strand']}, inplace=True)
    # subsetting original data for merge with origin dataframes
    pos_data = juncCounts[juncCounts["strand"]=="+"]
    neg_data = juncCounts[juncCounts["strand"]=="-"]
    # merging junction data with origin data
    pos_data_wOrigin = pos_data.merge(counts_by_start, on=['chrom', 'start', 'strand'], how='left')
    neg_data_wOrigin = neg_data.merge(counts_by_end, on=['chrom', 'end', 'strand'], how='left')
    fullData_wOrigin = pd.concat([pos_data_wOrigin, neg_data_wOrigin], ignore_index=True)
    fullData_wOrigin = fullData_wOrigin.sort_values(by=['chrom', 'start', 'end', 'strand'])
    # Clean up memory
    del counts_by_start, counts_by_end, pos_data, neg_data, pos_data_wOrigin, neg_data_wOrigin
    gc.collect()
    return(fullData_wOrigin)




# Stitch Extract Junctions
## stitch together the singular output files from extract junctions into JunctionCounts.txt
def stitch_extractJunctions(control_bams, treated_bams, outdir):
    # get correct order for extractJunc file to be joined - for following steps
    control_names = [os.path.basename(path).replace(".bam", "") for path in control_bams]
    control_beds = [outdir+sample+".ExtractJunctions.bed" for sample in control_names]
    treated_names = [os.path.basename(path).replace(".bam", "") for path in treated_bams]
    treated_beds = [outdir+sample+".ExtractJunctions.bed" for sample in treated_names]
    sample_beds = control_beds+treated_beds
    merged_df = pd.DataFrame()
    for sample_bed in sample_beds:
        juncCounts = pd.read_csv(sample_bed, delimiter='\t') 
        if merged_df.empty:
            merged_df = juncCounts
        else:
            merged_df = pd.merge(merged_df, juncCounts, on=['chrom', 'start', 'end', 'strand', 'juncID'], how='outer')
    # Reindex DataFrame with the juncCounts cols first then originCounts
    junc_counts_cols = [col for col in merged_df.columns if '_juncCounts' in col]
    origin_counts_cols = [col for col in merged_df.columns if '_originCounts' in col]
    reordered_columns = ['chrom', 'start', 'end', 'strand'] + junc_counts_cols + origin_counts_cols + ['juncID']
    merged_df = merged_df[reordered_columns]
    # Fill na with 0
    merged_df=merged_df.fillna(0)        
    # Save the merged DataFrame to a new file
    if not merged_df.empty:
        merged_df.to_csv(outdir+'JunctionCounts.txt', sep="\t", index=False)
        [os.remove(f) if os.path.exists(f) else None for f in sample_beds]
        return(1)
    else:
        return(0)
