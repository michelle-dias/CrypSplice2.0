# ExtractJunctions.py

import sys, os
import pandas as pd
from collections import defaultdict
import concurrent.futures as cf
import subprocess

sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities

# Extracting junction locations and counts using pysam #


# Get Junction Counts
## master function for obtaining junctions using pysam
def get_junction_counts(outdir, processors, control, treated, strand, gtf_path, length_cutoff):
    # adjust strand for regtools #
    if strand == 1:
        regtools_strand = 2
    elif strand == 2:
        regtools_strand = 1
    else:
        regtools_strand = 0
    samples=control+treated
    pl=[]
    for sample in samples:
        pl.append([sample,outdir,regtools_strand, length_cutoff])
    with cf.ProcessPoolExecutor(max_workers=int(processors)) as (executor):
        result = list(executor.map(run_regtools, pl))
    if 0 in result:
        return(0)
    else:
        count_list = []
        for sample in samples:
            sample_name=sample.split("/")[-1].replace(".bam","")
            junction_bed=outdir+sample_name+".regtools.juncExtract.bed"
            junction_counts = pd.read_csv(junction_bed, sep="\t", header=None)
            junction_counts = junction_counts.copy()
            junction_counts[[12, 13]] = junction_counts[10].str.split(pat=',', n=1, expand=True)
            # adding/subtracting block sizes due to regtools formatting - check regtools docs for explanation
            junction_counts[1]=junction_counts[1].astype(int)+junction_counts[12].astype(int)
            junction_counts[2]=junction_counts[2].astype(int)+junction_counts[13].astype(int)
            junction_counts = junction_counts[[0,1,2,5,4]]
            junction_counts.columns = ["chrom", "start", "end", "strand", sample_name]
            count_list.append(junction_counts)
            #os.system("rm "+junction_bed)
        merged_juncCounts = pd.concat(count_list).groupby(["chrom", "start", "end", "strand"]).sum().reset_index()
        merged_juncCounts = merged_juncCounts.fillna(0)
        junc_dict = create_dictionary(merged_juncCounts)
        junction_counts = get_junction_origin_counts(merged_juncCounts, control+treated, junc_dict)
        junction_counts['juncID'] = junction_counts.apply(lambda row: ':'.join([str(row['chrom']), str(row['start'])])+'-'+ str(row['end'])+"("+str(row['strand'])+")", axis=1)
        if junction_counts.empty:
            return(0)
        else:
            junction_counts.to_csv(outdir+"JunctionCounts.txt", sep="\t", index=None)
            return(1)





# Run Regtools
## using regtools to extract junctions from each bam file 
def run_regtools(file):
    input_bam=file[0]; output_bed=file[1]+file[0].split("/")[-1].replace(".bam",".regtools.juncExtract.bed")
    strand=str(file[2])
    length_cutoff=str(file[3])
    log = subprocess.run(['regtools','junctions','extract','-s',strand,'-m',length_cutoff,'-o',output_bed,input_bam],stderr=subprocess.DEVNULL,shell=False)
    return(1)


# Create Dictionary 
## create junction counts dictionary from junction counts dataframe to find origin sums
def create_dictionary(juncCount_data):
    junction_dict = {}
    sample_cols = [col for col in juncCount_data.columns if col not in ['chrom', 'start', 'end', 'strand']]
    for column in sample_cols:
        junction_dict[column] = juncCount_data.set_index(['chrom', 'start', 'end', 'strand'])[column].to_dict()
    return(junction_dict)


# Sum Origin Counts
## summing counts for all junctions starting at the same origin point (strand-specific)
def sum_origin_counts(all_regions_dict):
    origin_counts = {
        key: defaultdict(int) for key in all_regions_dict
    }
    for sample, sample_regions in all_regions_dict.items():
        for key, value in sample_regions.items():
            chrom, start, end, strand = key
            if strand == '+':
                new_key = (chrom, start, strand)
            elif strand == '-':
                new_key = (chrom, end, strand)
            else:
                new_key = (chrom, start, end, strand)
            origin_counts[sample][new_key] += value
    return origin_counts


# Getting Junction Origin Counts Per Row
## getting junction origin counts for each individual junction
def get_junction_origin_counts_perRow(row, sample, origin_sums_dict):
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    strand = row['strand']
    if strand == '+':
        key = (chrom, start, '+')
    elif strand == '-':
        key = (chrom, end, '-')
    else:
        key = (chrom, start, end, strand)
    value = origin_sums_dict[sample].get(key, 0)
    return value


# Getting Junction Origin Counts 
## overhead function to get the origin counts   
def get_junction_origin_counts(junctions, sample_paths, all_junctions):
    all_junctions_origin_sums = sum_origin_counts(all_junctions)
    #all_exons_origin_sums = sum_origin_counts(all_exons)
    samples = [sample_path.split("/")[-1].replace(".bam", "") for sample_path in sample_paths]
    for sample in samples:
        junctions[sample+"_originCounts"] = junctions.apply(get_junction_origin_counts_perRow, args=(sample, all_junctions_origin_sums,), axis=1)
    return junctions

