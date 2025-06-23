# AnnotateJunctions.py

import os, sys
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities

import numpy as np
import concurrent.futures as cf
import pandas as pd
import pybedtools as pb
import subprocess


# Annotating junctions using regtools keys #

# Annotate Junctions
## master function to add annotations to each junction
def annotate_junctions(junction_counts_path, gtf, fasta, outdir):
    junction_counts = pd.read_csv(junction_counts_path, sep="\t")
    ## scaffolds mess with bed sorting and regtools annotate cuts off at scaffolds without doing chroms
    junction_counts = junction_counts[junction_counts['chrom'].str.contains('chr')]
    junction_df = junction_counts[['chrom', 'start', 'end', 'strand', 'juncID']].copy()
    # reformatting to Bed12 for regtools
    junction_df = bed12_reformat(junction_df)
    junction_bed = pb.BedTool.from_dataframe(junction_df)
    junction_bed_sorted = junction_bed.sort()
    junction_bed_path = outdir+"junctions.tmp.bed"
    junction_bed_sorted.saveas(junction_bed_path)
    regtools_output_path = outdir+"annotations.tmp.txt"
    command = ['regtools', 'junctions', 'annotate', '-S', '-o', regtools_output_path, junction_bed_path, fasta, gtf]
    exit_code = subprocess.run(['regtools','junctions','annotate','-S','-o',regtools_output_path,junction_bed_path,fasta,gtf],stderr=subprocess.DEVNULL,shell=False)
    # 0 means regtools worked
    if exit_code.returncode == 0:
        regtools_annotations = pd.read_csv(regtools_output_path, sep="\t")
        regtools_annotations = regtools_annotations[["name","splice_site","anchor"]]
        regtools_annotations.rename(columns={'name': 'juncID', 'anchor':'annotation'}, inplace=True)
        regtools_annotations = regtools_annotations.drop_duplicates()
        # dropping duplcated IDs
        duplicated_juncIDs = regtools_annotations[regtools_annotations["juncID"].duplicated()]["juncID"].unique()
        duplicated_df = regtools_annotations[regtools_annotations["juncID"].isin(duplicated_juncIDs)]
        retained_dups = duplicated_df[duplicated_df["annotation"]!="N"]
        regtools_annotations_noDups = regtools_annotations[~regtools_annotations['juncID'].isin(duplicated_juncIDs)]
        regtools_annotations = pd.concat([regtools_annotations_noDups, retained_dups])
        regtools_annotations = regtools_annotations.drop_duplicates(subset=['juncID'])
        # merging annotations with junction counts
        merged_df = pd.merge(junction_counts, regtools_annotations[['juncID', 'splice_site', 'annotation']], on='juncID', how='left')
        junction_counts[['splice_site', 'annotation']] = merged_df[['splice_site', 'annotation']]
        merged_df.to_csv(outdir+"Annotated_JunctionCounts.txt", sep="\t", index=None)
        # removing temp files
        os.system("rm "+regtools_output_path)
        os.system("rm "+junction_bed_path)
        return(1)
    else:
        return(0)


# Bed12 Reformatting
## converting junction counts temp bed to bed 12 format for regtools
def bed12_reformat(junction_df):
    for index in [3,4,6,7,9]:
        junction_df['hold_'+str(index)]=0
    junction_df["hold_8"]="255,0,0"
    junction_df["hold_10"]="0,0"
    diff = junction_df["end"] - junction_df["start"]
    diff = "0," + diff.astype(str)
    junction_df["hold_11"]=diff
    junction_df = junction_df[["chrom","start","end","juncID","hold_4","strand","hold_6","hold_7","hold_8","hold_9","hold_10","hold_11"]]
    return(junction_df)