# ExtractJunctions.py

import sys, os
import pandas as pd
from collections import defaultdict
import concurrent.futures as cf
import pysam

sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities

# Extracting junction locations and counts using pysam #


# Get Junction Counts
## master function for obtaining junctions using pysam
def get_junction_counts(outDir, processors, control, treated, strand, gtf_path, read_cutoff, length_cutoff):
    ## reading in gtf data and getting list of uniue chromosomes
    gtf_data = Utilities.read_gtf_file(gtf_path)
    chromosomes = gtf_data["chrom"].unique()
    #all_junctions, all_exons = parallelize_extract_counts(control+treated, strand, chromosomes, processors)
    #all_junctions = parallelize_extract_counts(control+treated, strand, chromosomes, processors)
    all_junctions = parallelize_extract_counts(control+treated, strand, chromosomes, processors)
    junctions = parallelize_merge_filter_counts(all_junctions, control, treated, chromosomes, read_cutoff, processors)
    # drop junctions that are on both strands (technical mapping issues)
    junctions = junctions[~junctions.duplicated(subset=["chrom","start","end"], keep=False)]
    # retain junctions with length at least 25 bp 
    junctions = junctions[junctions["end"]-junctions["start"]>=length_cutoff]
    junctions = get_junction_origin_counts(junctions, control+treated, all_junctions)
    junctions['juncID'] = junctions.apply(lambda row: ':'.join([str(row['chrom']), str(row['start'])]) + '-' + str(row['end']), axis=1)
    junctions.to_csv(outDir+".JunctionCounts.txt", sep="\t", index=None)
    return 0 if junctions.empty else 1


# Parallelize Extract Counts
## parallelization of extract_counts for increased efficiency - split by chromosome
def parallelize_extract_counts(samples, strand, chromosomes, processors):
    all_junctions = defaultdict(int)
    #all_exons = defaultdict(int)
    parameter_list=[]
    for sample_path in samples:    
        sampleID = sample_path.split("/")[-1].replace(".bam","")
        all_junctions[sampleID] = defaultdict(int)
        #all_exons[sampleID] = defaultdict(int)
        for chromosome in chromosomes:
            parameter_list.append([sample_path, sampleID, strand, chromosome])
    with cf.ProcessPoolExecutor(max_workers=int(processors)) as executor:
        #for sampleID, junctions, exons in executor.map(extract_counts, parameter_list):
        for sampleID, junctions in executor.map(extract_counts, parameter_list):
            all_junctions[sampleID].update(junctions)
            #all_exons[sampleID].update(exons)
    return all_junctions#, all_exons


# Extract Counts
## extracting counts for each junction using pysam - also defining strand information 
def extract_counts(parameter_list):
    bam_file_path = parameter_list[0]
    sampleID = parameter_list[1]
    strand = parameter_list[2]
    chromosome = parameter_list[3]
    junctions = defaultdict(int)
    #exons = defaultdict(int)
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch(chromosome):
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
            #if read.is_unmapped:
                continue
            #chrom = bam_file.get_reference_name(read.reference_id)
            blocks = read.get_blocks()
            if strand == 0:
                read_strand = 0
            if strand == 2:
                if read.is_read1 and read.is_reverse:
                    read_strand = "+"
                elif read.is_read1:
                    read_strand = "-"
                elif read.is_read2 and read.is_reverse:
                    read_strand = "-"
                elif read.is_read2:
                    read_strand = "+"
            #for block in blocks:
            #    exons[(chrom, block[0], block[1], read_strand)] += 1
            if "N" in read.cigarstring:
                for i in range(1, len(blocks)):
                    junction_start = blocks[i - 1][1]
                    junction_end = blocks[i][0]
                    junctions[(chromosome, junction_start, junction_end, read_strand)] += 1
    return sampleID, junctions#, exons


# Parallelize Merge and Filter Counts 
## parallelizing filtering and merging of junctions counts extracted from each chromosome
def parallelize_merge_filter_counts(all_counts, control, treated, chromosomes, read_cutoff, processors):
    merged_filtered_counts = defaultdict(int)
    parameter_list = []
    for chromosome in chromosomes:
        all_counts_perChrom = {key: {sub_key: sub_value for sub_key, sub_value in value.items() if sub_key[0] == chromosome} for key, value in all_counts.items()}
        parameter_list.append([all_counts_perChrom, control, treated, chromosome, read_cutoff])
    with cf.ProcessPoolExecutor(max_workers=int(processors)) as executor:
        for filtered_counts_perChrom in executor.map(merge_filter_counts, parameter_list):
            merged_filtered_counts.update(filtered_counts_perChrom)
    sampleIDs = [sample.split("/")[-1].replace(".bam","") for sample in control+treated]
    merged_filtered_counts = pd.DataFrame(merged_filtered_counts).T.reset_index()
    #merged_filtered_counts.columns = ["chrom", "start", "end", "strand"] + sampleIDs
    merged_filtered_counts.columns = ["chrom", "start", "end", "strand"] + [sampleID + "_juncCounts" for sampleID in sampleIDs]

    return merged_filtered_counts


# Merge and Filter Counts 
## merging and filtering the junction counts extracted for each chromosome and sample
def merge_filter_counts(parameter_list):
    all_counts_perChrom = parameter_list[0]
    control = parameter_list[1]
    treated = parameter_list[2]
    chromosome = parameter_list[3]
    read_cutoff = parameter_list[4]
    merged_counts_perChrom = defaultdict(list)
    keys = set().union(*all_counts_perChrom.values())
    sampleIDs = [sample.split("/")[-1].replace(".bam","") for sample in control+treated]
    for sampleID in sampleIDs:
        for key in keys:
            value = all_counts_perChrom[sampleID].get(key, 0)
            merged_counts_perChrom[key].append(value)
    filtered_counts_perChrom = {key: value for key, value in merged_counts_perChrom.items() if all(x > read_cutoff for x in value[len(control):]) or all(x > read_cutoff for x in value[:len(control)])}
    return filtered_counts_perChrom


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
