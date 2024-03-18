# AnnotateJunctions.py

import os, sys
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities

import numpy as np
import concurrent.futures as cf
import pandas as pd

# Annotating junctions using regtools keys #

# Annotate Junctions
## master function to add annotations to each junction
def annotate_junctions(junction_counts, gtf_path, outDir, strand, processors):
    gtf_data = Utilities.read_gtf_file(gtf_path)
    exons_gtf = gtf_data[gtf_data["feature"]=="exon"]
    transcript_exons = build_transcript_exons(exons_gtf)
    all_known_junctions = set()
    for transcript_id, exons in transcript_exons.items():
        for i in range(len(exons) - 1):
            chrom, start, end, strand = exons[i][0], exons[i][2], exons[i + 1][1], exons[i][3]
            all_known_junctions.add((chrom, start, end, strand))
    global all_known_junctions_array
    all_known_junctions_array = np.array(list(all_known_junctions))
    split_range = int(len(junction_counts) / int(processors))
    split_junctions_df = np.array_split(junction_counts, range(split_range, len(junction_counts), split_range))
    results = []
    with cf.ProcessPoolExecutor(max_workers=int(processors)) as executor:
        for array in executor.map(create_annotated_junction_arrays, split_junctions_df):
            if len(results) == 0:
                results = array
            else:
                results = np.concatenate((results, array), axis=0)
    junction_counts_annotated = adding_annotations(junction_counts, results)
    junction_counts_annotated.to_csv(outDir+".Annotated_JunctionCounts.txt", sep="\t", index=None)
    return 0 if junction_counts_annotated.empty else 1


# Build Transcripts Exons
## building the exon dictionary from gtf information
def build_transcript_exons(gtf_data):
    transcript_exons = {}
    attributes_col = "attributes" if "attributes" in gtf_data.columns else 8
    for _, row in gtf_data[gtf_data["feature"] == "exon"].iterrows():
        chrom, start, end, strand, transcript_id = row["chrom"], int(row["start"]) - 1, int(row["end"]), row["strand"], row[attributes_col].split("transcript_id ")[1].split(';')[0].replace('"', '')
        if transcript_id not in transcript_exons:
            transcript_exons[transcript_id] = []
        transcript_exons[transcript_id].append((chrom, start, end, strand))
    for transcript_id in transcript_exons:
        transcript_exons[transcript_id].sort(key=lambda x: (x[1], x[2]))
    return transcript_exons


# Find Annotated Junctions
## finding the junctions that correspond to annotated exon ends
def find_annotated_junctions(junctions_array_elem):
    chrom_annotated = all_known_junctions_array[:, 0] == junctions_array_elem[0]
    start_annotated = all_known_junctions_array[:, 1] == str(junctions_array_elem[1])
    end_annotated = all_known_junctions_array[:, 2] == str(junctions_array_elem[2])
    strand_annotated = all_known_junctions_array[:, 3] == junctions_array_elem[3]
    annotated_junction_bool = np.any(np.all([chrom_annotated, start_annotated, end_annotated, strand_annotated], axis=0))
    annotated_start_bool = np.any(np.all([chrom_annotated, start_annotated, strand_annotated], axis=0))
    annotated_end_bool = np.any(np.all([chrom_annotated, end_annotated, strand_annotated], axis=0))
    junctions_array_elem_annotated = np.append(junctions_array_elem, [annotated_junction_bool, annotated_start_bool, annotated_end_bool])
    return junctions_array_elem_annotated


# Create Annotated Junction Arrays
## create junction arrays to send into parallelization
def create_annotated_junction_arrays(split_junctions_df):
    split_junctions_array = np.array(split_junctions_df) 
    annotated_arrays = np.apply_along_axis(find_annotated_junctions, axis=1, arr=split_junctions_array)
    return annotated_arrays


# Adding Annotations
## adding annotations per regtools keys (DA, D, A, NDA, N)
def adding_annotations(junction_counts, results):
    junction_counts_headers = junction_counts.columns.tolist()
    junction_counts_annotated = pd.DataFrame(results)
    junction_counts_headers = junction_counts_headers + ["annotated_junction", "annotated_start", "annotated_end"]
    junction_counts_annotated.columns = junction_counts_headers
    junction_counts_annotated.loc[junction_counts_annotated['annotated_junction'] == True, 'annotation'] = 'DA'
    junction_counts_annotated.loc[((junction_counts_annotated['annotated_junction'] == False) & (junction_counts_annotated['annotated_start'] == True) & (junction_counts_annotated["annotated_end"] == True)), 'annotation'] = 'NDA'
    junction_counts_annotated.loc[((junction_counts_annotated['annotated_junction'] == False) & (junction_counts_annotated['annotated_start'] == True) & (junction_counts_annotated["annotated_end"] == False) & (junction_counts_annotated["strand"] == "+")), 'annotation'] = 'D'
    junction_counts_annotated.loc[((junction_counts_annotated['annotated_junction'] == False) & (junction_counts_annotated['annotated_start'] == False) & (junction_counts_annotated["annotated_end"] == True) & (junction_counts_annotated["strand"] == "-")), 'annotation'] = 'D'
    junction_counts_annotated.loc[((junction_counts_annotated['annotated_junction'] == False) & (junction_counts_annotated['annotated_start'] == False) & (junction_counts_annotated["annotated_end"] == True) & (junction_counts_annotated["strand"] == "+")), 'annotation'] = 'A'
    junction_counts_annotated.loc[((junction_counts_annotated['annotated_junction'] == False) & (junction_counts_annotated['annotated_start'] == True) & (junction_counts_annotated["annotated_end"] == False) & (junction_counts_annotated["strand"] == "-")), 'annotation'] = 'A'
    junction_counts_annotated['annotation'].fillna('N', inplace=True)
    junction_counts_annotated.drop(columns=["annotated_junction", "annotated_end", "annotated_start"], inplace=True)
    return junction_counts_annotated

