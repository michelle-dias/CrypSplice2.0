#!/usr/bin/env python3


#!/usr/bin/env python3

import os
import sys
import stat
import shutil
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import pandas.testing as pdt


# Resolve Tool
## function to import extract-as from classification_imports folder
def resolve_tool(name: str) -> str:
    """
    Resolves the path to an executable located in classification_imports folder
    relative to this script.
    """
    here = Path(__file__).resolve().parent
    tool_path = here / "classification_imports" / name

    if not tool_path.exists():
        raise FileNotFoundError(f"Tool {name} not found at {tool_path}")
    
    # Make sure it's executable
    mode = tool_path.stat().st_mode
    if not (mode & stat.S_IXUSR):
        tool_path.chmod(mode | stat.S_IXUSR)

    return str(tool_path)




# Classify Junctions
## master function to run stringtie and classify junctions 
def classify_junctions(bfiles, outDir, processors, gtf, fasta, novel_junctions_path):
    # running stringtie
    # if run_stringtie(bfiles,outDir,processors,gtf) == 1:
    run=1
    if run == 1:
    #     # merging stringtie gtfs
    #     merge_stringtie(bfiles, outDir, gtf, processors)
        #annotate merged gtf
        annotate_gtf(fasta, outDir)
        # add annotations to novel junctions
        novel_junctions = add_annotations(outDir, novel_junctions_path)
        novel_junctions.to_csv(novel_junctions_path, sep="\t", index=None)
        saved_df = pd.read_csv(novel_junctions_path, sep="\t")
        pdt.assert_frame_equal(saved_df, novel_junctions, check_dtype=False, check_like=True)        
        return_val = 1
    else:
        return_val = 0
    return(return_val)



# Run Stringtie
## running stringtie to create de-novo gtfs for each sample
def run_stringtie(bfiles,outDir,processors,gtf):
    try:
        STRINGTIE = resolve_tool("stringtie")
        for bfile in bfiles:
            gtf_out = f"{outDir}/{Path(bfile).stem}.gtf"
            command = f"{STRINGTIE} {bfile} -G {gtf} -p {processors} -o {gtf_out}"
            os.system(command)
    except:
        return(0)
    return(1)



# Merge GTFs with Stringtie
## running stringtie to create de-novo gtfs for each sample
def merge_stringtie(bfiles, outDir, gtf, processors):
    gtf_list_path = outDir+"merge_GTF_list.txt"
    with open(gtf_list_path, "w") as f:
        for bfile in bfiles:
            gtf_path = f"{outDir}/{bfile.split('/')[-1].replace('.bam', '.gtf')}"
            f.write(gtf_path + "\n") 
    command = f"stringtie --merge -G {gtf} -p {processors} -o {outDir}/merged_GTF.gtf {gtf_list_path}"
    os.system(command)
    os.system("rm "+outDir+"merge_GTF_list.txt")



    
# Annotate GTF
## function to annotate the merged GTF
def annotate_gtf(fasta, outDir):
    mergedGTF_path=outDir+'/merged_GTF.gtf'
    cleanGTF_path=outDir+'/merged_GTF.clean.gtf'
    gtf_columns = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    mergedGTF = pd.read_csv(mergedGTF_path, sep="\t", skiprows=2, header=None, names=gtf_columns)
    # cleaning gtf
    cleanGTF = mergedGTF[mergedGTF["feature"]=="exon"]
    cleanGTF = cleanGTF[cleanGTF["seqname"].astype(str).str.startswith("chr")]
    cleanGTF.to_csv(cleanGTF_path, sep="\t", header=False, index=False, quoting=3)
    # annotating the gtf
    EXTRACT_AS = resolve_tool("extract-as")
    command = f"{EXTRACT_AS} {cleanGTF_path} {fasta} > {outDir}/annotate_merged_GTF.gtf"
    os.system(command)
    os.system("rm "+cleanGTF_path+" "+mergedGTF_path)
    

# Add Annotations to Novel Junctions 
def add_annotations(outDir, novel_junctions):
    # read in annotated GTF
    annotated_GTF = outDir+"/annotate_merged_GTF.gtf"
    annotated_GTF = pd.read_csv(annotated_GTF, sep="\t")
    # read in junctions
    novel_junctions = pd.read_csv(novel_junctions, sep="\t")
    novel_junctions["event_type"] = np.nan
    # adding TSS, TTS, and AE annotations
    novel_junctions_TSS = TSS_TTS_AE_eventLabels("TSS", novel_junctions, annotated_GTF)
    novel_junctions_TSS_TTS = TSS_TTS_AE_eventLabels("TTS", novel_junctions_TSS, annotated_GTF)
    novel_junctions_TSS_TTS_AE = TSS_TTS_AE_eventLabels("AE", novel_junctions_TSS_TTS, annotated_GTF)
    # converting all XSKIP_ON to SKIP_ON and deleting duplicates 
    annotated_GTF["event_type"] = annotated_GTF["event_type"].replace("XSKIP_ON", "SKIP_ON")
    annotated_GTF = annotated_GTF.drop_duplicates()
    # adding SKIP and MSKIP annotations
    novel_junctions_TSS_TTS_AE["skipped_exons"] = np.nan
    novel_junctions_TSS_TTS_AE_SKIP = SKIP_eventLabels("SKIP_ON", novel_junctions_TSS_TTS_AE, annotated_GTF)
    novel_junctions_TSS_TTS_AE_SKIP_MSKIP = SKIP_eventLabels("MSKIP_ON", novel_junctions_TSS_TTS_AE_SKIP, annotated_GTF)
    return novel_junctions_TSS_TTS_AE_SKIP_MSKIP
    
    
# Corresponding junctions with TSS, TTS, AE
def TSS_TTS_AE_eventLabels(event, novel_junctions, annotated_GTF):
    # EVENT Dataframe 
    EVENTDat=annotated_GTF[annotated_GTF["event_type"]==event]
    ## EVENT exon start: + junction end = event_start
    novel_junctions = (
        novel_junctions
        .merge(
            EVENTDat.loc[EVENTDat["strand"]=="+", ["chrom","event_start","event_end"]]
            .drop_duplicates()
            .assign(event_type=event+" (exon start)"),
            how="left",
            left_on=["chrom", "end"],
            right_on=["chrom", "event_start"]
        ).assign(
            event_type=lambda df: df[["event_type_x", "event_type_y"]]
            .apply(lambda x: ",".join(x.dropna().unique()), axis=1)
            .replace("", np.nan)
        )
        .drop(columns=["event_start", "event_end"], errors='ignore'))
    ## EVENT exon end: + junction start = event_end
    novel_junctions = (
        novel_junctions
        .merge(
            EVENTDat.loc[EVENTDat["strand"]=="+", ["chrom","event_start","event_end"]]
            .drop_duplicates()
            .assign(event_type=event+" (exon end)"),
            how="left",
            left_on=["chrom", "start"],
            right_on=["chrom", "event_end"]
        )
        .assign(
            event_type=lambda df: df[["event_type_x", "event_type_y"]]
            .apply(lambda x: ",".join(x.dropna().unique()), axis=1)
            .replace("", np.nan)
        )
        .drop(columns=["event_type_x", "event_type_y", "event_start", "event_end"], errors='ignore'))
    ## EVENT exon start: - junction start = event_end
    novel_junctions = (
        novel_junctions
        .merge(
            EVENTDat.loc[EVENTDat["strand"]=="-", ["chrom","event_start","event_end"]]
            .drop_duplicates()
            .assign(event_type=event+" (exon start)"),
            how="left",
            left_on=["chrom", "start"],
            right_on=["chrom", "event_end"]
        )
        .assign(
            event_type=lambda df: df[["event_type_x", "event_type_y"]]
            .apply(lambda x: ",".join(x.dropna().unique()), axis=1)
            .replace("", np.nan)
        )
        .drop(columns=["event_type_x", "event_type_y", "event_start", "event_end"], errors='ignore')
    )
    ## EVENT exon end: - junction end = event_start
    novel_junctions = (
        novel_junctions
        .merge(
            EVENTDat.loc[EVENTDat["strand"]=="-", ["chrom","event_start","event_end"]]
            .drop_duplicates()
            .assign(event_type=event+" (exon end)"),
            how="left",
            left_on=["chrom", "end"],
            right_on=["chrom", "event_start"]
        )
        .assign(
            event_type=lambda df: df[["event_type_x", "event_type_y"]]
            .apply(lambda x: ",".join(x.dropna().unique()), axis=1)
            .replace("", np.nan)
        )
        .drop(columns=["event_type_x", "event_type_y", "event_start", "event_end"], errors='ignore')
    )
    return novel_junctions




# Corresponding junctions with SKIP and MSKIP
def SKIP_eventLabels(event, novel_junctions, annotated_GTF):
    EVENTDat = annotated_GTF[annotated_GTF["event_type"] == event].copy()
    EVENTDat["junc_start"] = pd.to_numeric(
        EVENTDat["event_pattern"].str.split(",").str[0], errors="coerce")
    EVENTDat["junc_end"] = pd.to_numeric(
        EVENTDat["event_pattern"].str.split(",").str[-1], errors="coerce")
    EVENTDat["skipped_exons"] = (
        EVENTDat["event_pattern"].str.split(",")
        .apply(lambda x: ",".join(x[1:-1]) if len(x) > 2 else np.nan))
    EVENTDat["junc_end"] = EVENTDat["junc_end"] - 1 
    novel_junctions = novel_junctions.assign(
        start=pd.to_numeric(novel_junctions["start"], errors="coerce"),
        end=pd.to_numeric(novel_junctions["end"], errors="coerce"),
        chrom=lambda d: d["chrom"].astype(str).str.strip(),
    )
    EVENTDat["chrom"] = EVENTDat["chrom"].astype(str).str.strip()
    out = (
        novel_junctions
        .merge(
            EVENTDat.loc[:, ["chrom","junc_start","junc_end","skipped_exons","event_type"]]
                    .drop_duplicates()
                    .assign(event_type=event.split("_")[0]),
            how="left",
            left_on=["chrom","start","end"],
            right_on=["chrom","junc_start","junc_end"]
        )
        .assign(
            event_type=lambda df: df.filter(like="event_type")
                                     .apply(lambda x: ",".join(pd.unique(x.dropna().astype(str))), axis=1)
                                     .replace("", np.nan),
            skipped_exons=lambda df: df.filter(like="skipped_exons")
                                       .apply(lambda x: ",".join(pd.unique(x.dropna().astype(str))), axis=1)
                                       .replace("", np.nan)
        )
        .pipe(lambda df: df.drop(
            columns=[*df.filter(like="event_type_").columns,
                     *df.filter(like="skipped_exons_").columns,
                     "junc_start","junc_end"],
            errors="ignore"
        ))
    )
    return out




    
# controls="/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Ctrl_1/Ctrl_1.sorted.bam,/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Ctrl_2/Ctrl_2.sorted.bam,/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Ctrl_3/Ctrl_3.sorted.bam"
# treated="/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Cherp_1/Cherp_1.sorted.bam,/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Cherp_2/Cherp_2.sorted.bam,/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Antonia_Data/Alignment//Cherp_3/Cherp_3.sorted.bam"
# gtf="/mnt/localstorage/michelle/data/References/gencode_v28_hg38/gencode.v28.annotation.gtf"
# processors=50
# outDir="/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypSplice_Editing/junction_classification_10_20_25/JuncClass_Output/"
# junctions_path="/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Results/CrypticJunctions/U2Surp/U2Surp.Novel_Junctions.txt"
# fasta="/mnt/localstorage/michelle/data/References/gencode_v28_hg38/GRCh38.p12.genome.fa"
# novel_junctions_path="/mnt/localstorage/michelle/data/Projects/CrypSplice/Venkata_WriteUp/Results/CrypticJunctions/U2Surp/U2Surp.Novel_Junctions.txt"



# import os
# import sys
# import pandas as pd
# import numpy as np 

# control_list = "".join(controls).replace(" ", "").split(",")
# treated_list = "".join(treated).replace(" ", "").split(",")


# bfiles = control_list+treated_list


        
# classify_junction(bfiles, outDir, processors, gtf, novel_junctions_path)
    
#     