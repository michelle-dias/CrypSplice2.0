# DifferentialUsage.py

import sys, os
import subprocess
import pandas as pd


# Testing for differentially used junctions using the countdata R package #


# Run bbTest
## calling the countdata R function to calculate the significance
def run_bbTest(junctions, control_num, treated_num, processors, test_type):
    Rscript_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "RunBBTest.R")   
    arguments = [junctions, str(control_num), str(treated_num), str(processors), str(test_type)]
    command = ["Rscript", Rscript_path] + arguments
    result = subprocess.check_output(command, universal_newlines=True)
    return 1 if "Done" in result else 0


# Calculate Junction Strength
## quantifying junction strength by dividing the average junction counts per condition over the average total counts
## ex: (C1 junction average / (C1 origin counts average)
def calc_junction_strength(sjdfBB_path, nc, nt, start_idx):
    sjdfBB = pd.read_csv(sjdfBB_path, sep="\t")
    ## get junction average counts ##
    sjdfBB['C1_Cavg'] = sjdfBB[sjdfBB.columns[start_idx:start_idx+nc]].mean(axis=1)
    sjdfBB['C2_Cavg'] = sjdfBB[sjdfBB.columns[start_idx+nc:start_idx+nc+nt]].mean(axis=1)
    sjdfBB['C1_Tavg'] = sjdfBB[sjdfBB.columns[start_idx+nc+nt:start_idx+nc+nt+nc]].mean(axis=1)
    sjdfBB['C2_Tavg'] = sjdfBB[sjdfBB.columns[start_idx+nc+nt+nc:start_idx+nc+nt+nc+nt]].mean(axis=1)
    ## get junction strength columns and diff ##
    sjdfBB['C1_CS'] = sjdfBB['C1_Cavg']/sjdfBB['C1_Tavg']
    sjdfBB['C2_CS'] = sjdfBB['C2_Cavg']/sjdfBB['C2_Tavg']
    sjdfBB['C1_CS'] = sjdfBB['C1_CS'].fillna(0)
    sjdfBB['C2_CS'] = sjdfBB['C2_CS'].fillna(0)
    sjdfBB['CS_diff']=sjdfBB['C2_CS']-sjdfBB['C1_CS']
    if 'C1_CS' in sjdfBB.columns and 'C2_CS' in sjdfBB.columns and not sjdfBB[['C1_CS', 'C2_CS']].empty:
        sjdfBB.to_csv(sjdfBB_path, sep="\t", index=None)
        return(1)
    else:
        return(0)




