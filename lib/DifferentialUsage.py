# DifferentialUsage.py

import sys, os
import subprocess
import pandas as pd


# Testing for differentially used junctions using the countdata R package #


# Run bbTest
## calling the countdata R function to calculate the significance
def run_bbTest(junctions, control_num, treated_num, processors, test_type):
	Rscript_path = "/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib/RunBBTest.R"
	arguments = [junctions, str(control_num), str(treated_num), str(processors), str(test_type)]
	command = ["Rscript", Rscript_path] + arguments
	result = subprocess.check_output(command, universal_newlines=True)
	return 1 if "Done" in result else 0


# Calculate Junction Strength
## quantifying junction strength by dividing the average junction counts per condition over the average total counts
## ex: (C1 junction average / (C1 origin counts average)
def calc_junction_strength(sjdfBB_path, nc, nt):
    sjdfBB = pd.read_csv(sjdfBB_path, sep="\t")
    ## get junction average counts ##
    sjdfBB['C1_Javg'] = sjdfBB[sjdfBB.columns[4:4+nc]].mean(axis=1)
    sjdfBB['C2_Javg'] = sjdfBB[sjdfBB.columns[4+nc:4+nc+nt]].mean(axis=1)
    sjdfBB['C1_Tavg'] = sjdfBB[sjdfBB.columns[4+nc+nt:4+nc+nt+nc]].mean(axis=1)
    sjdfBB['C2_Tavg'] = sjdfBB[sjdfBB.columns[4+nc+nt+nc:4+nc+nt+nc+nt]].mean(axis=1)
    ## get junction strength columns and diff ##
    sjdfBB['C1_JS'] = sjdfBB['C1_Javg']/sjdfBB['C1_Tavg']
    sjdfBB['C2_JS'] = sjdfBB['C2_Javg']/sjdfBB['C2_Tavg']
    sjdfBB['C1_JS'] = sjdfBB['C1_JS'].fillna(0)
    sjdfBB['C2_JS'] = sjdfBB['C2_JS'].fillna(0)
    sjdfBB['JS_diff']=sjdfBB['C2_JS']-sjdfBB['C1_JS']
    print(sjdfBB)
    if 'C1_JS' in sjdfBB.columns and 'C2_JS' in sjdfBB.columns and not sjdfBB[['C1_JS', 'C2_JS']].empty:
        sjdfBB.to_csv(sjdfBB_path, sep="\t", index=None)
        return(1)
    else:
        return(0)




