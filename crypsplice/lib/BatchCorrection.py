# Batch Correcting using the ComBat-Seq R package #




import sys, os
import subprocess
import pandas as pd




# Run Batch Corrections
## calling the countdata R function to run combat-seq to batch correct
def run_batchCorr(junctions_path, bc_meta):
	Rscript_path = "/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib/RunBatchCorrection.R"
	arguments = [junctions_path, bc_meta]
	command = ["Rscript", Rscript_path] + arguments
	result = subprocess.check_output(command, universal_newlines=True)
	return int(result.strip().split("\n")[-1])



