# Setup.py
import os, subprocess, sys

# Setup environment for CrypSplice run #

# Check Dependencies
## checking for dependencies and installing if necessary
def check_dependencies():
    dependency = ['pandas', 'numpy', 'pybedtools']
    
    for package in dependency:
        try:
            # Try to import the package to check if it's installed
            __import__(package)
            print(f"Package {package} found.")
        except ImportError:
            print(f"Package {package} not found. Installing ....\n")
            try:
                subprocess.run([sys.executable, '-m', 'pip', 'install', package], 
                             stderr=subprocess.DEVNULL, shell=False, check=True)
                print(f"Successfully installed {package}")
            except subprocess.CalledProcessError:
                print(f"Failed installing {package} ....\n")
                return package
    
    return None

# Check Directories
## checking for output directories and creating/cleaning if necessary
def check_directories(outDir):
    if os.path.isdir(outDir):
        return_code = subprocess.run(['rm', '-R', outDir], stdout=subprocess.PIPE, shell=False)
        return_code = subprocess.run(['mkdir', outDir], stdout=subprocess.PIPE, shell=False)
    else:
        return_code = subprocess.run(['mkdir', outDir], stdout=subprocess.PIPE, shell=False)
    return return_code

### Check Files
## checking that all files provided in the list are accessible
def check_files(file_list):
    for file in file_list:
        if os.path.exists(file):
            pass
        else:
            return file
    return None

# Initialize Samples
## splitting the user-provided comma delimited paths of samples into vectors
def initialize_samples(controls, treated):
    control_list = "".join(controls).replace(" ", "").split(",")
    treated_list = "".join(treated).replace(" ", "").split(",")
    control_num = len(control_list)
    treated_num = len(treated_list)
    return control_list, treated_list, control_num, treated_num

# CL Clust Initialize Samples
## CrypticLoad Clust specific - splitting the user-provided comma delimited paths of samples into vectors
def CLclust_initialize_samples(samples):
    LoadSamples = "".join(samples).replace(" ", "").split(",")
    if (len(LoadSamples) % 2) == 0:
        control_list = LoadSamples[0:int(len(LoadSamples)/2)]
        treated_list = LoadSamples[int(len(LoadSamples)/2):]
    else:
        control_list = LoadSamples[0:int((len(LoadSamples)/2)+1)]
        treated_list = LoadSamples[int((len(LoadSamples)/2)+1):]
    
    control_num = len(control_list)
    treated_num = len(treated_list)
    return control_list, treated_list, control_num, treated_num