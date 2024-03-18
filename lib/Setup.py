# Setup.py

import pkg_resources, os, subprocess


# Setup environment for CrypSplice run #


# Check Dependencies
## checking for dependencies and installing if necessary
def check_dependencies():
    dependency=['pandas','numpy','pybedtools']
    installed_packages = " ".join([str(d) for d in pkg_resources.working_set])
    for package in dependency:
        if package in installed_packages:
            pass
        else:
            print("Package "+package+" not found. Installing ....\n")
            try:        
                subprocess.run(['pip3','install',package],stderr=subprocess.DEVNULL,shell=False)
            except:
                print("Failed installing "+package+" ....\n")
                return(package)      


# Check Directories
## checking for output directories and creating/cleaning if necessary
def check_directories(outDir):
    if (os.path.isdir(outDir)):
        return_code=subprocess.run(['rm','-R',outDir],stdout=subprocess.PIPE,shell=False)
        return_code=subprocess.run(['mkdir',outDir],stdout=subprocess.PIPE,shell=False)
    else:
        return_code=subprocess.run(['mkdir',outDir],stdout=subprocess.PIPE,shell=False)
    return(return_code)


### Check Files
## checking that all files provided in the list are accessible 
def check_files (file_list):
    for file in file_list:
        if os.path.exists(file):
            pass
        else:
            return(file)
       

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
    LoadSamples="".join(samples).replace(" ","").split(",")
    if (len(LoadSamples) % 2) == 0:
        control_list=LoadSamples[0:int(len(LoadSamples)/2)]
        treated_list=LoadSamples[int(len(LoadSamples)/2):]
    else:
        control_list=LoadSamples[0:int((len(LoadSamples)/2)+1)]
        treated_list=LoadSamples[int((len(LoadSamples)/2)+1):]
    control_num = len(control_list)
    treated_num = len(treated_list)
    return control_list, treated_list, control_num, treated_num







   
