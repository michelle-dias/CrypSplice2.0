#!/usr/bin/env python3
"""
CrypSplice 2.1 - Main functionality
"""

# Standard library imports
import os, sys, glob, time
import pandas as pd, numpy as np, argparse
import concurrent.futures as cf, subprocess
import warnings, logging
from collections import defaultdict
import pysam

# Import modules from lib package
try:
    # When installed as a package, use relative imports
    from .lib import Setup, LogFile
    from .lib import ExtractJunctions
    from .lib import AddGenes
    from .lib import AnnotateJunctions
    from .lib import DifferentialUsage
    from .lib import FilterJunctions
    from .lib import CrypticLoad
    from .lib import PlotJunctions
    from .lib import BatchCorrection
except ImportError:
    # When running directly (not as package), use the old method
    sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
    import Setup, LogFile
    import ExtractJunctions
    import AddGenes
    import AnnotateJunctions
    import DifferentialUsage
    import FilterJunctions
    import CrypticLoad
    import PlotJunctions
    import BatchCorrection

# Configure warnings
warnings.filterwarnings(action="ignore", category=UserWarning, module="pysam")
warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None



#################### (Main Function) ####################
def main():
    
    
    ######## CrypSplice2.1 ########


    #### Initialization 


    ### Argument Parser
    # creating argument parser
    usage = '\r{}\nUsage:           %(prog)s <command> [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))
    parser = argparse.ArgumentParser(prog='CrypSplice2.1',usage=usage, formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('''Command:         CrypticJunctions        Infer cryptic splice junctions (across two conditions)
                 CrypticLoad             Infer cryptic load at gene/sample level and stratify samples
                 PlotJunctions           Visualize junctions of interest'''), epilog= " ")
    # creating sub-parser for main commands [CrypticJunctions, PlotJunctions, CrypticLoad] 
    sub_parsers = parser.add_subparsers(title="Commands",dest="command", help=argparse.SUPPRESS)
    
    ## CrypticJunctions sub-parser
    CJusage = '\r{}\nUsage:           %(prog)s [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))
    CJ = sub_parsers.add_parser(prog='CrypSplice2.1 CrypticJunctions', help='Infer cryptic splice junctions',name="CrypticJunctions", usage=CJusage, epilog=" ",formatter_class=argparse.RawDescriptionHelpFormatter)
    optional = CJ._action_groups.pop()
    required = CJ.add_argument_group('required arguments')
    CJ._action_groups.append(optional)
    # required arguments
    required.add_argument('-c1',help='Comma-separated list of condition1 files [control files]', nargs='+', required='True', type=str)
    required.add_argument('-c2',help='Comma-separated list of condition2 files [treated files]', nargs='+', required='True', type=str)
    required.add_argument('-gtf',help='Reference gtf file',required='True',type=str)
    required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
    required.add_argument('-s',help='Strand specific data 0: unstranded, 1: frd, 2: rev', choices=[0,1,2], required='True', type=int)
    required.add_argument('-o',help='Out put directory',type=str)
    # optional arguments
    optional.add_argument('-prefix',help='Output file prefix', default="CrypSplice",type=str)
    optional.add_argument('-p',help='No. of processors to use',type=int,default=10)
    optional.add_argument('-annotated',help='Infer annotated junction changes', action='store_true')
    # optional arguments - filters
    optional.add_argument('-j',help='Read cutoff. Junctions with readcounts below this threshold are ignored',type=int,default=10)
    optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=50)
    optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.5)
    optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=15)
    optional.add_argument('-g', help='Gene filter. Junctions spanning multiple genes or not attributed to a gene are removed from analysis', action='store_true')
    optional.add_argument('-b', help="Path to text file containing batch information (1st column: 'Sample', 2nd column: 'Batch')")

    ## CrypticLoad sub-parser
    CLusage = '\r{}\nUsage:           %(prog)s <sub-command> [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))
    CL = sub_parsers.add_parser(prog='CrypSplice2.1 CrypticLoad',name="CrypticLoad", usage=CLusage, formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('''Sub-Command:     Diff         Differential CrypticLoad across two conditions
                 Clust         Stratify samples based on CrypticLoad'''), epilog=" ")
    ## CrypticLoad sub-commands [Diff, Clust] sub-parser
    CLsub = CL.add_subparsers(help=argparse.SUPPRESS,dest="CLcommand")
    ## CrypticLoad Diff sub-parser
    CLDiffusage = '\r{}\nUsage:           %(prog)s [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))   
    CLdiff = CLsub.add_parser(prog='CrypSplice2.1 CrypticLoad Diff', help='Differential CrypticLoad across two conditions', name="Diff", usage=CLDiffusage, epilog=" ",formatter_class=argparse.RawDescriptionHelpFormatter)
    CLdiff_optional = CLdiff._action_groups.pop()
    CLdiff_required = CLdiff.add_argument_group('Required arguments')
    CLdiff._action_groups.append(CLdiff_optional)
    # required arguments
    CLdiff_required.add_argument('-c1',help='Comma-separated list of condition1 files', nargs='+', required='True', type=str)
    CLdiff_required.add_argument('-c2',help='Comma-separated list of condition2 files', nargs='+', required='True', type=str)
    CLdiff_required.add_argument('-gtf',help='Reference gtf file',required='True',type=str)
    CLdiff_required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
    CLdiff_required.add_argument('-s',help='Strand specific data 0: unstranded, 1: frd, 2: rev', choices=[0,1,2], required='True', type=int)
    CLdiff_required.add_argument('-o',help='Out put directory',type=str)
    # optional arguments
    CLdiff_optional.add_argument('-prefix',help='Output file prefix', default="CrypSplice",type=str)
    CLdiff_optional.add_argument('-p',help='No. of processors to use',type=int,default=10)
    CLdiff_optional.add_argument('-b', help="Path to text file containing batch information (1st column: 'Sample', 2nd column: 'Batch')")
    # optional arguments - filters
    CLdiff_optional.add_argument('-j',help='Read cutoff. Junctions with readcounts below this threshold are ignored',type=int,default=10)
    CLdiff_optional.add_argument('-d',help='Read cutoff. Junctions with origincounts below this threshold are ignored',type=int,default=0)
    CLdiff_optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=50)
    CLdiff_optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.5)
    CLdiff_optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=0)

    #CLdiff_optional.add_argument('-g', help='Gene filter. Junctions spanning multiple genes or not attributed to a gene are removed from analysis', action='store_true')
    # ## CrypticLoad Clust sub-parser
    CLclstusage = '\r{}\nUsage:           %(prog)s [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))       
    CLclst = CLsub.add_parser(prog='CrypSplice2.1 CrypticLoad Clust', help='Stratify samples based on CrypticLoad', name="Clust", usage=CLclstusage, epilog=" ", formatter_class=argparse.RawDescriptionHelpFormatter)
    CLclst_optional = CLclst._action_groups.pop()
    CLclst_required = CLclst.add_argument_group('required arguments')
    CLclst._action_groups.append(CLclst_optional)
    # required arguments
    CLclst_required.add_argument('-samples',help='Comma-separated list of files. Index files are expected at the same location', nargs='+', required='True', type=str)
    CLclst_required.add_argument('-gtf',help='Reference gtf file',required='True',type=str)
    CLclst_required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
    CLclst_required.add_argument('-s',help='Strand specific data 0: unstranded, 1: frd, 2: rev', choices=[0,1,2], required='True', type=int)
    CLclst_required.add_argument('-o',help='Out put directory',type=str)
    # optional arguments
    CLclst_optional.add_argument('-r',help='User-provided rank (NMF cluster) instead of automatically inferring rank',type=int)
    CLclst_optional.add_argument('-i',help='Iterations of NMF to infer best-fit rank',type=int, default=50)
    CLclst_optional.add_argument('-prefix',help='Output file prefix', default="CrypSplice",type=str)
    CLclst_optional.add_argument('-p',help='No. of processors to use',type=int,default=10)
    CLclst_optional.add_argument('-b', help="Path to text file containing batch information (1st column: 'Sample', 2nd column: 'Batch')")

    # optional arguments - filters
    CLclst_optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=50)
    CLclst_optional.add_argument('-j',help='Read cutoff. Junctions with readcounts bellow this threshold are ignored',type=int,default=10)
    CLclst_optional.add_argument('-d',help='Read cutoff. Junctions with origincounts bellow this threshold are ignored',type=int,default=0)
    CLclst_optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.5)
    CLclst_optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=0)

    ## PlotJunctions sub-parser
    PJusage = '\r{}\nUsage:           %(prog)s [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))
    PJ = sub_parsers.add_parser(prog='CrypSplice2.1 PlotJunctions', help='Visualize cryptic junctions of interest',name="PlotJunctions", usage=PJusage, epilog=" ",formatter_class=argparse.RawDescriptionHelpFormatter)
    PJ_optional = PJ._action_groups.pop()
    PJ_required = PJ.add_argument_group('required arguments')
    PJ._action_groups.append(PJ_optional)
     # required arguments
    PJ_required.add_argument('-o',help='Output directory where plots will be', required='True', type=str)
    PJ_required.add_argument('-pj',help='Path to output file from Crypticjunctions command containing junctions to be plotted [Novel_Junctions.txt or Annotated_Junctions.txt]', required='True', type=str)
    PJ_required.add_argument('-c1',help='Comma-separated list of condition1 files [control files]', nargs='+', required='True', type=str)
    PJ_required.add_argument('-c2',help='Comma-separated list of condition2 files [treated files]', nargs='+', required='True', type=str)
    PJ_required.add_argument('-gtf',help='Reference gtf file',required='True',type=str)
    PJ_required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
    PJ_required.add_argument('-s',help='Strand specific data 0: unstranded, 1: frd, 2: rev', choices=[0,1,2], required='True', type=int)
    # optional arguments
    PJ_optional.add_argument('-b', help="Path to text file containing batch information (1st column: 'Sample', 2nd column: 'Batch')")
    PJ_optional.add_argument('-prefix',help='Output file prefix', default="CrypSplice",type=str)
    PJ_optional.add_argument('-p',help='No. of processors to use',type=int,default=10)
    PJ_optional.add_argument('-top_n', help='Number of top significant junction ids from Novel_Junctions.txt to plot',type=int, default=10)
    PJ_optional.add_argument('-j_ids', help='Comma-separated list of specific junction ids from Novel_Junctions.txt input to plot', type=str, default=None)
    PJ_optional.add_argument('--merge', help='Enable merging of bigWigs for C1 and C2', action='store_true', default=False)
    PJ_optional.add_argument('-c1_names', help='Comma-separated list of c1 names in the same order as the c1 paths, if merge enabled pass single string to represent c1',  default=None)
    PJ_optional.add_argument('-c2_names', help='Comma-separated list of c2 names in the same order as the c2 paths, if merge enabled pass single string to represent c2', default=None)
    PJ_optional.add_argument('-c1_color', help='Pass hex code for color of C1 bigWigs', type=str, default='#2387de')
    PJ_optional.add_argument('-c2_color', help='Pass hex code for color of C1 bigWigs', type=str, default='#e04196')
    # optional arguments - filters
    PJ_optional.add_argument('-j',help='Read cutoff. Junctions with readcounts below this threshold will not be plotted',type=int,default=10)
    PJ_optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will not be plotted',type=int,default=50)
    

    # getting arguments passed by user
    args=parser.parse_args()


    ### Help Message 
    if len(sys.argv)==1:
        print(" ")
        parser.print_help(sys.stderr) 
        sys.exit(1)
    print(" ")
    if args.command == "CrypticLoad" and args.CLcommand is None:
        print(" ")
        CL.print_help(sys.stderr)
        sys.exit(1)

    # Adding delimiter to prefix if passed
    args.prefix += "." if args.prefix else ""

    ### Setup [Module] : Initial setup of environment for run
    
    # check_dependencies: checking for dependencies and installing if necessary
    missing_dependency = Setup.check_dependencies()
    if missing_dependency is not None:
        print("Unable to install "+missing_dependency)
        exit()

    # check_directories: checking for output directories and creating if necessary (will "clean" or remove directory if already existing and replace)
    args.o += "/" if not args.o.endswith("/") else ""
    if Setup.check_directories(args.o) == 1:
        print("Failed to find "+args.o)
        exit()

    # initialize_samples: splitting the comma delimited paths of samples into sample vectors
    if args.command == 'CrypticLoad' and args.CLcommand == 'Clust':
        controls, treated, control_num, treated_num = Setup.CLclust_initialize_samples(args.samples)        
    else:
        controls, treated, control_num, treated_num = Setup.initialize_samples(args.c1, args.c2)
    
    # check_files: checking that required files are accessible
    if args.command == 'PlotJunctions':
        checkfiles=controls+treated+[args.fasta]+[args.gtf]+glob.glob(args.pj)
    ## check controls, treated, fasta, and gtf for all other commands ##
    else:
        checkfiles = controls+treated+[args.fasta]+[args.gtf]
    missing_file = Setup.check_files(checkfiles)
    if missing_file is not None:
        print("Can not access "+missing_file)
        exit()
   



    ### LogFile [Module] : Initializing log file
    # write_arguments: initializing logfile with transcription of arguments passed by user
    LogFile.write_arguments(args.o+args.prefix, args)
    logfile_path = str(args.o+args.prefix+'log.txt')


    # execution of sub-commands [CrypticJunctions, CrypticLoad, etc. ]
    try:

        #### CrypticJunctions or CrypticLoad
        if args.command in ['CrypticJunctions', 'CrypticLoad', 'PlotJunctions']:


            ### ExtractJunctions [Module] : Extracting junction counts using pysam 
            if ExtractJunctions.get_junction_counts(args.o + args.prefix, args.p, controls, treated, args.s, args.gtf, args.l) == 1:
                LogFile.log_message(logfile_path, "Completed junction extraction")
                #junction_counts = pd.read_csv(args.o + args.prefix + "JunctionCounts.txt", sep="\t", dtype={"chrom": str})
            else:
                LogFile.log_message(logfile_path, "Failed junction extraction")
                LogFile.log_message(logfile_path, "Terminating...")
                exit()
            

            ### FilterJunctions [Module] : Extra junction filter to improve capture of significant junctions             
            # setting junctions below the junction cutoff to 0 in the regtools extraction files
            if args.command in ['CrypticLoad']:
                if FilterJunctions.juncCut_filter(args.o+args.prefix, args.p, args.j, args.j) == 1:
                    LogFile.log_message(logfile_path, "Completed filtering for junctions where junction counts and origin counts are less than the junction cutoff : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed filtering for junctions where junction counts and origin counts are less than the junction cutoff : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
            # # this is necessary for big dataframes (ROSMAP) where pandas can't handle all the data for too many samples 
            # # setting junctions below 5 to 0 in the regtools extraction files
            # if (treated_num + control_num) >= 30:
            #     if FilterJunctions.bigData_filter(args.o+args.prefix, args.p, 5, 5) == 1:
            #         LogFile.log_message(logfile_path, "Completed filtering for junctions where junction counts and origin counts are less than 5 [Big Data Filter] : ")
            #         pass
            #     else:
            #         LogFile.log_message(logfile_path, "Failed filtering for junctions where junction counts and origin counts are less than 5 [Big Data Filter] : ")
            #         LogFile.log_message(logfile_path, "Terminating ............... : ")
            #         exit()
            
            ### ExtractJunctions [Module] : Extracting junction counts using pysam 
            # stitch extract junction count files together
            if ExtractJunctions.stitch_extractJunctions(controls, treated, args.o + args.prefix) == 1:
                junction_counts = pd.read_csv(args.o + args.prefix + "JunctionCounts.txt", sep="\t")
                LogFile.log_message(logfile_path, "Completed stitching extracted junction counts together : ")
                pass
            else:
                LogFile.log_message(logfile_path, "Failed stitching extracted junction counts together : ")
                LogFile.log_message(logfile_path, "Terminating ............... : ")
                exit()


            ### BatchCorrection [Module] :
            if args.b:
                if BatchCorrection.run_batchCorr(args.o+args.prefix+"JunctionCounts.txt", args.b)==1:
                    junction_counts = pd.read_csv(args.o + args.prefix + "JunctionCounts.txt", sep="\t")
                    LogFile.log_message(logfile_path, "Completed batch correcting junction counts and origin counts : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed batch correcting junction counts and origin counts : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()





            ### FilterJunctions [Module] : Extra junction filter to improve capture of significant junctions 
            # filtering at this step for CrypticJunctions
            if args.command in ['CrypticJunctions']:
                # PoverAM filter
                junction_counts = FilterJunctions.PoverAM_filter(junction_counts, control_num, treated_num, args.pa_p, args.pa_a, args.j)
                if isinstance(junction_counts, pd.DataFrame):
                    LogFile.log_message(logfile_path, "Completed PoverA filtering : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed PoverA filtering : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                junction_counts.to_csv(args.o+args.prefix+"JunctionCounts.txt", sep="\t", index=None)


            ### AnnotateJunctions [Module] : Assign genes (ids/names) to each junction
            if AnnotateJunctions.annotate_junctions(args.o+args.prefix+"JunctionCounts.txt", args.gtf, args.fasta, args.o+args.prefix) == 1:                
                LogFile.log_message(logfile_path, "Completed annotating junctions: ")
                junction_counts = pd.read_csv(args.o+args.prefix+"Annotated_JunctionCounts.txt", sep="\t", dtype={"chrom": str})
                pass
            else:
                LogFile.log_message(logfile_path, "Failed to annotate junctions: ")
                LogFile.log_message(logfile_path, "Terminating ............... : ")
                exit()
            junction_counts = pd.read_csv(args.o+args.prefix+"Annotated_JunctionCounts.txt", sep="\t", dtype={"chrom": str})


            ### AddGenes [Module] : Assign genes (ids/names) to each junction
            if AddGenes.add_genes(junction_counts, args.gtf, args.p, args.o+args.prefix) == 1:
                LogFile.log_message(logfile_path, "Completed adding junction gene ids : ")
                junction_counts = pd.read_csv(args.o + args.prefix + "Annotated_JunctionCounts.txt", sep="\t", dtype={"chrom": str})
                pass
            else:
                LogFile.log_message(logfile_path, "Failed adding junction gene annotations : ")
                LogFile.log_message(logfile_path, "Terminating ............... : ")
                exit()


            # ### FilterJunctions [Module] : Extra junction filter to improve capture of significant junctions 
            # # gene filter (if enabled by user) or Cryptic Load
            # if (args.command in ['CrypticLoad'] or args.g):
            #     junction_counts = FilterJunctions.gene_filter(junction_counts)
            #     if isinstance(junction_counts, pd.DataFrame):
            #         LogFile.log_message(logfile_path, "Completed filtering junctions spanning multiple genes or not attributed to any gene: ")
            #         junction_counts.to_csv(args.o + args.prefix + "Annotated_JunctionCounts.txt", sep="\t", index=None)
            #         pass
            #     else:
            #         LogFile.log_message(logfile_path, "Failed filtering junctions spanning multiple genes or not attributed to any gene: ")
            #         LogFile.log_message(logfile_path, "Terminating ............... : ")
            #         exit()


          
            #### CrypticJunctions 
            if args.command in ['CrypticJunctions']:

                ### DifferentialUsage [Module] : Run beta-binomial test (R countdata) to identify differentially used junctions
                ### DifferentialUsage [Module] : Calculate junction strength (PSI) (junction counts)/(total counts originating from junction origin)
                # novel junctions
                novel_junctions = junction_counts[junction_counts['annotation']!="DA"]
                novel_junctions.to_csv(args.o+args.prefix+"Novel_Junctions.txt", sep="\t", index=None)
                if DifferentialUsage.run_bbTest(args.o+args.prefix+"Novel_Junctions.txt", control_num, treated_num, args.p, "CrypticJunctions") == 1:
                    LogFile.log_message(logfile_path, "Completed BB test on novel junctions : ")
                    if DifferentialUsage.calc_junction_strength(args.o+args.prefix+"Novel_Junctions.txt", control_num, treated_num, 4) == 1:
                        LogFile.log_message(logfile_path, "Completed computing junction strength on novel junctions : ")
                        novel_junctions = pd.read_csv(args.o+args.prefix+"Novel_Junctions.txt", sep="\t")
                        pass
                    else:
                        LogFile.log_message(logfile_path, "Failed computing junction strength on novel junctions : ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()
                else:
                    LogFile.log_message(logfile_path, "Failed BB test on novel junctions : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()

                # annotated junctions (if enabled by user)
                if args.annotated:
                    annotated_junctions = junction_counts[junction_counts['annotation']=="DA"]
                    annotated_junctions.to_csv(args.o+args.prefix+"Annotated_Junctions.txt", sep="\t", index=None)
                    if DifferentialUsage.run_bbTest(args.o+args.prefix+"Annotated_Junctions.txt", control_num, treated_num, 5, "CrypticJunctions") == 1:
                        LogFile.log_message(logfile_path, "Completed BB test on annotated junctions : ")
                        if DifferentialUsage.calc_junction_strength(args.o+args.prefix+"Annotated_Junctions.txt", control_num, treated_num, 4) == 1:
                            LogFile.log_message(logfile_path, "Completed computing junction strength on annotated junctions : ")
                            novel_junctions = pd.read_csv(args.o+args.prefix+"Annotated_Junctions.txt", sep="\t")
                            pass
                        else:
                            LogFile.log_message(logfile_path, "Failed computing junction strength on annotated junctions : ")
                            LogFile.log_message(logfile_path, "Terminating ............... : ")
                            exit()
                    else:
                        LogFile.log_message(logfile_path, "Failed BB test on annotated junctions : ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()

                # clean up intermediate files 
                intermediate_files = [args.o+args.prefix+"Annotated_JunctionCounts.txt", args.o+args.prefix+"JunctionCounts.txt"]
                [os.unlink(file) for file in intermediate_files]



            #### CrypticLoad 
            if args.command in ['CrypticLoad']:
                ### CrypticLoad [Module] : Calculate gene and sample-level load (equivalent of PSI)
                # gene-level
                if CrypticLoad.get_geneLoad(args.o+args.prefix+"Annotated_JunctionCounts.txt", control_num, treated_num, args.o+args.prefix) == 1:
                    LogFile.log_message(logfile_path, "Completed calculating gene-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed calculating gene-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # sample-level
                if CrypticLoad.get_sampleLoad(args.o+args.prefix+"Annotated_JunctionCounts.txt", control_num, treated_num, args.o+args.prefix) == 1:
                    LogFile.log_message(logfile_path, "Completed calculating sample-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed calculating sample-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()

                ### FilterJunctions [Module] : Extra junction filter to retain only genes that could have cryptic significance to improve significances
                # greater than 0 (a) in 0.5 (p) percent of the reads
                gene_loads = pd.read_csv(args.o+args.prefix+"GeneLoad.txt", sep="\t")
                gene_loads = FilterJunctions.PoverA_filter_CL(gene_loads, control_num, treated_num, args.pa_p, args.pa_a)
                if isinstance(gene_loads, pd.DataFrame):
                    LogFile.log_message(logfile_path, "Completed PoverA filtering for gene loads : ")
                    gene_loads.to_csv(args.o+args.prefix+"GeneLoad.txt", sep="\t", index=None)
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed PoverA filtering for gene loads : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()


            ### DifferentialUsage [Module] : Run beta-binomial test (R countdata) to identify genes/samples with differential loads between treatment groups
            if args.command in ['CrypticLoad'] and args.CLcommand =="Diff":
                # gene-level
                if DifferentialUsage.run_bbTest(args.o+args.prefix+"GeneLoad.txt", control_num, treated_num, args.p, "CrypticLoad") == 1:
                    LogFile.log_message(logfile_path, "Completed BB test on gene-level cryptic loads : ")
                    if DifferentialUsage.calc_junction_strength(args.o+args.prefix+"GeneLoad.txt", control_num, treated_num, 1) == 1:
                        LogFile.log_message(logfile_path, "Completed computing cryptic strength on gene-level cryptic loads : ")
                        pass
                    else:
                        LogFile.log_message(logfile_path, "Failed computing cryptic strength on gene-level cryptic loads : ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()
                else:
                    LogFile.log_message(logfile_path, "Failed BB test on gene-level cryptic loads : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # sample-level
                if DifferentialUsage.run_bbTest(args.o+args.prefix+"SampleLoad.txt", control_num, treated_num, args.p, "CrypticLoad") == 1:
                    LogFile.log_message(logfile_path, "Completed BB test on sample-level cryptic loads : ")
                    if DifferentialUsage.calc_junction_strength(args.o+args.prefix+"SampleLoad.txt", control_num, treated_num, 1) == 1:
                        LogFile.log_message(logfile_path, "Completed computing cryptic strength on sample-level cryptic loads : ")
                        pass
                    else:
                        LogFile.log_message(logfile_path, "Failed computing cryptic strength on sample-level cryptic loads : ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed BB test on sample-level cryptic loads : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # clean up intermediate files 
                intermediate_files = [args.o+args.prefix+"Annotated_JunctionCounts.txt", args.o+args.prefix+"JunctionCounts.txt"]
                [os.unlink(file) for file in intermediate_files]


            ### CrypticLoad [Module] : Clustering gene loads 
            if args.command in ['CrypticLoad'] and args.CLcommand =="Clust":
                # gene-level
                if CrypticLoad.cluster_geneLoad(args.o+args.prefix+"GeneLoad.txt", args.i, args.p, args.o+args.prefix+"GeneLevel.", args.r) == 1:
                    LogFile.log_message(logfile_path, "Completed clustering gene-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed clustering gene-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()

                # sample-level
                if CrypticLoad.cluster_sampleLoad(args.o+args.prefix+"SampleLoad.txt", args.o+args.prefix+"SampleLevel.", args.r) == 1:
                    LogFile.log_message(logfile_path, "Completed clustering sample-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed clustering sample-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # clean up intermediate files 
                intermediate_files = [args.o+args.prefix+"Annotated_JunctionCounts.txt", args.o+args.prefix+"JunctionCounts.txt"]
                [os.unlink(file) for file in intermediate_files]




            ### PlotJunctions [Module] : Plotting junctions of interest 
            if args.command in ['PlotJunctions']:
                if PlotJunctions.bam_to_bigWig(controls, treated, args.fasta, args.o+args.prefix, args.p, args.merge)==1:
                    LogFile.log_message(logfile_path, "Completed converting bams to bigWigs : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed converting bams to bigWigs : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                if PlotJunctions.create_sashimis(control_num, treated_num, args.o+args.prefix, args.merge)==1:
                    LogFile.log_message(logfile_path, "Completed creating sashimi files : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed to create sashimi files : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                return_val = PlotJunctions.create_plots(controls, treated, args.gtf, args.o+args.prefix, args.merge, args.pj, args.j_ids, args.top_n, args.c1_color, args.c2_color, args.c1_names, args.c2_names)
                if return_val == 1:    
                    LogFile.log_message(logfile_path, "Completed plotting junctions : ")
                    pass
                elif return_val == 2:
                    LogFile.log_message(logfile_path, "No significant junctions founds, please provide specific IDs : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                else: 
                    LogFile.log_message(logfile_path, "Failed to plot junctions : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()

                # clean up intermediate files 
                intermediate_files = [args.o+args.prefix+"Annotated_JunctionCounts.txt", args.o+args.prefix+"JunctionCounts.txt"]
                [os.unlink(file) for file in intermediate_files]











       









        









#         # Tidyup 
#         # del_files=glob.glob(args.o+"*.bed")+glob.glob(args.o+"*.arcs")+glob.glob(args.o+"*.out.tab")+glob.glob(args.o+args.prefix+"_Junctions_*.txt")
#         # log = subprocess.run(['rm']+del_files,stderr=subprocess.DEVNULL,shell=False)    
#         # localdate = time.strftime('%a %m/%d/%Y')
#         # localtime = time.strftime('%H:%M:%S')
#         # logfile.write('# Finished CrypSplice run: '+localdate+' at: ' + localtime+' \n')
#         # logfile.close()
    
    finally:
        print("Fin")
        # #pass
        # # Clean if terminated #
        # del_files=glob.glob(args.o+"*.bed")+glob.glob(args.o+"*.arcs")+glob.glob(args.o+"*.out.tab")+glob.glob(args.o+args.prefix+"_Junctions_*.txt")
        # for file in del_files:
        #   if os.path.exists(file):
        #       subprocess.run(['rm']+[file],stderr=subprocess.DEVNULL,shell=False)
        #       pass

if __name__ == "__main__":
    try:
        main()
    except ImportError as error:
        print(f"\nPython 3 or a more recent version is required: {error}\nTerminating ...\n")