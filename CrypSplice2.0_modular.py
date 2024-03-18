# CrypSplice.py



#################### (Main Function) ####################
def main():
    
    
    ######## CrypSplice2.1 ########


    #### Initialization 


    ### Argument Parser
    # creating argument parser
    usage = '\r{}\nUsage:           %(prog)s <command> [arguments]'.format('Program:         CrypSplice\nVersion:         2.1'.ljust(len('usage:')))
    parser = argparse.ArgumentParser(prog='CrypSplice2.1',usage=usage, formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('''Command:         CrypticJunctions        Infer cryptic splice junctions (across two conditions)
                 CrypticLoad             Infer cryptic load at gene/sample level and stratify samples'''), epilog= " ")
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
    optional.add_argument('-annotated',help='Infer annotated junction changes',type=int,choices=[0,1],default=0)
    # optional arguments - filters
    optional.add_argument('-j',help='Read cutoff. Junctions with readcounts bellow this threshold are ignored',type=int,default=10)
    optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=25)
    optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.5)
    optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=15)
    optional.add_argument('-g', help='Gene filter. Junctions spanning multiple genes or not attributed to a gene are removed from analysis', action='store_true')
    
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
    # optional arguments - filters
    # #CLdiff_optional.add_argument('-j',help='Read cutoff. Junctions with readcounts below this threshold are ignored',type=int,default=10)
    # CLdiff_optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=25)
    # CLdiff_optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.5)
    # CLdiff_optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=15)
    # CLdiff_optional.add_argument('-g', help='Gene filter. Junctions spanning multiple genes or not attributed to a gene are removed from analysis', action='store_true')
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
    CLclst_optional.add_argument('-prefix',help='Output file prefix', default="CrypSplice_Load",type=str)
    CLclst_optional.add_argument('-p',help='No. of processors to use',type=int,default=10)
    # optional arguments - filters
    CLclst_optional.add_argument('-j',help='Read cutoff. Junctions with readcounts bellow this threshold are ignored',type=int,default=10)
    CLclst_optional.add_argument('-l',help='Junction length cutoff. Junctions with intron length < cutoff will be dropped',type=int,default=25)
    
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



    ### Setup [Module] : Initial setup of environment for run
    
    # check_dependencies: checking for dependencies and installing if necessary
    missing_dependency = Setup.check_dependencies()
    if missing_dependency is not None:
        print("Unable to install "+missing_dependency)
        exit()

    # check_directories: checking for output directories and creating/cleaning if necessary
    # args.o += "/" if not args.o.endswith("/") else ""
    # if Setup.check_directories(args.o) == 1:
    #     print("Failed to find "+args.o)
    #     exit()

    # initialize_samples: splitting the comma delimited paths of samples into sample vectors
    if args.command == 'CrypticLoad' and args.CLcommand == 'Clust':
        controls, treated, control_num, treated_num = Setup.CLclust_initialize_samples(args.samples)        
    else:
        controls, treated, control_num, treated_num = Setup.initialize_samples(args.c1, args.c2)
    
    # check_files: checking that required files are accessible
    if args.command == 'PlotJunctions':
        checkfiles=controls+treated+[args.gtf]+glob.glob(args.cj+"*_Novel.txt")+glob.glob(args.cj+"*_Annotated.txt")
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
    logfile_path = str(args.o+args.prefix+'.CrypSplice.log.txt')


    # execution of sub-commands [CrypticJunctions, CrypticLoad, etc. ]
    try:

        #### CrypticJunctions or CrypticLoad
        if args.command in ['CrypticJunctions', 'CrypticLoad']:

            ### ExtractJunctions [Module] : Extracting junction counts using pysam 
            # turning off filtering at this step for CrypticLoad Clust
            if args.command in ['CrypticLoad']:
                junction_extraction_args = (args.o + args.prefix, args.p, controls, treated, args.s, args.gtf, 0, args.l)
            else:
                junction_extraction_args = (args.o + args.prefix, args.p, controls, treated, args.s, args.gtf, args.j, args.l)
            if ExtractJunctions.get_junction_counts(*junction_extraction_args) == 1:
                LogFile.log_message(logfile_path, "Completed junction extraction")
                junction_counts = pd.read_csv(args.o + args.prefix + ".JunctionCounts.txt", sep="\t", dtype={"chrom": str})
            else:
                LogFile.log_message(logfile_path, "Failed junction extraction")
                LogFile.log_message(logfile_path, "Terminating...")
                exit()

        
            ### AddGenes [Module] : Assign genes (ids/names) to each junction
            if AddGenes.add_genes(junction_counts, args.gtf, args.p, args.o+args.prefix) == 1:
                LogFile.log_message(logfile_path, "Completed adding junction gene annotations : ")
                junction_counts = pd.read_csv(args.o + args.prefix + ".JunctionCounts.txt", sep="\t", dtype={"chrom": str})
                pass
            else:
                LogFile.log_message(logfile_path, "Failed adding junction gene annotations : ")
                LogFile.log_message(logfile_path, "Terminating ............... : ")
                exit()


            ### FilterJunctions [Module] : Extra junction filter to improve capture of significant junctions 
            # special filtering at this step for CrypticLoad Clust
            if args.command in ['CrypticLoad']:
                junction_counts = FilterJunctions.gene_filter(junction_counts)
                if isinstance(junction_counts, pd.DataFrame):
                    LogFile.log_message(logfile_path, "Completed filtering junctions spanning multiple genes or not attributed to any gene: ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed filtering junctions spanning multiple genes or not attributed to any gene: ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                junction_counts = FilterJunctions.load_quality_filter(junction_counts, args.j)
                if isinstance(junction_counts, pd.DataFrame):
                    LogFile.log_message(logfile_path, "Completed filtering for quality junctions : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed filtering for quality junctions : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
            else:
                # PoverA filter
                junction_counts = FilterJunctions.PoverA_filter(junction_counts, control_num, treated_num, args.pa_p, args.pa_a)
                if isinstance(junction_counts, pd.DataFrame):
                    LogFile.log_message(logfile_path, "Completed PoverA filtering : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed PoverA filtering : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # gene filter (if enabled by user)
                if args.g:
                    junction_counts = FilterJunctions.gene_filter(junction_counts)
                    if isinstance(junction_counts, pd.DataFrame):
                        LogFile.log_message(logfile_path, "Completed filtering junctions spanning multiple genes or not attributed to any gene: ")
                        pass
                    else:
                        LogFile.log_message(logfile_path, "Failed filtering junctions spanning multiple genes or not attributed to any gene: ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()


            ### AnnotateJunctions [Module] : Assign genes (ids/names) to each junction
            if AnnotateJunctions.annotate_junctions(junction_counts, args.gtf, args.o+args.prefix, args.s, args.p) == 1:                
                LogFile.log_message(logfile_path, "Completed annotating junctions: ")
                junction_counts = pd.read_csv(args.o+args.prefix+".Annotated_JunctionCounts.txt", sep="\t", dtype={"chrom": str})
                pass
            else:
                LogFile.log_message(logfile_path, "Failed to annotate junctions: ")
                LogFile.log_message(logfile_path, "Terminating ............... : ")
                exit()

            

            #### CrypticJunctions 
            if args.command in ['CrypticJunctions']:

                ### DifferentialUsage [Module] : Run beta-binomial test (R countdata) to identify differentially used junctions
                ### DifferentialUsage [Module] : Calculate junction strength (PSI) (junction counts)/(total counts originating from junction origin)
                # novel junctions
                novel_junctions = junction_counts[junction_counts['annotation']!="DA"]
                novel_junctions.to_csv(args.o+args.prefix+".Novel_Junctions.txt", sep="\t", index=None)
                if DifferentialUsage.run_bbTest(args.o+args.prefix+".Novel_Junctions.txt", control_num, treated_num, args.p, "CrypticJunctions") == 1:
                    LogFile.log_message(logfile_path, "Completed BB test on novel junctions : ")
                    if DifferentialUsage.calc_junction_strength(args.o+args.prefix+".Novel_Junctions.txt", control_num, treated_num) == 1:
                        LogFile.log_message(logfile_path, "Completed computing junction strength on novel junctions : ")
                        novel_junctions = pd.read_csv(args.o+args.prefix+".Novel_Junctions.txt", sep="\t")
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
                if args.annotated == 1:
                    annotated_junctions = junction_counts[junction_counts['annotation']=="DA"]
                    annotated_junctions.to_csv(args.o+args.prefix+".Annotated_Junctions.txt", sep="\t", index=None)
                    if DifferentialUsage.run_bbTest(args.o+args.prefix+".Annotated_Junctions.txt", control_num, treated_num, args.p, "CrypticJunctions") == 1:
                        LogFile.log_message(logfile_path, "Completed BB test on annotated junctions : ")
                        if DifferentialUsage.calc_junction_strength(args.o+args.prefix+".Annotated_Junctions.txt", control_num, treated_num) == 1:
                            LogFile.log_message(logfile_path, "Completed computing junction strength on annotated junctions : ")
                            novel_junctions = pd.read_csv(args.o+args.prefix+".Annotated_Junctions.txt", sep="\t")
                            pass
                        else:
                            LogFile.log_message(logfile_path, "Failed computing junction strength on annotated junctions : ")
                            LogFile.log_message(logfile_path, "Terminating ............... : ")
                            exit()
                    else:
                        LogFile.log_message(logfile_path, "Failed BB test on annotated junctions : ")
                        LogFile.log_message(logfile_path, "Terminating ............... : ")
                        exit()


            #### CrypticLoad 
            if args.command in ['CrypticLoad']:
                ### CrypticLoad [Module] : Calculate gene and sample-level load (equivalent of PSI)
                # gene-level
                if CrypticLoad.get_geneLoad(args.o+args.prefix+".Annotated_JunctionCounts.txt", control_num, treated_num, args.o+args.prefix+".") == 1:
                    LogFile.log_message(logfile_path, "Completed calculating gene-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed calculating gene-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                # sample-level
                if CrypticLoad.get_sampleLoad(args.o+args.prefix+".Annotated_JunctionCounts.txt", control_num, treated_num, args.o+args.prefix+".") == 1:
                    LogFile.log_message(logfile_path, "Completed calculating sample-level cryptic load : ")
                    pass
                else:
                    LogFile.log_message(logfile_path, "Failed calculating sample-level cryptic load : ")
                    LogFile.log_message(logfile_path, "Terminating ............... : ")
                    exit()
                
            # ### DifferentialUsage [Module] : Run beta-binomial test (R countdata) to identify genes/samples with differential loads between treatment groups
            # if args.command in ['CrypticLoad'] and args.CLcommand =="Diff":
            #     # gene-level
            #     if DifferentialUsage.run_bbTest(args.o+args.prefix+".GeneLoad.txt", control_num, treated_num, args.p, "CrypticLoad") == 1:
            #         LogFile.log_message(logfile_path, "Completed BB test on gene-level cryptic loads : ")
            #         pass
            #     else:
            #         LogFile.log_message(logfile_path, "Failed BB test on gene-level cryptic loads : ")
            #         LogFile.log_message(logfile_path, "Terminating ............... : ")
            #         exit()
            #     # sample-level
            #     if DifferentialUsage.run_bbTest(args.o+args.prefix+".SampleLoad.txt", control_num, treated_num, args.p, "CrypticLoad") == 1:
            #         LogFile.log_message(logfile_path, "Completed BB test on sample-level cryptic loads : ")
            #         pass
            #     else:
            #         LogFile.log_message(logfile_path, "Failed BB test on sample-level cryptic loads : ")
            #         LogFile.log_message(logfile_path, "Terminating ............... : ")
            #         exit()

            # ### CrypticLoad [Module] : Clustering gene loads 
            # if args.command in ['CrypticLoad'] and args.CLcommand =="Clust":
            #     # gene-level
            #     if CrypticLoad.cluster_load(args.o+args.prefix+".GeneLoad.txt", args.i, args.p, args.o, args.r) == 1:
            #         LogFile.log_message(logfile_path, "Completed clustering gene-level cryptic load : ")
            #         pass
            #     else:
            #         LogFile.log_message(logfile_path, "Failed clustering gene-level cryptic load : ")
            #         LogFile.log_message(logfile_path, "Terminating ............... : ")
            #         exit()









                # get gene-level cryptic load matrix - cryptic counts (num), junction counts total (denom), load (proportion)
                # get sample-level cryptic load
                # ["Diff"]
                    # calculate differential load for gene-level 
                    # calculate differential load for sample-level
                # ["Clust"]
                    # cluster the gene-loads using NMF clustering 
                    # see if you can find groups that separate the treatmeents and memeberships of those







       









        








            

#         #     # Correct for batch effects using passed batch info  #
#         #     if args.batch != None:
#         #         batch = "".join(args.batch).replace(" ","").split(",")
#         #         sjdf = batchCorrect(sjdf,nc,nt,batch)
#         #         if isinstance(sjdf, pd.DataFrame):
#         #             sjdf.to_csv(args.o+args.prefix+"_Junctions_C_A_BatchCorr.txt", sep='\t',header=True, index=False)
#         #             logMessage(logfile, "Completed batch corrections : ")
#         #             pass
#         #         else:
#         #             logMessage(logfile, "Failed batch corrections : ")
#         #             logMessage(logfile, "Terminating ............... : ")
#         #             exit()

#             junction_counts = pd.read_csv(args.o + "/PSI_Annotated_Junction_Counts.txt", sep="\t")

#         #### CrypticJunctions ####
#         if args.command == 'CrypticJunctions': 
#             # BB testing for novel junctions #
#             junction_counts = pd.read_csv(args.o + "/PSI_Annotated_Junction_Counts.txt", sep="\t")


#             novel_junctions=junction_counts[junction_counts['annotation']!="DA"]
#             novel_junctions.to_csv(args.o+"Novel_Junctions.txt", sep="\t", index=None)
#             if run_countdata_bbTest(args.o+"Novel_Junctions.txt", controls, treated, args.p)==1:
#                 logMessage(logfile, "Completed BB test on novel junctions : ")
#             else:
#                 logMessage(logfile, "Failed BB test on novel junctions : ")
#                 logMessage(logfile, "Terminating ............... : ")
#                 exit()

#         #     # PoverAM filter only BB test not for computing load #
#         #     sjdfBB=PoverAM(sjdf,args.o+args.prefix,nc,nt,args.pa_p,args.pa_a,args.pa_m)
#         #     if isinstance(sjdfBB, pd.DataFrame):
#         #         logMessage(logfile, "Completed PoverAM filtering : ")
#         #         sjdfBB.to_csv(args.o+args.prefix+"_PoverA.txt", sep='\t',header=True, index=False)

#         #         pass
#         #     else:
#         #         logMessage(logfile, "Failed PoverAM filtering : ")
#         #         logMessage(logfile, "Terminating ............... : ")
#         #         exit()

#         #     # Calculating junction strength #
#         #     sjdfBB=junction_strength(sjdfBB,nc,nt)
#         #     if isinstance(sjdfBB, pd.DataFrame):
#         #         sjdfBB.to_csv(args.o+args.prefix+"_Junctions_C_A_BB.txt", sep='\t',header=True, index=False)
#         #         logMessage(logfile, "Completed computing junction strength : ")
#         #         pass
#         #     else:
#         #         logMessage(logfile, "Failed computing junction strength : ")
#         #         logMessage(logfile, "Terminating ............... : ")

        
         
#         #     # BB testing for annotated junctions (optional) #
#         #     if args.annotated == 0:
#         #         asjdf=sjdfBB[sjdfBB['anchor']=="DA"]
#         #         asjdf.to_csv(args.o.rstrip("/")+"/"+args.prefix+"_Annotated.txt",sep="\t",header=True,index=False)

#         #     if args.annotated == 1:
#         #         asjdf=sjdfBB[sjdfBB['anchor']=="DA"]
#         #         if runCountdata_BBtest(asjdf,nc,nt,args.o.rstrip("/")+"/"+args.prefix+"_Annotated.txt",args.p):
#         #             logMessage(logfile, "Completed BB test on annotated junctions : ")
#         #         else:
#         #             logMessage(logfile, "Failed BB test on annotated junctions : ")
#         #             logMessage(logfile, "Bypassing ............................ : ")

#         # #### PlotJunctions ####
#         # if args.command =='PlotJunctions':
#         #     logMessage(logfile, "Starting junction plots : ")


#         #     if bam2bw(controls+treated,args.o,args.p):
#         #         logMessage(logfile, "Finished converting bams to bigwigs : ")

#         #         # Make link and bed files #
#         #         if setUp_plots(args.cj, args.o, args.p) == 1:
#         #             logMessage(logfile, "Finished making arc and bed files : ")
#         #             pass
#         #         else:
#         #             logMessage(logfile, "Failed making arc and bed files : ")
#         #             logMessage(logfile, "Skipping plots ............................ : ")

#         #         # Cryptic junction plots #
#         #         if args.t in ['all','cryptic']:
#         #             if cryptic_plots(args.cj,args.n,args.o,nc,nt,args.gtf,args.p,args.ID) == 1:
#         #                 logMessage(logfile, "Finished creating cryptic junction plots : ")
#         #                 pass
#         #             else:
#         #                 logMessage(logfile, "Failed creating cryptic junction plots : ")
#         #                 logMessage(logfile, "Skipping plots ............................ : ")
            

#         #         # Annotated junction plots #
#         #         if args.t in ['all','annotated']:
#         #             if annotated_plots(args.cj,args.n,args.o,nc,nt,args.gtf,args.p,args.ID) == 1:
#         #                 logMessage(logfile, "Finished creating annotated junction plots : ")
#         #                 pass
#         #             else:
#         #                 logMessage(logfile, "Failed creating annotated junction plots : ")
#         #                 logMessage(logfile, "Skipping plots ............................ : ")

#         #     else:
#         #         logMessage(logfile, "Failed converting bams to bigwigs : ")
#         #         logMessage(logfile, "Skipping plots ............................ : ")

        
#         # #### CrypticLoad Diff ####
#         # if args.command =='CrypticLoad' and args.CLcommand =="Diff":
#         #     # Sample-level load #
#         #     if get_sampleLoad(sjdf,controls,treated,args.o.rstrip("/")+"/"+args.prefix,True):
#         #         logMessage(logfile, "Completed computing differential sample loads : ")
#         #     else:
#         #         logMessage(logfile, "Failed computing differential sample loads : ")
#         #         exit()
#         #     # Gene-level load #
#         #     sjdf.to_csv(args.o+"geneLoad_sjdf.txt", sep="\t")

#         #     if get_geneLoad(sjdf,controls,treated,args.o.rstrip("/")+"/"+args.prefix,FileFormat,True,args.p,nc,nt):
#         #         logMessage(logfile, "Completed computing differential gene loads : ")
#         #     else:
#         #         logMessage(logfile, "Failed computing differential gene loads : ")
#         #         exit()
  


#         # #### CrypticLoad Clst ####
#         # if args.command =='CrypticLoad' and args.CLcommand =="Clst":
#         #     # Sample-level load #
#         #     if get_sampleLoad(sjdf,controls,treated,args.o.rstrip("/")+"/"+args.prefix,False):
#         #         logMessage(logfile, "Completed computing sample loads : ")
#         #     else:
#         #         logMessage(logfile, "Failed computing sample loads : ")
#         #         exit()
           
#         #     # Gene-level load #
#         #     if get_geneLoad(sjdf,controls,treated,args.o.rstrip("/")+"/"+args.prefix,FileFormat,False,args.p,nc,nt):
#         #         logMessage(logfile, "Completed computing gene loads : ")

#         #         # NMF Clustering on gene load#
#         #         if clustering_plot(args.o.rstrip("/")+"/"+args.prefix+"_GeneLoad.txt",args.o.rstrip("/")+"/"+args.prefix,args.p,args.i):
#         #             logMessage(logfile, "Completed sample clustering on Gene loads : ")
#         #         else:
#         #             logMessage(logfile, "Failed sample clustering on Gene loads : ")
#         #             exit()
#         #     else:
#         #         logMessage(logfile, "Failed computing Gene loads : ")
#         #         exit()

#         # Tidyup 
#         # del_files=glob.glob(args.o+"*.bed")+glob.glob(args.o+"*.arcs")+glob.glob(args.o+"*.out.tab")+glob.glob(args.o+args.prefix+"_Junctions_*.txt")
#         # log = subprocess.run(['rm']+del_files,stderr=subprocess.DEVNULL,shell=False)    
#         # localdate = time.strftime('%a %m/%d/%Y')
#         # localtime = time.strftime('%H:%M:%S')
#         # logfile.write('# Finished CrypSplice run: '+localdate+' at: ' + localtime+' \n')
#         # logfile.close()
    
    finally:
        print("delete this at the end")
        # #pass
        # # Clean if terminated #
        # del_files=glob.glob(args.o+"*.bed")+glob.glob(args.o+"*.arcs")+glob.glob(args.o+"*.out.tab")+glob.glob(args.o+args.prefix+"_Junctions_*.txt")
        # for file in del_files:
        #   if os.path.exists(file):
        #       subprocess.run(['rm']+[file],stderr=subprocess.DEVNULL,shell=False)
        #       pass

if __name__ == "__main__":
    try:
        import os, sys, glob, time
        import pandas as pd, numpy as np, argparse
        sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
        import Setup, LogFile
        import ExtractJunctions
        import AddGenes
        import AnnotateJunctions
        import DifferentialUsage
        import FilterJunctions
        import CrypticLoad
        import concurrent.futures as cf, subprocess
        import pandas as pd, numpy as np, argparse
        import warnings,logging, pkg_resources
        from collections import defaultdict
        import pysam
        warnings.filterwarnings(action="ignore", category=UserWarning, module="pysam") # for pysam if bam index is older than bam
        warnings.simplefilter(action='ignore', category=FutureWarning) # for gtfparse
        warnings.filterwarnings(action="ignore", category=RuntimeWarning)  # for iNMF
        pd.options.mode.chained_assignment = None  # for pandas in junction_strangth function
        main()
    except ImportError as error:
        print("\nPython 3 or a more recent version is required ...\nTerminating ...\n".format(error.message[16:]))


