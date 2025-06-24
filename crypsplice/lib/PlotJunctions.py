# PlotJunctions.py
# bamCoverage, wiggletools need to be on the path 
# wigToBigWig utility needs to be in visualization imports file 


import sys, os
import glob
import pandas as pd
import pyBigWig

sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib/visualization_imports/")




# Bam to BigWig
## master function for converting bams to bigWigs and merging bigWigs if enabled by user
def bam_to_bigWig(control_bams, treated_bams, fasta_file, outDir, processors, merge):
    #converting all bams to bigWigs
    if run_bamCoverage(control_bams+treated_bams,outDir,processors) == 1:
        if merge==True:
            if create_genomeSizes(fasta_file, outDir) == 1:
                return_val = merge_bigWigs(control_bams, treated_bams, outDir)
            else:
                return_val = 0
        else:
            return_val = 1
    return(return_val)

# Run bamCoverage
## running bamCoverage to convert bams to bigWigs
def run_bamCoverage(bfiles,outDir,processors):
    try:
        for bfile in bfiles:
            command = f"bamCoverage -b {bfile} -bs 20 -p {processors} --normalizeUsing None --skipNonCoveredRegions --smoothLength 60 --centerReads -o {outDir}{bfile.split('/')[-1].replace('.bam', '.bw')}"
            os.system(command)
    except:
        return(0)
    return(1)

# Create Genome Sizes
## creating the genome sizes file from gtf, it's necessary to merge bigWigs
def create_genomeSizes(fasta_file, outDir):
    # setting output file name
    output_file=outDir+"genome_sizes.txt"
    # Create genome sizes file
    cmd="awk '/^>/ { if (seqlen) { print seqname, seqlen } seqname = substr($1, 2); seqlen = 0; next } { seqlen += length($0) } END { print seqname, seqlen }' "+fasta_file+" > "+output_file
    os.system(cmd)
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        return(1)
    else:
        return(0)

# Merge BigWigs
## merging the C1 and C2 bigWigs
def merge_bigWigs(c1_bams, c2_bams, outDir):
    # getting bigWig paths
    c1_samples = [bam.split("/")[-1].replace(".bam","") for bam in c1_bams]
    c2_samples = [bam.split("/")[-1].replace(".bam","") for bam in c2_bams]
    c1_bigWigs = [f"{outDir}{sample}.bw" for sample in c1_samples]
    c2_bigWigs = [f"{outDir}{sample}.bw" for sample in c2_samples]
    # combining C1 and C2 files with wiggletools (mean)
    c1_wig_path = outDir+"C1.wig"
    c2_wig_path = outDir+"C2.wig"
    c1_wiggletools_cmd = "wiggletools write " + c1_wig_path + " mean " + " ".join(c1_bigWigs)
    os.system(c1_wiggletools_cmd)
    c2_wiggletools_cmd = "wiggletools write " + c2_wig_path + " mean " + " ".join(c2_bigWigs)
    os.system(c2_wiggletools_cmd)
    # converting wig files to bigWigs
    genomeSizes_file = outDir+"genome_sizes.txt"
    c1_wig2bw_cmd = "wigToBigWig "+c1_wig_path+" "+genomeSizes_file+" "+outDir+"C1.bw"
    os.system(c1_wig2bw_cmd)
    c2_wig2bw_cmd = "wigToBigWig "+c2_wig_path+" "+genomeSizes_file+" "+outDir+"C2.bw"
    os.system(c2_wig2bw_cmd)
    # check that the bigWig files are there and not empty
    for bigwig_file in [f"{outDir}C1.bw", f"{outDir}C2.bw"]:
        if os.path.isfile(bigwig_file) and os.path.getsize(bigwig_file) > 0:
            return_val = 1
        else:
            return_val = 0
    # removing wig files
    os.system("rm "+c1_wig_path+" "+c2_wig_path)
    return(return_val)
    





# Create Sashimis
## master function for creating sashimi file input to pygenometracks for each sample or the merged samples 
def create_sashimis(nc, nt, outDir, merge):
    # get annotated and cryptic junction data
    annotated_junctions_path = glob.glob(outDir+"*Annotated_JunctionCounts.txt")[0]
    annotated_junctions = pd.read_csv(annotated_junctions_path, sep="\t")
    junction_data = annotated_junctions.iloc[:,:4+nc+nt]
    junction_data = pd.concat([junction_data, annotated_junctions['annotation']], axis=1)
    # splitting cryptic and annotated data
    cryptic_data = junction_data[junction_data["annotation"]!="DA"]
    annotated_data = junction_data[junction_data["annotation"]=="DA"]
    # if merge is enabled
    if merge == True:
        return_val = create_merged_sashimis(cryptic_data, annotated_data, nc, nt, outDir)
    else:
        return_val = create_indiv_sashimis(cryptic_data, annotated_data, nc, nt, outDir)
    return(return_val)


# Create Merged Sashimis
## creating sashimi file for the merged bigWigs
def create_merged_sashimis(cryptic_data, annotated_data, nc, nt, outDir):
    # creating sashimi for merged samples
    # getting C1 and C2 junction averages
    cryptic_data["C1"]=cryptic_data.iloc[:, 4:nc+4].mean(axis=1)
    cryptic_data["C2"]=cryptic_data.iloc[:, nc+4:nc+nt+4].mean(axis=1)
    annotated_data["C1"]=annotated_data.iloc[:, 4:nc+4].mean(axis=1)
    annotated_data["C2"]=annotated_data.iloc[:, nc+4:nc+nt+4].mean(axis=1)
    # formatting C1 and C2 sashimi file
    for C_id in ["C1", "C2"]:
        cryptic_sashimi_data = cryptic_data[["chrom", "start", "start", "chrom", "end", "end", C_id]].copy()
        cryptic_sashimi_data["annotation"] = "cryptic"
        annotated_sashimi_data = annotated_data[["chrom", "start", "start", "chrom", "end", "end", C_id]].copy()
        annotated_sashimi_data["annotation"] = "annotated"
        sashimi_data = pd.concat([cryptic_sashimi_data, annotated_sashimi_data], ignore_index=True)
        sashimi_data.columns = ["chrom", "start", "start_2", "chrom_2", "end", "end_2", "counts", "annotation"]
        sashimi_data = sashimi_data.sort_values(by=['chrom', 'start', 'end'], ascending=[True, True, True])
        sashimi_data = sashimi_data[sashimi_data['counts'] != 0]
        #saving sashimi file
        file_path = outDir+C_id+".sashimi"
        sashimi_data.to_csv(file_path, sep="\t", index=None, header=None)
        if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
            return_val = 1
        else:
            return_val = 0
    return(return_val)

# Create Individual Sashimis
## creating sashimi file for the individual bigWigs
def create_indiv_sashimis(cryptic_data, annotated_data, nc, nt, outDir):
    junc_columns = cryptic_data.columns[4:(nc+nt+4)]
    for junc_col in junc_columns:
        # getting sample name
        sample_id = junc_col.split("_juncCounts")[0]
        # formatting sashimi file
        cryptic_sashimi_data = cryptic_data[["chrom", "start", "start", "chrom", "end", "end", junc_col]].copy()
        cryptic_sashimi_data["annotation"] = "cryptic"
        annotated_sashimi_data = annotated_data[["chrom", "start", "start", "chrom", "end", "end", junc_col]].copy()
        annotated_sashimi_data["annotation"] = "annotated"
        sashimi_data = pd.concat([cryptic_sashimi_data, annotated_sashimi_data], ignore_index=True)
        sashimi_data.columns = ["chrom", "start", "start_2", "chrom_2", "end", "end_2", "counts", "annotation"]
        sashimi_data = sashimi_data.sort_values(by=['chrom', 'start', 'end'], ascending=[True, True, True])
        sashimi_data = sashimi_data[sashimi_data['counts'] != 0]
        #saving sashimi file
        file_path=outDir+sample_id+".sashimi"
        sashimi_data.to_csv(file_path, sep="\t", index=None, header=None)
        if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
            return_val = 1
        else:
            return_val = 0
    return(return_val)





# Create Plots
## master function to generate pyGenome configs and plots
def create_plots(control_bams, treated_bams, gtf, outDir, merge, pj_file, junc_ids, top_n, c1_color, c2_color, c1_names, c2_names):
    plot_params = setup_config_inputs(control_bams, treated_bams, gtf, outDir, merge, pj_file, junc_ids, top_n, c1_color, c2_color, c1_names, c2_names)
    if plot_params == "no_sig":
        return_val = 2
    else:
        for params in plot_params:
            plot_junction(params)
            # check that plot output is not empty and delete temp file
            if os.path.isfile(params[15]) and os.path.getsize(params[15]) > 0:
                return_val = 1
            else:
                return_val = 0
            # delete temp directory
            cmd="rm -r "+params[16]
            os.system(cmd)
    return(return_val)

    

# Setup Config Inputs
## setting up inputs for configurations of each plot 
def setup_config_inputs(control_bams, treated_bams, gtf, outDir, merge, pj_file, junc_ids, top_n, c1_color, c2_color, c1_names, c2_names):
    # getting paths of control and treated bigWig files
    if merge == True:
        c1_bigWigs = [outDir+"C1.bw"]
        c2_bigWigs = [outDir+"C2.bw"]
    else:
        control_samples = [path.rsplit("/", 1)[1].replace(".bam", "") for path in control_bams]
        c1_bigWigs = [f"{outDir}{sample}.bw" for sample in control_samples]
        treated_samples = [path.rsplit("/", 1)[1].replace(".bam", "") for path in treated_bams]
        c2_bigWigs = [f"{outDir}{sample}.bw" for sample in treated_samples]

    # getting the c1 and c2 names 
    if c1_names and c2_names is not None:
        c1_names = c1_names.split(",")
        c2_names = c2_names.split(",")
    elif merge == True:
        c1_names = ["C1"]
        c2_names = ["C2"]
    else:
        c1_names = control_samples
        c2_names = treated_samples

    # getting junction ids of interest
    # reading in the novel junctions
    junctions = pd.read_csv(pj_file, sep="\t")
    # if the junction_id is not provided select the top n to plot
    if junc_ids == None:
        # sort the novel junctions file to get the top n junction_ids
        sig_cryptic_juncs = junctions[junctions["adj.pVal"] <= 0.05]
        if sig_cryptic_juncs.empty:
            return("no_sig")
        else:
            sig_cryptic_juncs = sig_cryptic_juncs.copy()
            sig_cryptic_juncs["CS_diff_abs"] = sig_cryptic_juncs["CS_diff"].abs()
            sig_cryptic_juncs_sort = sig_cryptic_juncs.sort_values(by='CS_diff_abs', ascending=False)
            junc_ids = sig_cryptic_juncs_sort['juncID'].head(top_n).tolist()
    # if the junction_ids are provided - parse the comma separated list
    else:
        junc_ids = junc_ids.split(",")
    # getting the junction specific parameters
    junc_plot_params = list()
    for junc_id in junc_ids:
        # parsing the junction id
        chrom = junc_id.split(":")[0]
        start = junc_id.split(":")[1].split("-")[0]
        end = junc_id.split(":")[1].split("-")[1].split("(")[0]
        strand = junc_id.split("(")[1].split(")")[0]
        # creating temp sashimi files with tagged junction of interest
        c1_sashimis, c2_sashimis, tmp_dir = create_temp_sashimis(outDir, chrom, start, end, merge, control_bams, treated_bams)
        # setting the region (not strand specific)
        region = chrom+":"+start+"-"+end
        # setting the data range for the bigWigs
        data_range=0
        for bigWig in c1_bigWigs+c2_bigWigs:
            bw = pyBigWig.open(bigWig)
            stats = bw.stats(chrom, int(start), int(end), type="max")
            bw.close()
            if stats[0] >= data_range:
                data_range=stats[0]
        data_range = data_range + 300
        # setting region to plot
        ## adding 1kb on either side
        plot_start = int(start)-1000
        plot_end = int(end)+1000
        region=chrom+":"+str(plot_start)+"-"+str(plot_end)
        # setting the number of bins in the plot
        bin_num = 2*(plot_end - plot_start)
        # getting the gene name
        gene_name = junctions[junctions["juncID"]==junc_id]["gene_name"].iloc[0]
        # setting output paths
        config_output = outDir+chrom+"_"+start+"_"+end+"_"+strand+".config.ini"
        plot_output = outDir+chrom+"_"+start+"_"+end+"_"+strand+".pdf"
        # packaging all plot params together
        params = [junc_id, c1_bigWigs, c2_bigWigs, c1_sashimis, c2_sashimis, c1_names, c2_names, c1_color, c2_color, gtf, region, data_range, bin_num, gene_name, config_output, plot_output, tmp_dir] 
        junc_plot_params.append(params)
    return(junc_plot_params)



# Create Temp Sashimis
# creating a temporary sashimi file with the junction of interest identified 
def create_temp_sashimis(outDir, chrom, start, end, merge, control_bams, treated_bams):
    # create temporary directory for new sashimi storage
    tmp_dir=outDir+chrom+"_"+str(start)+"_"+str(end)+"_tmp/"
    cmd="mkdir "+tmp_dir
    os.system(cmd)
    # grab sashimi files and alter annotation column of select cryptic junction
    sashimi_files = glob.glob(outDir+"*.sashimi")
    # identifying cryptic junction of interest in sashimi file
    for sashimi_path in sashimi_files:
        sashimi_data = pd.read_csv(sashimi_path, sep="\t", header=None)
        matching_rows = sashimi_data[((sashimi_data[0] == chrom) & (sashimi_data[1] == int(start)) & (sashimi_data[4] == int(end)))]
        if matching_rows.empty:
            new_row = [chrom, start, start, chrom, end, end, 0, "select"]
            sashimi_data.loc[len(sashimi_data)] = new_row
            sashimi_data = sashimi_data.sort_values(by=[0, 1, 4], ascending=[True, True, True])
        else:
            matching_rows = sashimi_data[((sashimi_data[0] == chrom) & (sashimi_data[1] == int(start)) & (sashimi_data[4] == int(end)))]                                        
            sashimi_data.loc[matching_rows.index[0], 7] = 'select'        
        sashimi_data.to_csv(tmp_dir+sashimi_path.split("/")[-1], sep="\t", header=None, index=None)
    c1_sashimis = []
    c2_sashimis = []
    if merge==True:
        c1_sashimis.append(glob.glob(tmp_dir+"*C1*")[0])
        c2_sashimis.append(glob.glob(tmp_dir+"*C2*")[0])
    else:
        control_samples = [path.rsplit("/", 1)[1].replace(".bam", "") for path in control_bams]
        for sample in control_samples:
            sashimi_file=glob.glob(tmp_dir+"*"+sample+"*")[0]
            c1_sashimis.append(sashimi_file)
        treated_samples = [path.rsplit("/", 1)[1].replace(".bam", "") for path in treated_bams]
        for sample in treated_samples:
            sashimi_file=glob.glob(tmp_dir+"*"+sample+"*")[0]
            c2_sashimis.append(sashimi_file)
    return(c1_sashimis, c2_sashimis, tmp_dir)

# Plot Junction
## writing the configuration file and plotting each junctions with pyGenomeTracks
def plot_junction(params):
    junc_id = params[0]
    c1_bigWigs = params[1]
    c2_bigWigs = params[2]
    c1_sashimis = params[3]
    c2_sashimis = params[4]
    c1_names = params[5]
    c2_names = params[6]
    c1_color = params[7]
    c2_color = params[8]
    gtf = params[9]
    region = params[10]
    data_range = params[11]
    bin_num = params[12]
    gene_name = params[13]
    config_output= params[14]
    plot_output=params[15]

    # initializing config file
    with open(config_output, 'w') as file:
        file.write("[spacer]\n")
        file.write("height = 0.25\n\n")
        file.write("[scale]\n")
        file.write("title = \nheight = 0.1\nwhere = right\nfontsize = 4\nfile_type = scalebar\n\n")
        file.write("[spacer]\nheight = 0.25\n\n")

    # writing in track info for each C1 bigWig
    for i in range(len(c1_bigWigs)):
        bigWig_path = c1_bigWigs[i]
        sashimi_path = c1_sashimis[i]
        name = c1_names[i]
        sample_id = bigWig_path.split("/")[-1].replace(".bw", "") 
        with open(config_output, "a") as file:
            file.write("["+sample_id+"]\n")
            file.write("file = "+bigWig_path+"\n")
            file.write("link_file = "+sashimi_path+"\n")
            file.write("title = "+name+"\n")
            file.write("height = 3\n")
            file.write("bw_color = "+c1_color+"\n")
            file.write("link_color = "+c1_color+"\n")
            file.write("line_style = solid\n")
            file.write("number_of_bins = "+str(bin_num)+"\n")
            file.write("nans_to_zero = true\n")
            file.write("summary_method = mean\n")
            file.write("show_data_range = true\n")
            file.write("max_value = "+str(data_range)+"\n")
            file.write("min_value = 0\n")
            file.write("fontsize = 6\n")
            file.write("show_number = true\n")
            file.write("file_type = sashimiBigWig\n\n")
            file.write("[spacer]\nheight = 0.3\n\n")

    # writing in track info for each C2 bigWig
    for i in range(len(c2_bigWigs)):
        bigWig_path = c2_bigWigs[i]
        sashimi_path = c2_sashimis[i]
        name = c2_names[i]    
        sample_id = bigWig_path.split("/")[-1].replace(".bw", "") 
        with open(config_output, "a") as file:
            file.write("["+sample_id+"]\n")
            file.write("file = "+bigWig_path+"\n")
            file.write("link_file = "+sashimi_path+"\n")
            file.write("title = "+name+"\n")
            file.write("height = 3\n")
            file.write("bw_color = "+c2_color+"\n")
            file.write("link_color = "+c2_color+"\n")
            file.write("line_style = solid\n")
            file.write("number_of_bins = "+str(bin_num)+"\n")
            file.write("nans_to_zero = true\n")
            file.write("summary_method = mean\n")
            file.write("show_data_range = true\n")
            file.write("max_value = "+str(data_range)+"\n")
            file.write("min_value = 0\n")
            file.write("fontsize = 6\n")
            file.write("show_number = true\n")
            file.write("file_type = sashimiBigWig\n\n")
            file.write("[spacer]\nheight = 0.3\n\n")

    # writing gtf info 
    with open(config_output, "a") as file:
        file.write("[Genes]\n")
        file.write("file = "+gtf+"\n")
        file.write("height = 1.5\n")
        file.write("merge_transcripts = true\n")
        file.write("prefered_name = gene_name\n")
        file.write("fontsize = 0\n")
        file.write("arrowhead_included = true\n")
        file.write("arrow_interval = 10\n")
        file.write("color = gray\n")
        file.write("labels = on\n")
        file.write("file_type = gtf\n")
        file.write("style = UCSC\n")
        file.write("gene_rows = 1\n\n")
        file.write("[x-axis]\nfontsize = 6\nwhere = bottom\n\n")
        file.write("[spacer]\nheight = 0.05\n")

    #pygenome command
    title = "'"+str(junc_id) +"  ["+str(gene_name)+"]"+"'"
    plot_cmd = "pyGenomeTracks --tracks "+config_output+" --region "+region+" -t "+title+" --width 24 --trackLabelFraction 0.01 -out "+plot_output+" --fontSize 6"
    os.system(plot_cmd)

