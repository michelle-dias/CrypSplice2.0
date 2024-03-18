# CrypticLoad.py


import concurrent.futures as cf
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import numpy as np 
import matplotlib.pyplot as plt
import os, sys
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities


# Calculating gene and sample level cryptic load #


# Extract Gene Counts
## extracting counts for the genes which junctions fall into for load calculation 
def extract_gene_counts(junction_counts_path, outDir, sample_paths, strand, processors, gtf_path):    
    junction_counts = pd.read_csv(junction_counts_path, sep="\t",dtype={"chrom": str})
    junction_counts, no_crypJunc = get_gene_sums(junction_counts)
    unique_genes=set(junction_counts["gene"])
    gtf_data=Utilities.read_gtf_file(gtf_path)
    gene_gtf = gtf_data[gtf_data["feature"]=="gene"]
    gene_gtf = gene_gtf.assign(gene_id=gene_gtf['attribute'].str.extract(r'gene_id "(.*?)"'))
    unique_gene_gtf = gene_gtf[gene_gtf["gene_id"].isin(unique_genes)]
    # subsetting gtf for junction genes
    # creating saf file
    saf_file = unique_gene_gtf[["gene_id","chrom","start","end","strand"]]
    saf_file.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
    saf_file_path = outDir+"saf_tmp.txt"
    saf_file.to_csv(saf_file_path, sep="\t", index=None)
    bam_paths=" ".join(sample_paths)
    output_file=outDir+"geneCounts.txt"
    featureCounts_cmd = "featureCounts -p -F 'SAF' -s "+str(strand)+" -T "+str(processors)+" -a "+saf_file_path+" -o "+output_file+" "+bam_paths
    exit_code = os.system(featureCounts_cmd)
    if exit_code==0:
        gene_counts = pd.read_csv(output_file, sep="\t", header=1)
        new_columns = [col.split('/')[-1].replace('.bam', '_geneCounts') for col in gene_counts.columns[6:]]
        gene_counts.columns = list(gene_counts.columns[:6]) + new_columns
        gene_counts.rename(columns={'Geneid': 'gene'}, inplace=True)
        gene_counts = gene_counts[['gene'] + list(gene_counts.columns[gene_counts.columns.get_loc('Length')+1:])]
        merged_df = pd.merge(junction_counts, gene_counts, on='gene', how='left', suffixes=('_junction', '_gene'))
        merged_df.to_csv(outDir+"GeneCounts.txt", sep="\t", index=None)
        os.system("rm "+output_file+" "+output_file+".summary "+saf_file_path)
        return(no_crypJunc)
    else:
        return(0)

# Get Gene-Load
## getting gene-level sums of cryptic junctions - nested in above function
def get_gene_sums(junction_counts):
    # filtering for cryptic counts 
    junction_counts = junction_counts[junction_counts["annotation"]!="DA"]
    no_crypJunc = len(junction_counts)
    cols_juncCounts = [col for col in junction_counts.columns if col.endswith('juncCounts')]
    cryptic_junctions = junction_counts[cols_juncCounts + ['gene']]
    crypJunc_sums = cryptic_junctions.groupby('gene').sum()
    crypJunc_sums["gene"] = crypJunc_sums.index
    crypJunc_sums.reset_index(drop=True, inplace=True)
    crypJunc_sums = crypJunc_sums[['gene'] + cols_juncCounts]
    return(crypJunc_sums, no_crypJunc)



# Calculate Gene-Load
## finding the gene-level cryptic load [100*(total cryptic junctions in gene / total junctions in gene)]
def calculate_geneLoad(gene_sums_path, outDir):
    geneLoad_df = pd.read_csv(gene_sums_path, sep="\t", dtype={"chrom": str})
    junc_columns = [col for col in geneLoad_df.columns if col.endswith('juncCounts')]
    gene_columns = [col for col in geneLoad_df.columns if col.endswith('geneCounts')]
    sample_prefixes = [col.replace('_juncCounts', '') for col in junc_columns]
    for prefix in sample_prefixes:
        geneLoad_df[prefix+"_crypticLoad"]=geneLoad_df[prefix+"_juncCounts"]/geneLoad_df[prefix+"_geneCounts"]
    load_columns = [col for col in geneLoad_df.columns if col.endswith('crypticLoad')]
    # replace all inf with 0 - gene counts were 0 -> 1/0=inf
    # replace all blanks with 0
    geneLoad_df.replace([np.inf, -np.inf, '', 'NaN'], 0, inplace=True)
    geneLoad_df.fillna(0, inplace=True)
    if not geneLoad_df.empty:
        geneLoad_df.to_csv(outDir+"GeneLoad.txt", sep="\t", index=None)
        return(1)
    else:
        return(0)




# Calculate Sample-Load
## finding the sample-level cryptic load [100*(total cryptic junctions in sample / total junctions in sample)]
def calculate_sampleLoad(gene_sums_path,outDir,n_crypticJunctions):
    junction_counts = pd.read_csv(gene_sums_path, sep="\t")
    # sample_names = junction_counts.columns[4:4+nc+nt]
    # # creating gene-load dataframe
    junc_columns = [col for col in junction_counts.columns if col.endswith('juncCounts')]
    gene_columns = [col for col in junction_counts.columns if col.endswith('geneCounts')]
    load_columns = [col.replace('_juncCounts', '_crypticLoad') for col in junc_columns]
    sample_prefixes = [col.replace('_juncCounts', '') for col in junc_columns]
    colNames = junc_columns + gene_columns + load_columns
    colNames.insert(0, 'total_junctions')
    sampleLoad_df = pd.DataFrame(columns = colNames)
    # getting totals for all retained junctions
    junc_sums = junction_counts[junc_columns].sum().values
    gene_sums = junction_counts[gene_columns].sum().values
    sample_loads = [0] * len(junc_columns)
    data = [[0] + junc_sums.tolist() + gene_sums.tolist() + sample_loads]
    sampleLoad_df = pd.DataFrame(data, columns=colNames)
    for prefix in sample_prefixes:
        sampleLoad_df[prefix+"_crypticLoad"] = sampleLoad_df[prefix+"_juncCounts"]/sampleLoad_df[prefix+"_geneCounts"]
    # finding total number of junctions in calculations
    sampleLoad_df["total_junctions"]=n_crypticJunctions
    if not sampleLoad_df.empty:
        sampleLoad_df.to_csv(outDir+"SampleLoad.txt", sep="\t", index=None)
        return(1)
    else:
        return(0)
	



# NMF Iterations
## do multiple iterations of NMF 
def NMF_iterations(nmf_df):
    rank_range = range(2, len(nmf_df.columns))
    silhouette_iter = []
    davies_bouldin_iter = []
    calinski_harabasz_iter = []
    reconstruction_error_iter = []
    rank_range = range(2, len(nmf_df.columns) + 1)
    for rank in rank_range:
        nmf = NMF(n_components=rank, max_iter=100000, init='random')
        W = nmf.fit_transform(nmf_df)
        labels = W.argmax(axis=1)
        silhouette_iter.append(silhouette_score(nmf_df, labels))
        davies_bouldin_iter.append(davies_bouldin_score(nmf_df, labels))
        calinski_harabasz_iter.append(calinski_harabasz_score(nmf_df, labels))
        reconstruction_error_iter.append(nmf.reconstruction_err_)
    return silhouette_iter, davies_bouldin_iter, calinski_harabasz_iter, reconstruction_error_iter


# Plot Rank Metrics
## plot the silhouette, davies_bouldin, calinski_harabasz, and reconstruction error over all ranks (clusters)
def plot_rank_metrics(silhouette_rank,davies_bouldin_rank,calinski_harabasz_rank,recon_error_rank, outdir):
    # Create a 2x2 panel of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 8),facecolor='white')
    fig.suptitle('Rank Evaluation Metrics', fontsize=16)
    # Plot 1: Silhouette
    axs[0, 0].plot(silhouette_rank, marker='o', linestyle='-', color='blue')
    axs[0, 0].set_title('Silhouette')
    # Plot 2: Davies-Bouldin
    axs[0, 1].plot(davies_bouldin_rank, marker='o', linestyle='-', color='green')
    axs[0, 1].set_title('Davies-Bouldin')
    # Plot 3: Calinski-Harabasz
    axs[1, 0].plot(calinski_harabasz_rank, marker='o', linestyle='-', color='orange')
    axs[1, 0].set_title('Calinski-Harabasz')
    # Plot 4: Reconstruction Error
    axs[1, 1].plot(recon_error_rank, marker='o', linestyle='-', color='red')
    axs[1, 1].set_title('Reconstruction Error')
    # Set common labels
    for ax in axs.flat:
        ax.set_ylabel('Score')
        ax.set_xlabel('Cluster Rank')
        ax.set_xticks(np.arange(0, len(calinski_harabasz_rank)))
        ax.set_xticklabels(np.arange(2, len(calinski_harabasz_rank) + 2))
    plt.tight_layout()
    # Save the plot
    plt.savefig(outdir+'RankEvalMetrics.png')



# Find Cluster Memberships
## create dataframes to show membership of each cluster  
def find_cluster_membership(rank, nmf_df, outDir, geneLoad_df):
    np.random.seed(2024)
    # Initialize NMF with the specified rank
    nmf = NMF(n_components=rank, max_iter=1000000, init='random', random_state=2024)
    # Fit and transform the data to obtain W matrix
    W = nmf.fit_transform(nmf_df)
    # Obtain the H matrix
    H = nmf.components_
    # Creating Sample Membership df 
    sample_names = nmf_df.columns
    sample_names = [name.replace("_crypticLoad", "") for name in sample_names]
    H_mat = pd.DataFrame(H, columns=sample_names)
    max_indices = H_mat.idxmax() + 1  # Adding 1 to convert from zero-based index to 1-based index
    sample_membership = pd.DataFrame({
        'Sample': max_indices.index,
        'Cluster': max_indices.values,
        'Coefficient': H_mat.max()
    })
    sample_membership = sample_membership.reset_index(drop=True)
    # Creating Gene Membership DataFrame
    gene_names = geneLoad_df["gene"]
    W_mat = pd.DataFrame(W, index=gene_names)
    # creating a new cluster for rows with all 0 -> could not be properly assigned to a cluster
    W_mat[rank] = (W_mat == 0).all(axis=1).astype(int)
    max_columns = W_mat.idxmax(axis=1)
    max_values = W_mat.max(axis=1)
    gene_membership = pd.DataFrame({
        'Gene': W_mat.index,
        'Cluster': max_columns + 1,  # Adding 1 to convert from zero-based index to 1-based index
        'Coefficient': max_values
    })
    # if cluster is rank+1 then put NA in cluster column
    gene_membership['Cluster'] = np.where(gene_membership['Cluster'] == rank + 1, np.nan, gene_membership['Cluster'])
    gene_membership = gene_membership.sort_values(by=['Cluster', 'Coefficient'], ascending=[True, False])
    gene_membership = gene_membership.reset_index(drop=True)
    return sample_membership, gene_membership, H_mat
    
    
# Cluster Load
## master function to select rank, generate rank plot, and cluster the loads
def cluster_load(geneLoad_path, iterations, processors, outdir, clusters):
    geneLoad_df = pd.read_csv(geneLoad_path, sep="\t")
    nmf_df = geneLoad_df.filter(regex='crypticLoad$', axis=1)
    if clusters:
        rank= clusters
    else:
        pl=[]
        for i in range(0,iterations):
            pl.append(nmf_df)
        with cf.ProcessPoolExecutor(max_workers=int(processors)) as (executor):
            metrics = list(executor.map(NMF_iterations, pl))
        metrics_avg = np.mean(np.array(metrics), axis=0)
        silhouette_rank = metrics_avg[0]
        davies_bouldin_rank = metrics_avg[1]
        calinski_harabasz_rank = metrics_avg[2]
        recon_error_rank = metrics_avg[3]
        # plotting the metrics
        plot_rank_metrics(silhouette_rank, davies_bouldin_rank, calinski_harabasz_rank, recon_error_rank, outdir)
        # automating rank selection by selecting the highest of silhouette, davies-bouldin, and calinski-harabasz (if not given by user)
        max_index_silhouette = np.argmax(silhouette_rank)+2
        min_index_davies_bouldin = np.argmin(davies_bouldin_rank)+2
        max_index_calinski_harabasz = np.argmax(calinski_harabasz_rank)+2
        indices = [max_index_silhouette, min_index_davies_bouldin, max_index_calinski_harabasz]
        rank = max(set(indices), key=lambda x: indices.count(x))
    # sample and gene membership matrices 
    sample_membership, gene_membership, H_mat = find_cluster_membership(rank, nmf_df, outdir, geneLoad_df)
    if sample_membership.empty and gene_membership.empty:
        return(0)
    else:
        # save entire H matrix to better be able to interpret clusters for now  
        H_mat.to_csv(outdir+"/Hmatrix_Sample_ClusterMembership.txt", sep="\t", index=True)
        gene_membership.to_csv(outdir+"/Gene_ClusterMembership.txt", sep="\t", index=False)
        sample_membership.to_csv(outdir+"/Sample_ClusterMembership.txt", sep="\t", index=False)
        return(1)
    

