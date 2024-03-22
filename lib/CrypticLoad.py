# CrypticLoad.py


import concurrent.futures as cf
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import numpy as np 
import matplotlib.pyplot as plt



# Calculating gene and sample level cryptic load #


# Get Gene-Load
## finding the gene-level cryptic load [100*(total cryptic junctions in gene / total junctions in gene)]
def get_geneLoad(junction_counts_path, nc, nt, outDir):
    junction_counts = pd.read_csv(junction_counts_path, sep="\t")
    sample_names = junction_counts.columns[4:4+nc+nt]
    # only retaining junctions attributed to a single gene and not NA
    junction_counts = junction_counts.dropna(subset=['gene_ids'])
    junction_counts = junction_counts[~junction_counts['gene_ids'].str.contains(',')]
    # creating gene-load dataframe
    crypCount_cols = [sample + "_cryptic_juncCounts" for sample in sample_names]
    juncCount_cols = [sample + "_total_juncCounts" for sample in sample_names]
    crypLoad_cols = [sample + "_crypticLoad" for sample in sample_names]
    colNames = crypCount_cols + juncCount_cols + crypLoad_cols
    colNames.insert(0, 'gene_ids')
    geneLoad_df = pd.DataFrame(columns = colNames)
    # iterating through each gene in the full junctions df
    grouped = junction_counts.groupby('gene_ids')
    for gene, group in grouped:
        cryptic_counts = group[group["annotation"]!="DA"].iloc[:, 4:(4+nc+nt)]
        junction_counts = group.iloc[:, 4:(4+nc+nt)]
        cryptic_sums = cryptic_counts.sum().values
        junction_sums = junction_counts.sum().values
        # calc gene load and avoid error when dividing Na in np array 
        valid_division = (junction_sums != 0) & ~np.isnan(junction_sums)
        gene_load = np.zeros(nc + nt)
        gene_load[valid_division] = cryptic_sums[valid_division] / junction_sums[valid_division]
        new_row = np.concatenate((cryptic_sums, junction_sums, gene_load)).tolist()
        new_row.insert(0, gene)
        geneLoad_df.loc[len(geneLoad_df)] = new_row
    # saving the dataframe
    if not geneLoad_df.empty:
        geneLoad_df.to_csv(outDir+"GeneLoad.txt", sep="\t", index=None)
        return(1)
    else:
        return(0)



# Get Sample-Load
## finding the sample-level cryptic load [100*(total cryptic junctions in sample / total junctions in sample)]
def get_sampleLoad(junction_counts_path, nc, nt, outDir):
    junction_counts = pd.read_csv(junction_counts_path, sep="\t")
    sample_names = junction_counts.columns[4:4+nc+nt]
    # creating gene-load dataframe
    crypCount_cols = [sample + "_cryptic_juncCounts" for sample in sample_names]
    juncCount_cols = [sample + "_total_juncCounts" for sample in sample_names]
    crypLoad_cols = [sample + "_crypticLoad" for sample in sample_names]
    colNames = crypCount_cols + juncCount_cols + crypLoad_cols
    colNames.insert(0, 'total_junctions')
    sampleLoad_df = pd.DataFrame(columns = colNames)
    # getting totals for all retained junctions
    cryptic_counts = junction_counts[junction_counts["annotation"]!="DA"].iloc[:, 4:(4+nc+nt)]
    junction_counts = junction_counts.iloc[:, 4:(4+nc+nt)]
    cryptic_sums = cryptic_counts.sum().values
    junction_sums = junction_counts.sum().values
    # calc sample load and avoid error when dividing Na in np array 
    valid_division = (junction_sums != 0) & ~np.isnan(junction_sums)
    sample_load = np.zeros(nc + nt)
    sample_load[valid_division] = cryptic_sums[valid_division] / junction_sums[valid_division]
    new_row = np.concatenate((cryptic_sums, junction_sums, sample_load)).tolist()
    # finding total number of junctions in calculations
    new_row.insert(0, len(junction_counts))
    sampleLoad_df.loc[int(len(sampleLoad_df))] = new_row
    # saving the dataframe
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
    gene_names = geneLoad_df["gene_ids"]
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
    

