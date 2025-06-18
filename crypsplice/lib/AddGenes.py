# AddGenes.py

import sys, os
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import Utilities
import numpy as np
import pandas as pd
import concurrent.futures as cf


# Add Genes
## master function to add gene names/ids to each junction 
def add_genes(annotated_junction_counts, gtf_path, processors, outDir):
	gtf_data = Utilities.read_gtf_file(gtf_path)
	gene_gtf = gtf_data[gtf_data["feature"] == "gene"]
	gene_gtf = gene_gtf.assign(gene_id=gene_gtf['attribute'].str.extract(r'gene_id "(.*?)"'))
	gene_gtf_wNames = gene_gtf.assign(gene_name=gene_gtf['attribute'].str.extract(r'gene_name "(.*?)"'))
	gene_gtf = gene_gtf[["chrom", "start", "end", "strand", "gene_id"]]
	global gene_gtf_array
	gene_gtf_array = np.array(gene_gtf)
	column_names = annotated_junction_counts.columns.tolist()
	column_names.append("gene")
	junctions_array = np.array(annotated_junction_counts)
	split_range = int(round(len(junctions_array) / int(processors), 0))
	split_junctions_array = np.array_split(junctions_array, range(split_range, len(junctions_array), split_range))
	results = []
	with cf.ProcessPoolExecutor(max_workers=int(processors)) as executor:
	    for array in executor.map(prep_junction_genes, split_junctions_array):
	        if len(results) == 0:
	            results = array
	        else:
	            results = np.concatenate((results, array), axis=0)
	junctions_withGenes = pd.DataFrame(results, columns=column_names)
	if junctions_withGenes['gene'].isna().all():
		return(0) 
	else:
		# adding gene names corresponding with gene ids
		gene_mapping = dict(zip(gene_gtf_wNames['gene_id'], gene_gtf_wNames['gene_name']))
		junctions_withGenes = junctions_withGenes.rename(columns={'gene': 'gene_id'})
		junctions_withGenes['gene_name'] = junctions_withGenes['gene_id'].map(lambda x: ','.join(gene_mapping.get(gid.strip(), '') for gid in x.split(',')))
		junctions_withGenes.to_csv(outDir+"Annotated_JunctionCounts.txt", sep="\t", index=None)
		return(1)



# Find Junction Genes
## finding the gene each junction falls within using the gtf coords
def find_junction_genes(junctions_array_elem):
    # junction chrom is the same as gtf gene chrom
    gene_chrom = gene_gtf_array[:, 0] == junctions_array_elem[0]
    # junction strand is the same as gtf gene strand
    gene_strand = gene_gtf_array[:, 3] == junctions_array_elem[3]
    # junction start is in between any gene start and stop coordinates
    gene_juncStarts = ((junctions_array_elem[1] >= gene_gtf_array[:, 1]) & (junctions_array_elem[1] <= gene_gtf_array[:, 2]))
    # junction end is in between any gene start and stop coordinates
    gene_juncEnds = ((junctions_array_elem[2] >= gene_gtf_array[:, 1]) & (junctions_array_elem[2] <= gene_gtf_array[:, 2]))
    # finding genes that the junction may be in
    juncStart_inGene = np.where(gene_chrom & gene_strand & gene_juncStarts)[0]
    juncEnd_inGene = np.where(gene_chrom & gene_strand & gene_juncEnds)[0]
    genes = np.union1d(juncStart_inGene, juncEnd_inGene)
    genes = gene_gtf_array[genes]
    gene_names = ",".join([row[4] for row in genes])
    junctions_array_elem = np.append(junctions_array_elem, gene_names)
    return junctions_array_elem


# Prep Junction Genes
# prep each array element for parallelization
def prep_junction_genes(junctions_array):
    junctions_with_geneInfo = np.apply_along_axis(find_junction_genes, axis=1, arr=junctions_array)
    return junctions_with_geneInfo






