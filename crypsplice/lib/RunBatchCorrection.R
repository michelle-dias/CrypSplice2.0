



#devtools::install_github("zhangyuqing/sva-devel")
library(sva)

# Receiving arguments from python
args <- commandArgs(trailingOnly=TRUE)

# Capturing individual components of arguments
junctions_df_path <- args[1]
#junctions_df_path <- "/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypSplice_Editing/Batch_Corrections_4_2_25/BC_test_1/CrypSplice.JunctionCounts.txt"
batchCorr_meta_path <- args[2]
#batchCorr_meta_path <- "/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypSplice_Editing/Batch_Corrections_4_2_25/BC_test_1/batchCorr_meta.txt"



junctions_df <- read.table(junctions_df_path, sep="\t", header=TRUE)

batchCorr_meta <- read.table(batchCorr_meta_path, sep="\t", header=TRUE)


# Read the data (assuming it's in a dataframe called df)



# subsetting jucntion counts columns
juncCounts_cols <- grep("juncCounts$", colnames(junctions_df))

# extract junction counts and batch correct
juncCounts <- junctions_df[, juncCounts_cols]
juncCounts_order <- gsub("_juncCounts$", "", colnames(juncCounts))
batchCorr_meta$Sample <- factor(batchCorr_meta$Sample, levels = juncCounts_order, ordered = TRUE)
batches <- batchCorr_meta$Batch
juncCounts_bc <- ComBat_seq(juncCounts, batch=batches, group=NULL)

# replacing original junction counts with batch corrected counts
junctions_df[, juncCounts_cols] <- juncCounts_bc

# subsetting jucntion counts columns
originCounts_cols <- grep("originCounts$", colnames(junctions_df))

# extract junction counts and batch correct
originCounts <- junctions_df[, originCounts_cols]
originCounts_order <- gsub("_juncCounts$", "", colnames(originCounts))
batchCorr_meta$Sample <- factor(batchCorr_meta$Sample, levels = originCounts_order, ordered = TRUE)
batches <- batchCorr_meta$Batch
originCounts_bc <- ComBat_seq(originCounts, batch=batches, group=NULL)

# replacing original junction counts with batch corrected counts
junctions_df[, originCounts_cols] <- originCounts_bc



# saving batch corrected data to path if not empty and returning 1 
if (!all(is.na(junctions_df[, originCounts_cols])) && !all(is.na(junctions_df[, juncCounts_cols]))) {
  write.table(junctions_df, junctions_df_path, sep="\t", row.names = FALSE)
  cat(1)
} else {
  cat(0)
}




