if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("sva", quietly=TRUE)) BiocManager::install("sva")

library(sva)

args <- commandArgs(trailingOnly=TRUE)
junctions_df_path   <- args[1]
#junctions_df_path <- "/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypticLoad_ROSMAP/Attempt_3_Batch/CrypticLoad_Clust/CrypSplice.JunctionCounts.txt"

batchCorr_meta_path <- args[2]
#batchCorr_meta_path <- "/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypticLoad_ROSMAP/Attempt_3_Batch/batchMeta_test.txt"


junctions_df   <- read.table(junctions_df_path,  sep="\t", header=TRUE, check.names=FALSE)
batchCorr_meta <- read.table(batchCorr_meta_path, sep="\t", header=TRUE, check.names=FALSE)

# --- Helper: extract, batch-correct, and return corrected matrix ---
batch_correct <- function(df, col_pattern, strip_pattern, meta) {
  cols   <- grep(col_pattern, colnames(df))
  counts <- df[, cols, drop=FALSE]
  order  <- gsub(strip_pattern, "", colnames(counts))
  
  # Reorder metadata to match column order of count matrix
  meta$Sample <- as.character(meta$Sample)
  meta <- meta[match(order, meta$Sample), ]
  
  if (any(is.na(meta$Sample))) {
    stop(paste("Sample mismatch! These samples are in counts but not in metadata:",
               paste(order[!order %in% meta$Sample], collapse=", ")))
  }
  
  batches   <- meta$Batch
  corrected <- ComBat_seq(as.matrix(counts), batch=batches, group=NULL)
  return(list(corrected=corrected, cols=cols))
}

# Batch correct junction counts
junc_result   <- batch_correct(junctions_df, "juncCounts$",   "_juncCounts$",   batchCorr_meta)
junctions_df[, junc_result$cols] <- junc_result$corrected

# Batch correct origin counts
origin_result <- batch_correct(junctions_df, "originCounts$", "_originCounts$", batchCorr_meta)
junctions_df[, origin_result$cols] <- origin_result$corrected

# Save and return status
junc_cols   <- junc_result$cols
origin_cols <- origin_result$cols

# Replace the final if/else block with this:
if (!all(is.na(junctions_df[, origin_cols])) && !all(is.na(junctions_df[, junc_cols]))) {
  write.table(junctions_df, junctions_df_path, sep="\t", row.names=FALSE)
  cat("1\n")
} else {
  cat("0\n")
}

