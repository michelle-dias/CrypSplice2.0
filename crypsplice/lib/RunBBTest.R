# RunBBTest.R

library(countdata)

# Receiving arguments from python
args <- commandArgs(trailingOnly=TRUE)

# Capturing individual components of arguments
junctions_df_path <- args[1]
n_controls <- as.integer(args[2])
n_treated <- as.integer(args[3])
processors <- as.integer(args[4])
## CrypticJunctions or CrypticLoad
test_type <- args[5]



junctions_df <- read.csv(junctions_df_path, sep="\t")

if (test_type=="CrypticLoad") {
  # sub-setting the matrix to get the junction counts (numerator) and origin counts (denominator)
  numerator_mat <- as.matrix(junctions_df[2:(1+n_controls+n_treated)])
  denominator_mat <- as.matrix(junctions_df[(2+n_controls+n_treated):(1+2*(n_controls+n_treated))])
} else { 
  # sub-setting the matrix to get the junction counts (numerator) and origin counts (denominator)
  numerator_mat <- as.matrix(junctions_df[5:(4+n_controls+n_treated)])
  denominator_mat <- as.matrix(junctions_df[(5+n_controls+n_treated):(4+2*(n_controls+n_treated))])
  }


# adding 0.0001 to the denominator to pad 
## bb.test does not accept denominator=0 or numerator=denominator 
denominator_mat <- denominator_mat + 0.1


# creating the grouping classifications 
groups <- c(rep("control", n_controls), rep("treated", n_treated))

# running the countdata bb test
pValues <- countdata::bb.test(numerator_mat, denominator_mat, groups, n.threads = processors)
adj_pValues <- p.adjust(pValues$p.value, method = "BH")

# adding p-values and adjusted p-values to the annotated junctions dataframe 
junctions_df["pVal"] <- pValues$p.value
junctions_df["adj.pVal"] <- adj_pValues


# checking that the p-value columns are not empty 
if (!all(is.na(junctions_df$pVal)) && !all(is.na(junctions_df$adj.pVal))) {
  # saving the junctions dataframe 
  write.table(junctions_df, junctions_df_path, sep="\t", row.names = FALSE)
  cat(1)
} else {
  cat(0)
}









