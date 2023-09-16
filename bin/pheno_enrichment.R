#!/usr/bin/env Rscript
library(readr)
library(tidyverse)
library(stats)
library(PRROC)
library(reticulate)

np <- import("numpy")

# Input parameteres:
# 1.) nmf matrix file
# 2.) masterlist file
# 3.) Output file
# Rscript pheno_enrichment.R arg1 arg2 arg3 output.csv

args = commandArgs(trailingOnly=TRUE)

# test whether there are enough arguments provided: if not, return error message
if (length(args) < 3) {
        stop("At least something input should be provided", call.=FALSE)
}

print("Reading input matrix")
DHS_feature_nmf <- np$load(args[1])
DHS_feature_nmf_df <- data.frame(t(DHS_feature_nmf))

# This going to take some time, suggestion: input nmf matrix + chunk_id
print("Reading masterlist file")
masterlist_df <- read.csv(args[2], sep='\t', header = TRUE, check.names = FALSE)

# Combined chunk_id with nmf, preserved the order
matrix_chunk_id <- cbind(chunk_id = masterlist_df$chunk_id, DHS_feature_nmf_df)

print("reading phenotype indicator file")
indicator_bed_file <- read.table(args[3])
colnames(indicator_bed_file) <- c("chr", "start", "end", "neglog10_P", "chunk_id", "indicator")

print("merging and append matrix")
merged_df <- merge(indicator_bed_file[c("chunk_id", "indicator")], matrix_chunk_id, by = "chunk_id", sort=FALSE)

print("Checking order preservation")
# checking whether in the same order, should abort if similar
same_order <- identical(merged_df$chunk_id, indicator_bed_file$chunk_id[match(merged_df$chunk_id, indicator_bed_file$chunk_id)])
if (!same_order) {
    print('error in order')
}

# Split dataset to training set and test set
#print("Split dataset to train and test")
#training_set <- subset(combine_df, chr != 'chr7', select = -chr)
#test_set_component <- subset(combine_df, chr == 'chr7', select = -c(chr, indicator))
#test_set_indicator <- subset(combine_df, chr == 'chr7', select = indicator)

print("Running logistic regression model on training set")
log_model = glm(indicator ~., data=subset(merged_df, select = -chunk_id), family=binomial(link="logit"))

coef_df <- as.data.frame(summary(log_model)['coefficients'])
coef_df['motif_id'] <- motif_id_args
names(coef_df)[1] <- "estimate"
names(coef_df)[2] <- "std_error"
names(coef_df)[3] <- "z_value"
names(coef_df)[4] <- "Pr(>|z|)"

# Save in .csv format
print("Writing roc_pr csv file")
csv_name1 <- paste(motif_id_args, '.', n_components, ".metrics.tsv", sep='')
write.table(roc_pr_df, file=csv_name1, quote=FALSE, sep='\t')

# Save in .csv format X1, X2, ..., X16, PRROC score, PRAUC score
csv_name2 <- paste(motif_id_args, '.', n_components, ".coeff.tsv", sep='')
print("Writing coefficient csv file")
write.table(coef_df, file=csv_name2, quote=FALSE, sep='\t')