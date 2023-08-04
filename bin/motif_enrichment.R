#!/usr/bin/env Rscript
library(readr)
library(tidyverse)
library(reticulate)
library(stats)
library(PRROC)


np <- import("numpy")

# Input parameteres:
# 1.) DHS_nmf file
# 2.) Motif indicator file
# 3.) Output prefix - prefix of output filenames (Use motif_id)
# 4.) Number of components
# 5.) Metadata file
# Rscript motif_enrichment.R arg1 arg2 arg3 output.csv

args = commandArgs(trailingOnly=TRUE)

# test whether there are enough arguments provided: if not, return error message
if (length(args) < 3) {
        stop("At least something input should be provided", call.=FALSE)
}

print("Reading input matrix")
DHS_feature_nmf <- np$load(args[1])
DHS_feature_nmf_df <- data.frame(t(DHS_feature_nmf))

# Include gc_count and gc_count^2
gc_dataset <- read.csv('/net/seq/data2/projects/afathul/motif_enhancement/regions_gc_annotated.bed.gz', sep='\t', header = FALSE)
colnames(gc_dataset) <- c('chr', 'str', 'end', 'total', 'gc_count', 'gc_content', 'total_unique', 'chunk', 'unknown')
DHS_feature_nmf_df$gc_count <- gc_dataset$gc_count
DHS_feature_nmf_df$gc_count2 <- (gc_dataset$gc_count)^2

print("Reading Motif Indicator")
motif_indicator <- read.table(args[2])
motif_count_sum <- sum(motif_indicator)
motif_indicator <- unlist(motif_indicator)

motif_id_args <- args[3]
n_components <- args[4]

print("Running logistic regression model")
log_model = glm(motif_indicator ~., data=DHS_feature_nmf_df, family=binomial(link="logit"))

# Set up dataframe for coefficents and rename column

coef_df <- as.data.frame(summary(log_model)['coefficients'])
coef_df['motif_id'] <- motif_id_args
names(coef_df)[1] <- "estimate"
names(coef_df)[2] <- "std_error"
names(coef_df)[3] <- "z_value"
names(coef_df)[4] <- "Pr(>|z|)"
coef_df['ncomponents'] <- n_components
coef_df['components'] <- row.names(coef_df)

# Calculate AUC & PR
pred_probs = predict(log_model, type = "response")
roc_score <- roc.curve(scores.class0 = pred_probs, weights.class0 = motif_indicator, curve = TRUE)
pr_score <- pr.curve(scores.class0 = pred_probs, weights.class0 = motif_indicator, curve = TRUE)

# save the list
roc_pr <- list(
        roc_auc = roc_score$auc, 
        pr_auc_integral = pr_score$auc.integral,
        pr_auc_davis = pr_score$auc.davis.goadrich,
        motif_count = motif_count_sum, 
        motif_id = motif_id_args
)

roc_pr_df <- as.data.frame(roc_pr)
roc_pr_df['ncomponents'] <- n_components

# Save in .csv format
print("Writing roc_pr csv file")
csv_name1 <- paste(motif_id_args, '.', n_components, ".metrics.tsv", sep='')
write.table(roc_pr_df, file=csv_name1, quote=FALSE, sep='\t')

# Save in .csv format X1, X2, ..., X16, PRROC score, PRAUC score
csv_name2 <- paste(motif_id_args, '.', n_components, ".coeff.tsv", sep='')
print("Writing coefficient csv file")
write.table(coef_df, file=csv_name2, quote=FALSE, sep='\t')