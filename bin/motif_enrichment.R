#!/usr/bin/env Rscript
library(readr)
library(tidyverse)
library(stats)
library(PRROC)
library(reticulate)


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
print("Combined dataset")
gc_dataset <- read.csv('/net/seq/data2/projects/afathul/motif_enhancement/regions_gc_annotated.bed.gz', sep='\t', header = FALSE)
colnames(gc_dataset) <- c('chr', 'str', 'end', 'total', 'gc_count', 'gc_content', 'total_unique', 'chunk', 'unknown')

print("Reading Motif Indicator")
motif_indicator <- read.table(args[2], col.names='indicator')
motif_count_sum <- sum(motif_indicator)
# motif_indicator <- unlist(motif_indicator)

motif_id_args <- args[3]
n_components <- args[4]

# Combining as one dataframe: Indicator, Matrix, GC
combine_df <- cbind(motif_indicator, DHS_feature_nmf_df,gc_dataset['gc_count'], gc_count2 = (gc_dataset$gc_count)^2, gc_dataset['chr'])
# gc_dataset['gc_count'], gc_count2 = (gc_dataset$gc_count)^2,

# Split dataset to training set and test set
print("Split dataset to train and test")
training_set <- subset(combine_df, chr != 'chr7', select = -chr)
test_set_component <- subset(combine_df, chr == 'chr7', select = -c(chr, indicator))
test_set_indicator <- subset(combine_df, chr == 'chr7', select = indicator)

print("Running logistic regression model on training set")
log_model = glm(indicator ~., data=training_set, family=binomial(link="logit"))

# Set up dataframe for coefficents and rename column
coef_df <- as.data.frame(summary(log_model)['coefficients'])
coef_df['motif_id'] <- motif_id_args
names(coef_df)[1] <- "estimate"
names(coef_df)[2] <- "std_error"
names(coef_df)[3] <- "z_value"
names(coef_df)[4] <- "Pr(>|z|)"
coef_df['ncomponents'] <- n_components
coef_df['components'] <- row.names(coef_df)

# Using AUPRG score from PRG packages in python
print("Running prg package")
library(prg)

print("Run actual and predicted_prob")
actual <- test_set_indicator$indicator
predicted_prob <- predict(log_model, test_set_component, type = "response")

print("prg_curve with create_prg_curve")
prg_curve <- create_prg_curve(actual, predicted_prob)
print("calculate prg score")
auprg_val <- calc_auprg(prg_curve)


# Calculate AUC & PR
print("Predict Train and Set model")
pred_probs_train = predict(log_model, type = "response")
pred_probs_test = predict(log_model, test_set_component, type = "response")

roc_score_train <- roc.curve(scores.class0 = pred_probs_train,
                                 weights.class0 = training_set$indicator, curve = TRUE)
pr_score_train <- pr.curve(scores.class0 = pred_probs_train,
                                 weights.class0 = training_set$indicator, curve = TRUE)

roc_score_test <- roc.curve(scores.class0 = pred_probs_test,
                                 weights.class0 = unlist(test_set_indicator), curve = TRUE)
pr_score_test <- pr.curve(scores.class0 = pred_probs_test,
                                 weights.class0 = unlist(test_set_indicator), curve = TRUE)

# save the list
roc_pr <- list(
        roc_auc_train = roc_score_train$auc, 
        pr_auc_integral_train = pr_score_train$auc.integral,
        pr_auc_davis_train = pr_score_train$auc.davis.goadrich,
        roc_auc_test = roc_score_test$auc, 
        pr_auc_integral_test = pr_score_test$auc.integral,
        pr_auc_davis_test = pr_score_test$auc.davis.goadrich,
        motif_count = motif_count_sum, 
        motif_id = motif_id_args
        auprg_score = auprg_val
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