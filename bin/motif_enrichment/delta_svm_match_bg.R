#!/usr/bin/env Rscript
library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(IRanges)

# Input:
# a bed file: motif_bed
# motif_id

# Output:
# bed file
# positive fasta file
# negative fasta file

# Workflow:
# Read both, filter, run function.


args = commandArgs(trailingOnly=TRUE)

# test whether there are enough arguments provided: if not, return error message
if (length(args) < 2) {
        stop("At least something input should be provided", call.=FALSE)
}


# BSgenome.Hsapiens.UCSC.hg38.masked
genNullSeqs(
  inputBedFN = args[1],
  genome = BSgenome.Hsapiens.UCSC.hg38.masked,
  outputBedFN = args[2],
  outputPosFastaFN = paste0('tmp.posSet.dhs.fa'),
  outputNegFastaFN = paste0('tmp.negSet.dhs.fa'),
    xfold = 1,
    repeat_match_tol = 1,
    GC_match_tol = 0.02,
    length_match_tol = 0.02,
    batchsize = 10000, # large: no problem with memory
    nMaxTrials = 20 # large: more sequences are matched. Vice versa
)
