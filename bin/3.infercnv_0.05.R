args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

cat("working on", path , "\n")

library(infercnv)
library(qs,lib.loc = "~/R_Library/4.5")

infercnvobject = CreateInfercnvObject(raw_counts_matrix=qread("countmat.qs"),
                                      annotations_file="anno.txt",
                                      delim="\t",
                                      gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                     ref_group_names=readRDS("normal_clust.Rds"))
if (!file.exists("infercnv_0.05")){
  dir.create("infercnv_0.05")
}

infercnvobject = infercnv::run(infercnvobject,
                                cutoff=0.05, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="infercnv_0.05",
                                analysis_mode = "cells",
                                cluster_by_groups= T,
                                denoise=F,
                                HMM=F, 
                                output_format = "pdf",
                                num_threads = 2)
cat("Done")