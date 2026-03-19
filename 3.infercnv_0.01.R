args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)
setwd("bined_ouput")
cat("working on", path , "bined_ouput\n")

library(infercnv)
library(qs,lib.loc = "~/R_Library/4.5")
srt <- qread("srt_fastCNV.qs")

refs <- list("1" = "4", 
            "2" = "2",
            "3" = "1",
            "4" = "3")

ref <- refs[arg]

infercnvobject = CreateInfercnvObject(raw_counts_matrix=GetAssayData(srt, assay = "Spatial.016um", layer = "counts"),
                                        annotations_file=data.frame(srt$cnv_clusters),
                                        delim="\t",
                                        gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                        ref_group_names=ref)
  if (dir.exists("infercnv_0.01")) {
    unlink("infercnv_0.01", recursive = TRUE)
    message("Folder deleted.")
  } else {
    message("Folder not found.")
  }

  if (!file.exists("infercnv_0.01")){
    dir.create("infercnv_0.01")
  }

  infercnvobject = infercnv::run(infercnvobject,
                                 cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir="infercnv_0.01",
                                 analysis_mode = "cells",
                                 cluster_by_groups= T,
                                 denoise=F,
                                 HMM=F,
                                 save_rds = F,
                                 write_phylo = F, 
                                 write_expr_matrix = F,
                                 output_format = "pdf",
                                 num_threads = 6)

  
if (dir.exists("infercnv_uncluster")) {
  unlink("infercnv_uncluster", recursive = TRUE)
  message("Folder deleted.")
 } else {
  message("Folder not found.")
 }
  
 if (!file.exists("infercnv_uncluster")){
   dir.create("infercnv_uncluster")
 }
  
  infercnvobject = infercnv::run(infercnvobject,
                                 cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir="infercnv_uncluster",
                                 analysis_mode = "cells",
                                 cluster_by_groups= F,
                                 denoise=F,
                                 HMM=F, 
                                 save_rds = F,
                                 write_phylo = T, 
                                 write_expr_matrix = T, 
                                 output_format = "pdf",
                                 num_threads = 6)
  
  
  cat("Done")