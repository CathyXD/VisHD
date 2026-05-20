# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(fastCNV)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")  
library(UCell, lib.loc = "~/R_Library/4.5")  
library(pals)
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")


args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
cat("working on", path , "\n")
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]


setwd(paths[[i]])
setwd("bined_ouput")
srt <- qs_read("srt_infercnv.qs2")
# srt <- qread("srt.qs")
# normal_modules <- intersect(unlist(all_marker), rownames(srt))

# srt <- AddModuleScore_UCell(srt, features = list(tumour_score = c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8"), normal_score = normal_modules), storeRanks = T)

# srt$tumour_cutoff <- binarise_expression(srt$tumour_score_UCell,min_expr = 0,  plot_out = "tumour_score_binarisation.png")
# srt$normal_cutoff <- binarise_expression(srt$normal_score_UCell,min_expr = 0,  plot_out = "normal_score_binarisation.png")
# srt$tumour_normal <- paste(srt$tumour_cutoff, srt$normal_cutoff)
# g <- ImageFeaturePlot(srt, c("tumour_score_UCell", "normal_score_UCell"))+ImageDimPlot(srt, group.by = "tumour_normal")
# ggsave(plot = g, "tumour_normal_score.png", width = 9,height = 3)
# qs_save(srt, "srt_infercnv.qs2")


library(infercnv)
infercnvobject = CreateInfercnvObject(raw_counts_matrix=as.matrix(GetAssayData(srt, assay = "Spatial.016um", layer = "counts")),
                                      annotations_file=data.frame(srt$tumour_normal),
                                      delim="\t",
                                      min_max_counts_per_cell = c(50, Inf), 
                                      chr_exclude = c('chrY', 'chrM'), 
                                      gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                      ref_group_names="0 1")

if (!file.exists("infercnv_uncluster")){
  dir.create("infercnv_uncluster")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_uncluster",
                               analysis_mode='cells',
                               cluster_by_groups= F,
                               denoise=F,
                               HMM=F,
                               save_rds = T, 
                               plot_steps = F, 
                               write_phylo = T, 
                               write_expr_matrix = F,
                               output_format = "pdf",
                               num_threads = 6, 
                               resume_mode = F)
