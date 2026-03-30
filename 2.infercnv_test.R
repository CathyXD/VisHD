# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(fastCNV)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs, lib.loc = "~/R_Library/4.5")  
library(UCell, lib.loc = "~/R_Library/4.5")  
library(pals)
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")
library(qs2)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
cat("working on", path , "\n")
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]


setwd(paths[[i]])
setwd("bined_ouput")
srt <- qread("srt.qs")
srt@meta.data <- srt@meta.data[, 1:5]
srt[["genomicScores"]] <- NULL
normal_modules <- intersect(unlist(all_marker), rownames(srt))

srt <- AddModuleScore_UCell(srt, features = list(tumour_score = c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8"), normal_score = normal_modules), storeRanks = T)

srt$tumour_cutoff <- binarise_expression(srt$tumour_score_UCell,min_expr = 0,  plot_out = "tumour_score_binarisation.png")
srt$normal_cutoff <- binarise_expression(srt$normal_score_UCell,min_expr = 0,  plot_out = "normal_score_binarisation.png")
srt$tumour_normal <- paste(srt$tumour_cutoff, srt$normal_cutoff)
g <- ImageFeaturePlot(srt, c("tumour_score_UCell", "normal_score_UCell"))+ImageDimPlot(srt, group.by = "tumour_normal")
ggsave(plot = g, "tumour_normal_score.png", width = 9,height = 3)
qs_save(srt, "srt_infercnv.qs2")
# library(fastCNVdata,  lib.loc = "~/R_Library/4.5")
# HDBreast <- load_HDBreast()
# HDBreast <- annotations8umTo16um(HDBreast, referenceVar = "annots_8um")
# ImageDimPlot(HDBreast, group.by = "projected_annots_8um")
# HDBreast$tumour_normal <- HDBreast$projected_annots_8um
# DefaultAssay(HDBreast) <- "Spatial.016um"
# HDBreast[["Spatial.008um"]] <- NULL
# 
# normal_HDBreast <- subset(HDBreast, cells = colnames(HDBreast)[HDBreast$tumour_normal == "NoTumor" & HDBreast$nFeature_Spatial.016um>100])
#normal_HDBreast <- qs_read("~/VisHD/HDBreast_normal.qs2")
#merged_hd <- merge(
#  x = srt,
 # y = normal_HDBreast,
#  add.cell.ids = c("LUT", "Breast")
#)

#merged_hd <- JoinLayers(merged_hdk)

# srt$cnv_cutoff <- binarise_expression(srt$cnv_fraction, plot_out = "cnv_fraction_binarisation.png")
# test <- as.data.frame.matrix(table(srt$cnv_clusters, srt$cnv_cutoff))
# test <- apply(test, 1, function(x) x/sum(x))
# ref <- names(which.max(test["0", ] ))
# 

library(infercnv)
infercnvobject = CreateInfercnvObject(raw_counts_matrix=GetAssayData(srt, assay = "Spatial.016um", layer = "counts"),
                                      annotations_file=data.frame(srt$tumour_normal),
                                      delim="\t",
                                      gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                      ref_group_names="0 1")

if (!file.exists("infercnv_test")){
  dir.create("infercnv_test")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_test",
                               analysis_mode = "cells",
                               cluster_by_groups= T,
                               denoise=F,
                               HMM=F,
                               save_rds = F,
                               write_phylo = F, 
                               write_expr_matrix = F,
                               output_format = "pdf",
                               num_threads = 6)
