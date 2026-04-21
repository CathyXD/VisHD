# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(fastCNV)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs, lib.loc = "~/R_Library/4.5")  
library(pals)
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")

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
  srt_cell <- qread("spanorm_srt.qs")
  srt_cell <- AddModuleScore(srt_cell, list(c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")))
  srt_cell$tumour_score <- srt_cell$Cluster1
  
  test <- ImageFeaturePlot(srt_cell, "tumour_score") + FeaturePlot(srt_cell, "tumour_score") + VlnPlot(srt_cell , "tumour_score", group.by = "SR_Cluster")
  ggsave(plot = test, "png/spanorm_tumour_score.png", width = 10, height = 4)
  srt_cell$tumour_cutoff <- binarise_expression(srt_cell$tumour_score,min_expr = min(srt_cell$tumour_score),  plot_out = "png/tumour_score_binarisation.png")
  normal_modules <- intersect(unlist(all_marker), rownames(srt_cell))
  srt_cell <- AddModuleScore(srt_cell, list(normal_modules))
  srt_cell$normal_score <- srt_cell$Cluster1
  ImageFeaturePlot(srt_cell, "normal_score") + FeaturePlot(srt_cell, "normal_score")
  srt_cell$normal_cutoff <- binarise_expression(srt_cell$normal_score,min_expr = 0,  plot_out = "png/normal_score_binarisation.png")
  srt_cell$tumour_normal <- paste(srt_cell$tumour_cutoff, srt_cell$normal_cutoff)
  ImageDimPlot(srt_cell, group.by = "tumour_normal")
  
  setwd("bined_ouput")
  srt <- qread("srt.qs")
  srt@meta.data <- srt@meta.data[, 1:5]
  srt[["genomicScores"]] <- NULL
  srt <- subset(srt, cells = colnames(srt)[srt$nFeature_Spatial.016um>quantile(srt$nFeature_Spatial.016um, 0.05)] )
  srt <- transfer_visiumhd_to_cells(
    srt,
    srt_cell,
    annotation_cols = c("SR_Cluster","tumour_cutoff" )
  )
  srt$ref <- paste(srt$tumour_cutoff, srt$K5_cluster)
  
  test <- as.matrix(table(srt$K5_cluster, srt$tumour_cutoff))
  test <- apply(test, 1 , function(x) x/sum(x))
  cat(test , "\n")
  ref <- names(which.max(test["0", ]))
  srt <- fastCNV_10XHD(srt, sampleName = i, referenceVar = "K5_cluster", referenceLabel =ref, printPlot = TRUE,
                       getCNVClusters = F,
                       saveGenomicWindows = T)
  srt <- CNVCluster(srt, k = 4)
  srt <- mergeCNVClusters(srt, mergeThreshold = 0.9)
  p <- ImageDimPlot(srt, group.by = "cnv_clusters", cols = "polychrome")
  f <- ImageFeaturePlot(srt, "cnv_fraction")
  ggsave("CNV_cluster_fraction.png", plot = p+f, width = 8, height = 4)
  plotCNVResultsHD(srt, referenceVar = "tumour_cutoff",raster_resize_mat =F)
  qsave(srt, "srt_fastCNV.qs")
  
  # srt$cnv_cutoff <- binarise_expression(srt$cnv_fraction, plot_out = "cnv_fraction_binarisation.png")
  # test <- as.data.frame.matrix(table(srt$cnv_clusters, srt$cnv_cutoff))
  # test <- apply(test, 1, function(x) x/sum(x))
  # ref <- names(which.max(test["0", ] ))
  # 

  library(infercnv)
  infercnvobject = CreateInfercnvObject(raw_counts_matrix=GetAssayData(srt, assay = "Spatial.016um", layer = "counts"),
                                        annotations_file=data.frame(srt$K5_cluster),
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

  
  if (dir.exists("infercnv_uncluster_K5anno")) {
   unlink("infercnv_uncluster_K5anno", recursive = TRUE)
   message("Folder deleted.")
  } else {
   message("Folder not found.")
  }
  
  
  if (!file.exists("infercnv_uncluster_K5anno")){
    dir.create("infercnv_uncluster_K5anno")
  }
  
  infercnvobject = infercnv::run(infercnvobject,
                                 cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir="infercnv_uncluster_K5anno",
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

  
  # if (dir.exists("infercnv_uncluster")) {
  #  unlink("infercnv_uncluster", recursive = TRUE)
  #  message("Folder deleted.")
  # } else {
  #  message("Folder not found.")
  # }
  
  # if (!file.exists("infercnv_uncluster")){
  #   dir.create("infercnv_uncluster")
  # }

# Try external reference not work
  # HDBreast <- load_HDBreast()
  # HDBreast <- annotations8umTo16um(HDBreast, referenceVar = "annots_8um")
  # ImageDimPlot(HDBreast, group.by = "projected_annots_8um")
  # HDBreast$ref <- HDBreast$projected_annots_8um
  # DefaultAssay(HDBreast) <- "Spatial.016um"
  # HDBreast <- fastCNV_10XHD(HDBreast, sampleName = "HDBreast_nocluster", referenceVar = "projected_annots_8um", referenceLabel = "NoTumor", printPlot = TRUE, getCNVClusters = F, getCNVPerChromosomeArm = F)
#   
# srt$ref <- srt$K5_cluster
# 
# seuratList <- c(srt, HDBreast)
# sampleNames <- c("LUT", "HDB")
# names(seuratList) <- sampleNames
# referencelabels <- c("1", "2", "NoTumor")
# 
# cnv <- fastCNV_10XHD(seuratList,
#                      sampleName = sampleNames, 
#                      referenceVar = "ref", 
#                      referenceLabel = referencelabels, 
#                      printPlot = TRUE,
#                      getCNVPerChromosomeArm = F,
#                      getCNVClusters = F, 
#                      pooledReference = T)


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
