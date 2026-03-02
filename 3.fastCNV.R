# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(fastCNV)
library(fastCNVdata, lib.loc = "~/R_Library/4.5")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs, lib.loc = "~/R_Library/4.5")  
library(pals, lib.loc = "~/R_Library/4.5")
.libPaths(c("~/R_Library/4.5", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/bined_ouput", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
cat("working on", path , "\n")
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]
ref_lab<- list(c("1", "2"), 
                      c("2", "3", "4", "5"), 
                      c("2", "3", "4"), 
                      c("1", "3","4"))
names(ref_lab) <- samples


  setwd(paths[[i]])
  srt <- qread("srt.qs")
  srt@meta.data <- srt@meta.data[, 1:5]
  srt[["genomicScores"]] <- NULL
  srt <- fastCNV_10XHD(srt, sampleName = i, referenceVar = "K5_cluster", referenceLabel =ref_lab[[i]], printPlot = TRUE, 
                       getCNVClusters = F, 
                       saveGenomicWindows = T)
  srt <- CNVCluster(srt)
  srt <- mergeCNVClusters(srt, mergeThreshold = 0.85)
  p <- ImageDimPlot(srt, group.by = "cnv_clusters", cols = "polychrome")
  f <- ImageFeaturePlot(srt, "cnv_fraction")
  ggsave("CNV_cluster_fraction.png", plot = p+f, width = 8, height = 4)
  plotCNVResultsHD(srt, referenceVar = "K5_cluster",raster_resize_mat =F)
  qsave(srt, "srt_fastCNV.qs")
  library(infercnv)
  infercnvobject = CreateInfercnvObject(raw_counts_matrix=GetAssayData(srt, assay = "Spatial.016um", layer = "counts"),
                                        annotations_file=data.frame(srt$K5_cluster),
                                        delim="\t",
                                        gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                        ref_group_names=ref_lab[[i]])
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
                                 output_format = "pdf",
                                 num_threads = 2)

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
                                 cluster_by_groups= T,
                                 denoise=T,
                                 HMM=F, 
                                 output_format = "pdf",
                                 num_threads = 4)
  
  
  cat("Done")



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



