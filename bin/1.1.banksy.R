library(Seurat)
library(ggplot2)
library(patchwork)
source("~/VisHD/functions.R")
library(dplyr)
library(qs)
library()
# 1. Parse command line arguments
# This allows the script to be run as a job array on HPC (e.g., Slurm), 
# where 'arg' corresponds to the index of the sample to process.
args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# 4. Define file paths
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
path <- paste0(path, "/analysis/segmented_outputs/")
setwd(path)
srt <-  qread(paste0(path, "banksy_srt.qs"))
# srt <- SeuratWrappers::RunBanksy(srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
#                                  assay = 'SpaNorm', slot = 'data',
#                                  k_geom = c(15), assay_name = 'BANKSY_0.2')
# qsave(srt, "banksy_srt.qs")
# 
srt <- RunPCA(srt, npcs = 30, features = rownames(srt), reduction.name = "banksy0.2.pca")
srt <- RunUMAP(srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name  = "banksy0.2.umap")
srt <-  FindNeighbors(srt, reduction = "banksy0.2.pca", dims = 1:20)
srt <-   FindClusters(srt,resolution = 1, algorithm = 4)
qsave(srt, "banksy_srt.qs")
# spatial_plot(srt, "png/", "Bansky_lam0.2")
# 
# 
# deg <- FindAllMarkers(srt, test.used = "MAST") # Using MAST test for single-cell DE
# deg <- deg %>% filter(p_val_adj < 0.05)
# 
# # Extract top 5 genes per cluster by log fold change
# top5_genes <- deg %>%
#   filter(p_val_adj < 0.05) %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 5)
# 
# saveRDS(deg, "banksy_lambda0.2_deg.Rds")
# 
# # Plot top markers
# f3 <- ImageFeaturePlot(srt, top5_genes$gene, size = 0.3)
# ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_BANKSY_0.2_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)
DefaultAssay(srt) <- "SpaNorm"
srt <- SeuratWrappers::RunBanksy(srt, lambda = 0.8, verbose = TRUE, use_agf = TRUE,
                                 assay = 'SpaNorm', slot = 'data',
                                 k_geom = c(15), assay_name = 'BANKSY_0.8')
qsave(srt, "banksy_srt.qs")

srt <- RunPCA(srt, npcs = 30, features = rownames(srt), reduction.name = "banksy0.8.pca") %>%
  RunUMAP(dims = 1:20, reduction = "banksy0.8.pca", reduction.name  = "banksy0.8.umap") %>%
  FindNeighbors(reduction = "banksy0.8.pca", dims = 1:20) %>%
  FindClusters(resolution = 1, algorithm = 4)
spatial_plot(srt, "png/", "Bansky_lam0.8")
qsave(srt, "banksy_srt.qs")

deg <- FindAllMarkers(srt, test.used = "MAST") # Using MAST test for single-cell DE
deg <- deg %>% filter(p_val_adj < 0.05)

# Extract top 5 genes per cluster by log fold change
top5_genes <- deg %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

saveRDS(deg, "banksy_lambda0.8_deg.Rds")

# Plot top markers
f3 <- ImageFeaturePlot(srt, top5_genes$gene, size = 0.3)
ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_BANKSY_0.8_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)


