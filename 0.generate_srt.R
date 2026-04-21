# ==============================================================================
# SETUP AND DIRECTORY MANAGEMENT
# ==============================================================================

# 1. Parse command line arguments
# This allows the script to be run as a job array on HPC (e.g., Slurm), 
# where 'arg' corresponds to the index of the sample to process.
args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

# 2. Load necessary libraries for spatial analysis, plotting, and data handling
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)      # For reading .h5 files
library(data.table) # For high-performance data manipulation
library(sf)         # For handling spatial vector data (GeoJSON)
library(qs)         # For fast reading/writing of R objects


# 3. Load custom internal functions
# 'functions.R' likely contains helper utils and 'do.spanorm'.
# 'banksy.R' contains the wrapper for the BANKSY spatial clustering algorithm.
source("~/VisHD/functions.R")


# 4. Define file paths
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/raw/LUT*", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
sample <- sapply(strsplit(path, split = "/"), "[", 6)

cat("===========working in", path, "=================\n")

# Create the working directory if it doesn't exist and set it as the working dir
if (!dir.exists(paste0("~/VisHD/", sample))) dir.create(paste0("~/VisHD/", sample), recursive = TRUE)
setwd(paste0("~/VisHD/", sample))

# ==============================================================================
# DATA LOADING AND PRE-PROCESSING
# ==============================================================================

# 5. Load Expression Data
localdir <- paste0(path, "/outs/segmented_outputs/")
# Read the 10X HDF5 gene expression matrix
counts <- Read10X_h5(paste0(localdir, "filtered_feature_cell_matrix.h5"), unique.features = T)

# Clean Cell IDs:
# Converts format "cellid_000042818-1" to integer "42818" to match GeoJSON IDs.
simple_ids <- as.integer(gsub("cellid_|-1", "", colnames(counts)))
colnames(counts) <- simple_ids

# 6. Load Spatial Segmentation Data (GeoJSON)
sf_use_s2(FALSE) # Turn off spherical geometry for planar spatial data
# Read cell boundaries/segmentations
geo_data <- st_read(paste0(localdir, "graphclust_annotated_cell_segmentations.geojson"))

# Create metadata frame linking cell IDs to initial graph cluster annotations
metadata_df <- data.frame(
  cell = geo_data$cell_id,         
  # Extract the cluster name from the JSON 'classification' string using regex
  SR_Cluster = sub('.*"name": "([^"]+)".*', '\\1', geo_data$classification)
)
rownames(metadata_df) <- metadata_df$cell

# Filter metadata to match only cells present in the count matrix
metadata_df <- metadata_df[colnames(counts), ]

# 7. Initialize Seurat Object
srt <- CreateSeuratObject(counts = counts, assay = "Spatial", meta.data = metadata_df)

# ==============================================================================
# SPATIAL COORDINATES SETUP
# ==============================================================================

# 8. Calculate Centroids
# Get X,Y coordinates of cell centers from the polygon geometries
centroids <- st_coordinates(st_centroid(geo_data))
coords_df <- data.frame(
  x = centroids[, 1],
  y = centroids[, 2]
)
rownames(coords_df) <- geo_data$cell_id

# 9. Add Field of View (FOV) to Seurat
# This stores spatial coordinates allowing for ImageFeaturePlots
srt[["fov"]] <- CreateFOV(
  coords = CreateCentroids(coords_df[colnames(srt), ]),
  type = "centroids",
  assay = "Spatial"
)

# ==============================================================================
# MORPHOLOGY METRICS (AREA CALCULATION)
# ==============================================================================

# 10. Define Shoelace Formula
# A mathematical algorithm to calculate the area of a polygon given its vertices.
calculate_area <- function(x, y) {
  n <- length(x)
  i <- 1:(n - 1)
  area <- 0.5 * abs(sum(x[i] * y[i + 1] - x[i + 1] * y[i]))
  return(area)
}

# 11. Calculate Cell Area
coords_dt <- as.data.table(st_coordinates(geo_data$geometry))
coords_dt$cell <- geo_data$cell_id[coords_dt$L2] # Map coordinates to cell IDs
# Apply Shoelace formula per cell
area_dt <- coords_dt[, .(area = calculate_area(X, Y)), by = cell]
cell_area <- setNames(area_dt$area, area_dt$cell)
srt$cell_area <- as.numeric(cell_area[colnames(srt)]) # Add to Seurat metadata

# 12. Calculate Nucleus Area
# Loads a separate GeoJSON for nuclei and repeats the area calculation
nucleus_geo <-  st_read(paste0(localdir, "nucleus_segmentations.geojson"))
coords_dt <- as.data.table(st_coordinates(nucleus_geo))
coords_dt$cell <- geo_data$cell_id[coords_dt$L2]
area_dt <- coords_dt[, .(area = calculate_area(X, Y)), by = cell]
nucleus_area <- setNames(area_dt$area, area_dt$cell)
srt$nucleus_area <- as.numeric(nucleus_area[colnames(srt)])

# Save raw object before filtering
qsave(srt, "raw_srt.qs")

# Build binned srt with only 16um bins==========================================
if (!dir.exists("bined_ouput")) dir.create("bined_ouput", recursive = TRUE)
srt <- Load10X_Spatial(data.dir = paste0(path, "/outs/"), bin.size = c(16))
qsave(srt, "bined_ouput/srt.qs")


# # QC AND FILTERING==============================================================
# 
# 
# # 13. Filter Low Quality Cells
# # Special handling for sample "LUT-245-20" (fixed cutoff > 25 genes)
# # All other samples use a dynamic cutoff (bottom 5% of cells removed)
# if(grepl("LUT-245-20", path)){
#   srt <- subset(srt, cells = colnames(srt)[srt$nFeature_Spatial > 25]) 
# } else {
#   srt <- subset(srt, cells = colnames(srt)[srt$nFeature_Spatial > quantile(srt$nFeature_Spatial, 0.05)])
# }
# 
# 
# # MODULE SCORING (TUMOUR VS NORMAL) ============================================
# 
# 
# # 14. Load Gene Signatures
# clean_module <- readRDS("~/CASCADEpaper/paper/Fig5_archetype/clean_module.Rds")
# names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis") # Renaming tumour modules
# source("~/CASCADEpaper/paper/normal_cells_202406/normal_markers.R") # Loads 'all_marker' (Normal cells)
# 
# if (!dir.exists("png")) dir.create("png", recursive = TRUE)
# 
# # 15. Normalize Data
# # Running custom normalization function (likely spatially aware or SCTransform-like)
# srt <- do.spanorm(srt)
# qsave(srt, "spanorm_srt.qs")
# 
# # 16. Generate Overview Plots
# spatial_plot(srt, outdir = "png/", name = "spanorm")
# 
# # 17. Integrate External Classification (CB)
# # Loads a CSV containing binary classifications (likely "Cell Binary")
# CB <- read.csv(paste0(path, "CB.csv"))
# CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode)) # Fix IDs to match integer format
# srt$category <- NA
# srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
# srt$category[srt$category == "CB NA"] <- "DT" # Handle NAs, labeling them "DT" (likely Doublet/Trash)
# 
# # 18. Calculate Aggregate Tumour vs Normal Scores
# tumour_modules <- intersect(unlist(clean_module), rownames(srt))
# normal_modules <- intersect(unlist(all_marker), rownames(srt))
# 
# srt <- AddModuleScore(srt, list(tumour_modules))
# srt$tumour_score <- srt$Cluster1 # Score for all tumour genes combined
# 
# srt <- AddModuleScore(srt, list(normal_modules))
# srt$normal_score <- srt$Cluster1 # Score for all normal genes combined
# 
# # Plot Tumour vs Normal scores
# ImageFeaturePlot(srt, c("tumour_score", "normal_score"), cols = c("white", "red"))
# ggsave("png/spanorm_ImageFeaturePlot_TN_score.png", width = 10, height = 5, dpi = 350)
# FeaturePlot(srt, c("tumour_score", "normal_score"), cols = c("white", "red"))
# ggsave("png/spanorm_FeaturePlot_TN_score.png", width = 10, height = 5, dpi = 350)
# 
# # ==============================================================================
# # DETAILED MODULE SCORING
# # ==============================================================================
# 
# # 19. Score Specific Normal Cell Types
# srt <- AddModuleScore(srt, all_marker)
# # Rename the resulting metadata columns (Cluster1, Cluster2...) to actual cell types
# colnames(srt@meta.data)[colnames(srt@meta.data) %in% paste0("Cluster", 1:13)] <- names(all_marker)
# 
# g <- ImageFeaturePlot(srt, names(all_marker), cols = c("white", "red"))
# ggsave("png/spanorm_ImageFeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)
# 
# g <- FeaturePlot(srt, names(all_marker), cols = c("white", "red"))
# ggsave("png/spanorm_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)
# 
# # 20. Score Specific Tumour Archetypes (AR, NE, etc.)
# srt <- AddModuleScore(srt, clean_module)
# colnames(srt@meta.data)[colnames(srt@meta.data) %in% paste0("Cluster", 1:6)] <- names(clean_module)
# 
# g <- ImageFeaturePlot(srt, names(clean_module), cols = c("white", "red"))
# ggsave("png/spanorm_ImageFeaturePlot_tumour_score.png", plot = g, width = 25, height = 15, dpi = 350)
# 
# g <- FeaturePlot(srt, names(clean_module), cols = c("white", "red"))
# ggsave("png/spanorm_FeaturePlot_tumour_score.png", plot = g, width = 25, height = 15, dpi = 350)
# 
# # ==============================================================================
# # DIFFERENTIAL EXPRESSION ANALYSIS (DEG)
# # ==============================================================================
# 
# # 21. Find Markers
# deg <- FindAllMarkers(srt, test.used = "MAST") # Using MAST test for single-cell DE
# deg <- deg %>% filter(p_val_adj < 0.05)
# 
# # Extract top 5 genes per cluster by log fold change
# top5_genes <- deg %>%
#   filter(p_val_adj < 0.05) %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 5)
# 
# saveRDS(deg, "deg.Rds")
# 
# # Plot top markers
# f3 <- ImageFeaturePlot(srt, top5_genes$gene, size = 0.3)
# ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)
# 
# # ==============================================================================
# # SPATIAL CLUSTERING (BANKSY) seperate to second script
# # ==============================================================================
# 
# # 22. Run BANKSY Algorithm
# # Banksy uses spatial neighborhood information to improve clustering
# # srt <- RunBanksy(srt, lambda = 0.8, verbose = TRUE, use_agf = TRUE,
# #                  assay = 'SCT', slot = 'data', features = 'all',
# #                  k_geom = 15, assay_name = 'BANKSY_0.8')
# # 
# # # 23. Standard Seurat Clustering on BANKSY Results
# # srt <- FindVariableFeatures(srt) %>% 
# #   ScaleData(srt) %>%
# #   RunPCA(npcs = 30) %>% 
# #   RunUMAP(dims = 1:20) %>% 
# #   FindNeighbors(reduction = "pca", dims = 1:20) %>% 
# #   FindClusters(resolution = 1, algorithm = 4) # Algorithm 4 = Leiden clustering
# # 
# # # 24. Final Output
# # spatial_plot(srt, "png/", "Bansky_lam0.8")
# # qsave(srt, "banksy_srt.qs")