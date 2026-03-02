
args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)
cat("working on", path , "\n")

library(Seurat)
library(BiocNeighbors)
library(Matrix)
library(qs, lib.loc = "~/R_Library/4.5")  
library(ggplot2)

# 1. Load your Segmented Visium HD Object
# Assumes you have a column in @meta.data called 'cell_type' (e.g., "Tumor", "Normal")
# and spatial coordinates in the 'spatial' or 'tissue_lowres' slot.
viz_hd <- qread("banksy_srt.qs")
DefaultAssay(viz_hd) <- "Spatial"
Tumour <- read.csv("Tumour.csv")
Tumour$cell <- as.integer(gsub("cellid_|-1", "", Tumour$Barcode))
if(!any(grepl("tumour|tumor", Tumour$Tumour))){
  Tumour$Tumour <- paste(Tumour$Tumour, "Tumour")
}

cellanno <- read.csv("K-Means 5 cellanno.csv")
cellanno$cell <-  as.integer(gsub("cellid_|-1", "", cellanno$Barcode))

viz_hd$anno <- NA
viz_hd$anno<-  cellanno[match(colnames(viz_hd), cellanno$cell), 2]

CB <- read.csv("CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
viz_hd$category <- NA
viz_hd$category <- paste("CB", CB$CB[match(viz_hd$cell, CB$cellid)])
viz_hd$category[viz_hd$category == "CB NA"] <- "DT"


DimPlot(viz_hd, group.by = "anno", cols = as.character(pals::trubetskoy(20)), split.by = "category")
ggsave("png/cell_anno.png", width = 8, height = 3)

# 2. Parameters
k_neighbors <- 20
group_col <- "anno" # Change this to your actual metadata column name

# 3. Function to perform group-restricted pseudobulking
pseudobulk_with_metadata <- function(obj, k, group_var) {
  # Extract coordinates, groups, and raw counts
  coords <- GetTissueCoordinates(obj)
  groups <- obj[[group_var, drop = TRUE]]
  counts <- GetAssayData(obj, layer = "counts", assay = "Spatial")
  
  # Initialize a sparse matrix to save memory
  # We use a list to store chunks and then combine them to avoid dense matrix overhead
  smoothed_list <- list()
  
  unique_groups <- unique(groups)
  
  for (g in unique_groups) {
    message("Processing: ", g)
    
    # Subset indices for the specific group
    group_indices <- which(groups == g)
    
    if (length(group_indices) <= k) {
      message("Group ", g, " is too small. Copying raw counts.")
      smoothed_list[[g]] <- counts[, group_indices, drop = FALSE]
      next
    }
    
    group_coords <- coords[group_indices, ]
    
    # Find k-NN within the group (k+1 includes the cell itself)
    nn_search <- findKNN(as.matrix(group_coords), k = k, warn_ties = FALSE)
    
    # Pre-allocate a matrix for this group
    group_smoothed <- matrix(0, nrow = nrow(counts), ncol = length(group_indices))
    colnames(group_smoothed) <- colnames(counts)[group_indices]
    
    for (i in seq_along(group_indices)) {
      # nn_search$index gives indices relative to the SUBSET (group_indices)
      neighbor_indices_relative <- nn_search$index[i, ]
      
      # Map those back to the original full counts matrix indices
      actual_neighbor_indices <- group_indices[neighbor_indices_relative]
      
      # Include the target cell itself
      all_indices <- c(group_indices[i], actual_neighbor_indices)
      
      # Calculate Mean
      group_smoothed[, i] <- rowSums(counts[, all_indices])
    }
    
    smoothed_list[[g]] <- as(group_smoothed, "dgCMatrix")
  }
  
  # Combine all groups back into one matrix and reorder to match original object
  combined_smoothed <- do.call(cbind, smoothed_list)
  combined_smoothed <- combined_smoothed[, colnames(obj)]
  
  # Create New Assay
  # obj[["Pseudobulk"]] <- CreateAssayObject(counts = combined_smoothed)
  
  # Ensure the metadata is explicitly linked (Seurat does this by cell name)
  # You can verify by looking at the first few rows
  message("Pseudobulk assay created. Metadata preserved for all ", ncol(obj), " cells.")
  
  return(combined_smoothed)
}

# 4. Execute
#viz_hd_smoothed <- pseudobulk_with_metadata(viz_hd, k = k_neighbors, group_var = group_col)
#colnames(viz_hd_smoothed) <- colnames(viz_hd)
#rownames(viz_hd_smoothed) <- rownames(viz_hd)
#qsave(viz_hd_smoothed, "pseudobulk_count.qs")
viz_hd_smoothed <- qread("pseudobulk_count.qs")

anno <- data.frame(viz_hd$anno)

# Run InferCNV on the pseudobulk
library(infercnv)
infercnvobject = CreateInfercnvObject(raw_counts_matrix=viz_hd_smoothed,
                                      annotations_file=anno,
                                      delim="\t",
                                      gene_order_file= readRDS("~/VisHD/gene_ord2.Rds"),
                                      ref_group_names= intersect(unique(viz_hd$anno), c("Immune","SM","Stromal", "SN","immune")))
if (!file.exists("infercnv_0.05_pseudobulk")){
  dir.create("infercnv_0.05_pseudobulk")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_0.05_pseudobulk",
                               analysis_mode = "cells",
                               cluster_by_groups= T,
                               denoise=F,
                               HMM=F, 
                               output_format = "pdf",
                               num_threads = 3)
cat("Done")


