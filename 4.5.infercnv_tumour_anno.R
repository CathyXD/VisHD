library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")

in_dir    <- path.expand("~/VisHD/4.4.integrate_tumour_anno")
full_path <- file.path(in_dir, "integrated_pearson_srt.qs2")
out_dir   <- file.path(in_dir, "infercnv_clustered_tumour_anno")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading full integrated srt:", full_path, "\n")
srt <- qs_read(full_path)

# InferCNV group label = <tumour_anno> <pearson_cluster> so reference groups
# can be selected by their Normal-prefixed cluster IDs.
srt@meta.data$infercnv_group <- paste(srt$tumour_anno, srt$pearson_clusters)

# Stratified subsample to 60 k cells by the InferCNV group so each reference
# group keeps adequate counts (rare groups get a 50-cell floor).
set.seed(42)
target_n <- 60000
if (ncol(srt) > target_n) {
  cells_by_group <- split(colnames(srt), srt@meta.data$infercnv_group)
  group_sizes    <- lengths(cells_by_group)
  alloc <- round(group_sizes / sum(group_sizes) * target_n)
  alloc <- pmin(alloc, group_sizes)
  small <- alloc == 0 & group_sizes > 0
  alloc[small] <- pmin(50, group_sizes[small])
  keep_cells <- unlist(Map(sample, cells_by_group, alloc))
  srt <- subset(srt, cells = keep_cells)
  cat("Subsampled from", sum(group_sizes), "to", ncol(srt),
      "cells, stratified by infercnv_group\n")
}

ref_clusters <- paste("Normal", c(0, 1, 4, 5, 2, 9, 22, 14, 20, 17))

annotations <- data.frame(
  cluster   = as.character(srt@meta.data$infercnv_group),
  row.names = colnames(srt)
)
cat("Annotation groups (count): ",
    paste(names(table(annotations$cluster)), table(annotations$cluster),
          sep = "=", collapse = ", "), "\n")

missing_refs <- setdiff(ref_clusters, unique(annotations$cluster))
if (length(missing_refs) > 0) {
  cat("[warn] reference groups absent after subsample:",
      paste(missing_refs, collapse = ", "), "\n")
  ref_clusters <- intersect(ref_clusters, unique(annotations$cluster))
}

counts <- as.matrix(GetAssayData(srt, assay = "Spatial", layer = "counts"))

infercnvobject <- CreateInfercnvObject(
  raw_counts_matrix       = counts,
  annotations_file        = annotations,
  delim                   = "\t",
  min_max_counts_per_cell = c(50, Inf),
  chr_exclude             = c("chrY", "chrM"),
  gene_order_file         = readRDS("~/VisHD/gene_ord2.Rds"),
  ref_group_names         = ref_clusters
)

infercnvobject <- infercnv::run(
  infercnvobject,
  cutoff             = 0.01,
  out_dir            = out_dir,
  analysis_mode      = "cells",
  cluster_by_groups  = TRUE,
  denoise            = FALSE,
  HMM                = FALSE,
  save_rds           = TRUE,
  plot_steps         = FALSE,
  write_phylo        = TRUE,
  write_expr_matrix  = FALSE,
  num_threads        = 8,
  resume_mode        = TRUE
)

plot_per_group(
  infercnvobject,
  on_references     = TRUE,
  on_observations   = TRUE,
  base_filename     = "infercnv_per_group",
  output_format     = "png",
  write_expr_matrix = FALSE,
  png_res           = 300,
  dynamic_resize    = 0,
  useRaster         = TRUE,
  out_dir           = out_dir
)
