library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")

in_dir  <- path.expand("~/VisHD/1.1.integrate_raw_cell")
full_path <- file.path(in_dir, "integrated_pearson_srt.qs2")
out_dir <- file.path(in_dir, "infercnv_clustered_sample")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading full integrated srt:", full_path, "\n")
srt <- qs_read(full_path )

# Drop SVEC cluster (cluster 6) before any downstream work — those cells
# aren't tumour and aren't a useful diploid reference for InferCNV.
svec_cluster <- "6"
n_before <- ncol(srt)
srt <- subset(srt, cells = colnames(srt)[as.character(srt$pearson_clusters) != svec_cluster])
cat("Dropped SVEC cluster", svec_cluster, ":", n_before, "->", ncol(srt), "cells\n")

# Stratified subsample by pearson_clusters: InferCNV at full ~320k cells
# OOM'd at 200 GiB (43.2 GiB sparse→dense coercion alone). 50 k keeps cluster
# proportions and reference-cluster counts intact.
set.seed(42)
target_n <- 50000
if (ncol(srt) > target_n) {
  cells_by_cluster <- split(colnames(srt), srt$pearson_clusters)
  cluster_sizes    <- lengths(cells_by_cluster)
  alloc <- round(cluster_sizes / sum(cluster_sizes) * target_n)
  alloc <- pmin(alloc, cluster_sizes)
  small <- alloc == 0 & cluster_sizes > 0
  alloc[small] <- pmin(50, cluster_sizes[small])
  keep_cells <- unlist(Map(sample, cells_by_cluster, alloc))
  srt <- subset(srt, cells = keep_cells)
  cat("Subsampled from", sum(cluster_sizes), "to", ncol(srt),
      "cells, stratified by pearson_clusters\n")
}

ref_clusters <- as.character(c(0, 1, 12, 19, 4, 3, 5, 20, 9, 17))

# Annotation: reference clusters keep their cluster ID (so InferCNV uses them
# as the diploid baseline). Every other cell (tumour) is labelled by its
# slide so subclonal CNV structure is grouped per sample.
cluster_chr <- as.character(srt$pearson_clusters)
slide_chr   <- as.character(srt$slide)
group       <- ifelse(cluster_chr %in% ref_clusters, cluster_chr, slide_chr)
annotations <- data.frame(
  cluster   = group,
  row.names = colnames(srt)
)
cat("Annotation groups (count): ",
    paste(names(table(group)), table(group), sep = "=", collapse = ", "), "\n")

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
