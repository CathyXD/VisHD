# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")

in_dir  <- path.expand("~/VisHD/1.integrate_raw_binned")
sub_path <- file.path(in_dir, "integrated_pearson_srt_sub30k.qs2")
full_path <- file.path(in_dir, "integrated_pearson_srt.qs2")
out_dir <- file.path(in_dir, "infercnv_clustered_all")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading full integrated srt:", full_path, "\n")
srt <- qs_read(full_path )

ref_clusters <- c("1", "2", "0", "7", "4", "5", "19", "21", "14", "11")

annotations <- data.frame(
  cluster = as.character(srt$pearson_clusters),
  row.names = colnames(srt)
)

counts <- as.matrix(GetAssayData(srt, assay = "Spatial.016um", layer = "counts"))

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
  analysis_mode      = "samples",
  cluster_by_groups  = TRUE,
  denoise            = FALSE,
  HMM                = FALSE,
  save_rds           = TRUE,
  plot_steps         = FALSE,
  write_phylo        = TRUE,
  write_expr_matrix  = FALSE,
  no_plot            = TRUE,
  num_threads        = 8,
  resume_mode        = TRUE
)

plot_per_group(
  infercnvobject,
  on_references     = TRUE,
  on_observations   = TRUE,
  base_filename     = "infercnv_per_group",
  output_format     = "pdf",
  write_expr_matrix = FALSE,
  png_res           = 300,
  dynamic_resize    = 0,
  useRaster         = TRUE,
  out_dir           = out_dir
)
