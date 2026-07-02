# 9.5.integrative_infercnv.r
# InferCNV on the combined tumour + normal integration from 10.0. InferCNV
# groups: normal cells use `final_annotation`; tumour cells use
# "<slide>_<subclone>". Each group is down-sampled to 1000 cells if larger.
# Normal annotation groups serve as the diploid reference. Run settings mirror
# 2.2.infercnv_test.R (cell-level Spatial assay here, not the 16um bins).
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")

in_path <- path.expand("~/VisHD/10.0.tumour_normal_integration/integrated_pearson_srt.qs2")
out_dir <- path.expand("~/VisHD/10.2.integrative_infercnv/infercnv_clustered_all")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading 10.0 integrated srt:", in_path, "\n")
srt <- qs_read(in_path)
srt <- JoinLayers(srt, assay = "Spatial")
cat("  ", ncol(srt), "cells x", nrow(srt), "genes\n")

# ── InferCNV group: normal -> final_annotation; tumour -> "<slide>_<subclone>" ─
is_tumour <- as.character(srt$final_annotation) == "Tumour"
subclone  <- as.character(srt$subclone)
subclone[is.na(subclone)] <- "unknown"
srt$infercnv_group <- ifelse(is_tumour,
                             paste0(srt$slide, "_", subclone),
                             as.character(srt$final_annotation))
cat("\nInferCNV groups (before down-sampling):\n")
print(table(srt$infercnv_group))

# ── Down-sample each group to 1000 cells if larger ────────────────────────────
set.seed(42)
keep_cells <- unlist(lapply(split(colnames(srt), srt$infercnv_group), function(cells) {
  if (length(cells) > 1000) sample(cells, 1000) else cells
}), use.names = FALSE)
srt <- subset(srt, cells = keep_cells)
cat("\nAfter down-sampling to <=1000/group:", ncol(srt), "cells\n")
print(table(srt$infercnv_group))

# Normal annotation groups are the diploid reference.
ref_clusters <- sort(unique(srt$infercnv_group[as.character(srt$final_annotation) != "Tumour"]))
cat("\nReference (normal) groups:\n"); print(ref_clusters)

annotations <- data.frame(
  cluster   = as.character(srt$infercnv_group),
  row.names = colnames(srt)
)

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
  analysis_mode      = "samples",
  cluster_by_groups  = TRUE,
  denoise            = FALSE,
  HMM                = FALSE,
  save_rds           = TRUE,
  plot_steps         = FALSE,
  write_phylo        = TRUE,
  write_expr_matrix  = FALSE,
  no_plot            = FALSE,
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
