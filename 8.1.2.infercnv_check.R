library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")

in_dir    <- path.expand("~/VisHD/8.1.normal_cell_integration")
full_path <- file.path(in_dir, "integrated_pearson_srt2.qs2")
out_dir   <- file.path(in_dir, "infercnv_check")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading srt:", full_path, "\n")
srt <- qs_read(full_path)

# Pick whichever celltype_annotation column is present, preferring no-batch.
ct_candidates <- c("celltype_annotation",
                   "celltype_annotation_nobatch",
                   "celltype_annotation_batch")
ct_col <- ct_candidates[ct_candidates %in% colnames(srt@meta.data)][1]
if (is.na(ct_col)) stop("No celltype_annotation* column found in srt@meta.data")
cat("Using celltype_annotation column:", ct_col, "\n")

# Test (observation) cells: tumour_score > 0
test_cells <- colnames(srt)[!is.na(srt$tumour_score) & srt$tumour_score > 0]
# Reference pool: tumour_score < 0 AND celltype is not SVEC
ref_pool   <- colnames(srt)[
  !is.na(srt$tumour_score) & srt$tumour_score < 0 &
  srt@meta.data[[ct_col]] != "SVEC"
]

set.seed(42)
n_ref <- 5000
if (length(ref_pool) > n_ref) {
  ref_cells <- sample(ref_pool, n_ref)
} else {
  ref_cells <- ref_pool
  cat("[warn] ref_pool has only", length(ref_pool),
      "cells (< requested", n_ref, ")\n")
}

keep_cells <- c(test_cells, ref_cells)
srt <- subset(srt, cells = keep_cells)
cat("Test cells:", length(test_cells),
    "| Reference cells:", length(ref_cells),
    "| Total:", ncol(srt), "\n")

annotations <- data.frame(
  cluster   = ifelse(colnames(srt) %in% ref_cells, "Reference", "Test"),
  row.names = colnames(srt)
)
cat("Annotation groups:",
    paste(names(table(annotations$cluster)), table(annotations$cluster),
          sep = "=", collapse = ", "), "\n")

counts <- as.matrix(GetAssayData(srt, assay = "Spatial", layer = "counts"))

infercnvobject <- CreateInfercnvObject(
  raw_counts_matrix       = counts,
  annotations_file        = annotations,
  delim                   = "\t",
  min_max_counts_per_cell = c(50, Inf),
  chr_exclude             = c("chrY", "chrM"),
  gene_order_file         = readRDS("~/VisHD/gene_ord2.Rds"),
  ref_group_names         = "Reference"
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
