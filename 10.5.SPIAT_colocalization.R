#!/usr/bin/env Rscript
# 10.5.SPIAT_colocalization.R   (per-sample, array 1-8)
# SPIAT cell-colocalisation workflow on the per-sample tumour+normal object
# (LUT-245-XX/tumour_normal_anno_srt.qs2): quantifies how much each normal
# `cell_type` colocalizes with Tumour cells (Tumour subclones 1/2 pooled into
# a single "Tumour" reference level, same as 10.3), using three independent
# metrics per non-Tumour cell type:
#   - CIN  (Cells In Neighbourhood): % of target cells within radius r of a
#           Tumour cell, at r = 50/100/200 px (~15/30/60 µm; same radii
#           convention as 10.1.Statial.R / 10.2.tumour_distance_cell_composition.R)
#   - MS/NMS (Mixing Score / Normalised Mixing Score): ratio of
#           Tumour-target adjacency pairs to Tumour-Tumour pairs, at the
#           same 3 radii
#   - Cross-K function (Ripley's cross-K, AUC + ring-crossing fraction):
#           spatial clustering (AUC > 0) vs. dispersion (AUC < 0) between
#           Tumour and the target cell type, evaluated out to
#           dist = max(radii) = 200
#
# Reference: https://trigosteam.github.io/SPIAT/articles/cell-colocalisation.html
#
# Unlike 10.3, this script does not need `alphahull` (that's only used for
# tissue-structure/margin detection, not colocalisation).
#
# Outputs under ~/VisHD/10.5.tumour_colocalization/:
#   per_sample_tables/  <slide>_CIN.csv            cell_type, radius_px, radius_um, CIN_pct
#                        <slide>_mixing_score.csv   raw mixing_score_summary() rows + radius_px/radius_um
#                        <slide>_crossK_summary.csv cell_type, crossK_AUC, crossK_crossing_frac
#   per_sample_crossK/   <slide>_<cell_type>_crossK.csv  (r, theo, border columns from the fv object)
#                        <slide>_<cell_type>_crossK.png  (base plot(df_cross))
#
#   Rscript 10.5.SPIAT_colocalization.R <sample-index 1-8>
#
# NOTE: SPIAT must be installed (already available in the spatial:4.5v6 image;
#       see 10.3.SPIAT_tumour_border.R for the BiocManager install line if not).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs2)
  library(SpatialExperiment)
  library(S4Vectors)
  library(SPIAT)
})

# ── CLI arg ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 10.5.SPIAT_colocalization.R <sample-index 1-8>")
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[arg]
i     <- basename(path)
cat("Working at", path, "\n")

out_dir    <- "~/VisHD/10.5.tumour_colocalization"
tbl_dir    <- file.path(out_dir, "per_sample_tables")
crossk_dir <- file.path(out_dir, "per_sample_crossK")
for (d in c(tbl_dir, crossk_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

radii <- c(50, 100, 200)   # FOV pixel units; ~0.29 µm/px -> ≈ 15/30/60 µm
UM_PER_PX <- 0.29          # same approximation as 10.1.Statial.R / 10.2...R

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load object, derive feature columns (same as 10.3.SPIAT_tumour_border.R)
# ══════════════════════════════════════════════════════════════════════════════
srt <- qs_read(file.path(path, "tumour_normal_anno_srt.qs2"))
if (length(SeuratObject::Layers(srt, assay = "Spatial")) > 1)
  srt <- JoinLayers(srt, assay = "Spatial")

meta <- srt@meta.data
meta$category_bin <- ifelse(meta$category == "DT", "DT",
                      ifelse(grepl("^CB", meta$category), "CB", NA))
# Pool Tumour subclones ("Tumour 1"/"Tumour 2") into a single reference level.
meta$cell_type <- ifelse(grepl("^Tumour", meta$cell_type), "Tumour", meta$cell_type)

coords <- GetTissueCoordinates(srt, which = "centroids")
idx    <- match(colnames(srt), coords$cell)
stopifnot(!anyNA(idx))
xy       <- as.matrix(coords[idx, c("x", "y")])
countmat <- GetAssayData(srt, assay = "Spatial", layer = "counts")

cat("  ", i, ":", ncol(srt), "cells,", sum(meta$cell_type == "Tumour"), "Tumour cells\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Build SPIAT-formatted SpatialExperiment
# ══════════════════════════════════════════════════════════════════════════════
formatted_image <- format_image_to_spe(
  format           = "general",
  intensity_matrix = as.matrix(countmat),
  phenotypes       = meta$cell_type,
  coord_x          = xy[, 1],
  coord_y          = xy[, 2]
)
colData(formatted_image)$cell_type        <- meta$cell_type
colData(formatted_image)$final_annotation <- meta$final_annotation
colData(formatted_image)$category         <- meta$category
colData(formatted_image)$category_bin     <- meta$category_bin
colData(formatted_image)$slide            <- i

# ══════════════════════════════════════════════════════════════════════════════
# 3. Colocalization: CIN, Mixing Score, Cross-K — per non-Tumour cell type
# ══════════════════════════════════════════════════════════════════════════════
non_tumour_types <- setdiff(unique(meta$cell_type), "Tumour")

cin_rows   <- list()
ms_rows    <- list()
crossk_rows <- list()

for (ct in non_tumour_types) {

  n_target <- sum(meta$cell_type == ct)
  cat(sprintf("[%s] %s: %d cells\n", i, ct, n_target))

  # -- CIN + Mixing score, at each radius --
  for (r in radii) {
    cin <- tryCatch(
      average_percentage_of_cells_within_radius(
        spe_object = formatted_image, reference_celltype = "Tumour",
        target_celltype = ct, radius = r, feature_colname = "cell_type"),
      error = function(e) {
        message(sprintf("  [%s] CIN failed for %s at r=%d: %s", i, ct, r, conditionMessage(e)))
        NA_real_
      })
    cin_rows[[length(cin_rows) + 1]] <- data.frame(
      slide = i, cell_type = ct, radius_px = r, radius_um = r * UM_PER_PX, CIN_pct = cin)

    ms <- tryCatch(
      mixing_score_summary(
        spe_object = formatted_image, reference_celltype = "Tumour",
        target_celltype = ct, radius = r, feature_colname = "cell_type"),
      error = function(e) {
        message(sprintf("  [%s] mixing score failed for %s at r=%d: %s", i, ct, r, conditionMessage(e)))
        NULL
      })
    if (!is.null(ms))
      ms_rows[[length(ms_rows) + 1]] <- cbind(slide = i, radius_px = r, radius_um = r * UM_PER_PX, ms)
  }

  # -- Cross-K function, out to the largest radius --
  df_cross <- tryCatch(
    calculate_cross_functions(formatted_image, method = "Kcross",
                               cell_types_of_interest = c("Tumour", ct),
                               feature_colname = "cell_type", dist = max(radii)),
    error = function(e) {
      message(sprintf("  [%s] Cross-K failed for %s: %s", i, ct, conditionMessage(e)))
      NULL
    })

  if (!is.null(df_cross)) {
    auc      <- tryCatch(AUC_of_cross_function(df_cross), error = function(e) NA_real_)
    crossing <- tryCatch(crossing_of_crossK(df_cross), error = function(e) NA_real_)
    crossk_rows[[length(crossk_rows) + 1]] <- data.frame(
      slide = i, cell_type = ct, crossK_AUC = auc, crossK_crossing_frac = crossing)

    ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
    write.csv(as.data.frame(df_cross),
              file.path(crossk_dir, sprintf("%s_%s_crossK.csv", i, ct_safe)), row.names = FALSE)
    tryCatch({
      png(file.path(crossk_dir, sprintf("%s_%s_crossK.png", i, ct_safe)), width = 800, height = 600)
      plot(df_cross, main = sprintf("%s — Tumour vs %s Cross-K", i, ct))
      dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off() })
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Save outputs
# ══════════════════════════════════════════════════════════════════════════════
if (length(cin_rows))
  write.csv(bind_rows(cin_rows), file.path(tbl_dir, paste0(i, "_CIN.csv")), row.names = FALSE)
if (length(ms_rows))
  write.csv(bind_rows(ms_rows), file.path(tbl_dir, paste0(i, "_mixing_score.csv")), row.names = FALSE)
if (length(crossk_rows))
  write.csv(bind_rows(crossk_rows), file.path(tbl_dir, paste0(i, "_crossK_summary.csv")), row.names = FALSE)

cat(i, ": done —", length(non_tumour_types), "non-Tumour cell types analysed\n")
