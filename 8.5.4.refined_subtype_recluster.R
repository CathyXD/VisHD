#!/usr/bin/env Rscript
# 8.5.4.refined_subtype_recluster.R
# Maps 8.5.3's per-cell refined subtype/status calls back onto the full
# 8-sample normal-cell object, visualizes them on the existing clustering,
# drops cells flagged contaminant/unresolved by either method, re-preprocesses
# (SCTransform -> batch-corrected Pearson PCA/UMAP, mirroring
# 8.4.final_clear_normal_integration.R Steps 3-4) and visualizes the new
# clustering + final_annotation + refined subtypes + curated marker genes.
#
# `status_pearson`/`status_rpca` == "contaminant" or "no-identity" (junk) on
# EITHER method flags a cell for removal (union — the more conservative QC
# call, since contamination missed by one clustering method is often caught
# by the other per METHOD_COMPARISON_summary.md).
#
# Loads:
#   8.4.final_clear_normal_integration/normal_srt_final_anno.qs2
#   8.5.normal_cell_subtypes/refined_tables/integrated_refined_annotation.csv
# Writes:
#   8.5.4.refined_subtype_recluster/normal_srt_refined_recluster.qs2
#   8.5.4.refined_subtype_recluster/png/*.png

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
})
source("~/VisHD/functions.R")   # do.pearson_pca()
stopifnot(requireNamespace("glmGamPoi", quietly = TRUE)) # fast SCTransform
options(future.globals.maxSize = 8 * 1024^3)

in_srt      <- path.expand("~/VisHD/8.4.final_clear_normal_integration/normal_srt_final_anno.qs2")
refined_csv <- path.expand("~/VisHD/8.5.normal_cell_subtypes/refined_tables/integrated_refined_annotation.csv")
out_dir     <- path.expand("~/VisHD/8.5.4.refined_subtype_recluster")
png_dir     <- file.path(out_dir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

RES <- 1.5   # clustering resolution on the new batch-corrected Pearson graph

marker_list <- list(
  CAF           = c("COL1A1", "COL1A2", "DCN", "LUM", "SFRP4"),
  Smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES"),
  Pericyte      = c("RGS5", "NOTCH3", "MCAM", "BCAM", "PDGFRB"),
  Macrophages   = c("CD68", "CD163", "CSF1R", "C1QA", "MS4A7"),
  T_cells       = c("TRAC", "TRBC2", "IL7R", "IL32", "IKZF1"),
  B_cells       = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  Plasma        = c("MZB1", "JCHAIN", "XBP1", "IGHG1", "TXNDC5"),
  Epithelial    = c("KRT5", "KRT15", "KRT17", "TACSTD2", "CLDN4"),
  SVEC          = c("SEMG1", "SEMG2", "MUC6", "PIP", "PATE1"),
  Neurons       = c("S100B", "NRXN1", "SCN7A", "PMP22", "PTGDS")
)

# "Unknown" (NA) always drawn light grey; other levels use a polychrome palette.
plot_anno <- function(srt, group, title, reduction = "pearsonbatchumap", label = TRUE) {
  v <- as.character(srt@meta.data[[group]]); v[is.na(v)] <- "Unknown"
  srt@meta.data[[group]] <- v
  others <- setdiff(sort(unique(v)), "Unknown")
  pal <- setNames(as.vector(polychrome())[seq_along(others)], others)
  if ("Unknown" %in% v) pal <- c(pal, Unknown = "lightgrey")
  DimPlot(srt, reduction = reduction, group.by = group, cols = pal,
          label = label, label.size = 2, repel = TRUE) +
    ggtitle(title) + theme(legend.text = element_text(size = 6))
}

# ── 1. Load + map back refined subtypes from 8.5.3 ──────────────────────────
cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
cat("  ", ncol(srt), "cells\n")

refined <- read.csv(refined_csv, stringsAsFactors = FALSE)
idx <- match(colnames(srt), refined$barcode)
n_unmatched <- sum(is.na(idx))
if (n_unmatched > 0)
  cat(sprintf("WARNING: %d/%d cells have no match in %s (left as NA/Unknown)\n",
              n_unmatched, ncol(srt), refined_csv))

srt$refined_subtype_pearson <- refined$refined_subtype_pearson[idx]
srt$status_pearson          <- refined$status_pearson[idx]
srt$refined_subtype_rpca    <- refined$refined_subtype_rpca[idx]
srt$status_rpca             <- refined$status_rpca[idx]

# ── 2. Visualize refined subtypes on the existing clustering (pre-removal) ──
p_pearson <- plot_anno(srt, "refined_subtype_pearson", "Refined subtype (pearson) — before removal")
p_rpca    <- plot_anno(srt, "refined_subtype_rpca",    "Refined subtype (rpca) — before removal")
ggsave(file.path(png_dir, "1_refined_subtype_pearson_preremoval.png"), p_pearson,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(png_dir, "2_refined_subtype_rpca_preremoval.png"), p_rpca,
       width = 10, height = 8, dpi = 300, bg = "white")

# ── 3. Remove contaminant / unresolved cells (flagged by either method) ─────
bad_status <- c("no-identity", "unclassified")
is_bad <- (!is.na(srt$status_pearson) & srt$status_pearson %in% bad_status)
is_bad_pearson <- (!is.na(srt$status_pearson) & srt$status_pearson %in% bad_status)
srt$flagged_removal <- is_bad_pearson
plot_anno(srt, "status_pearson", "Flagged for status in pearson")
plot_anno(srt, "status_rpca",    "Flagged for status in rpca", reduction = "integrated.rpca_umap")
cat(sprintf("Removing %d / %d cells flagged contaminant/unresolved by either method\n",
            sum(is_bad), ncol(srt)))
srt <- subset(srt, cells = colnames(srt)[!is_bad])
cat("Kept", ncol(srt), "cells\n")

# ── 4. Re-preprocess: SCTransform -> HVGs -> batch-corrected Pearson PCA/UMAP
# (mirrors 8.4.final_clear_normal_integration.R Steps 3-4 on the cleaned set) ─
DefaultAssay(srt) <- "Spatial"
srt <- SCTransform(srt, assay = "Spatial", method = "glmGamPoi",
                   variable.features.n = 3000, verbose = FALSE)
hvg_sct <- VariableFeatures(srt, assay = "SCT")
spatial_genes <- rownames(GetAssayData(srt, assay = "Spatial", layer = "counts"))
hvg <- intersect(hvg_sct, spatial_genes)
VariableFeatures(srt, assay = "Spatial") <- hvg
cat("SCT HVGs:", length(hvg_sct), "-> usable on Spatial:", length(hvg), "\n")

srt <- do.pearson_pca(srt, batch_variable = "slide", assay = "Spatial",
                      find_hvgs = FALSE, reduction_prefix = "pearsonbatch",
                      clusters_col = "pearson_clusters_batch", resolution = RES)
n_clu <- nlevels(factor(srt$pearson_clusters_batch))
cat("New batch-corrected Pearson clustering (res", RES, "):", n_clu, "clusters\n")

# ── 5. Visualize new clustering, final_annotation, refined subtypes ─────────
p_clu   <- plot_anno(srt, "pearson_clusters_batch",
                     sprintf("New Pearson clusters (res %.1f, n=%d)", RES, n_clu))
p_final <- plot_anno(srt, "final_annotation",        "Final annotation")
p_sub_p <- plot_anno(srt, "refined_subtype_pearson", "Refined subtype (pearson)")
p_sub_r <- plot_anno(srt, "refined_subtype_rpca",    "Refined subtype (rpca)")

ggsave(file.path(png_dir, "3_new_pearson_clusters.png"),    p_clu,   width = 9,  height = 7, dpi = 300, bg = "white")
ggsave(file.path(png_dir, "4_final_annotation.png"),        p_final, width = 9,  height = 7, dpi = 300, bg = "white")
ggsave(file.path(png_dir, "5_refined_subtype_pearson.png"), p_sub_p, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(png_dir, "6_refined_subtype_rpca.png"),    p_sub_r, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(png_dir, "7_overlay_grid.png"),
       (p_clu | p_final) / (p_sub_p | p_sub_r), width = 20, height = 14, dpi = 300, bg = "white")

# ── 6. Feature plots of curated marker_list genes ────────────────────────────
for (grp in names(marker_list)) {
  feats <- intersect(marker_list[[grp]], rownames(srt))
  if (!length(feats)) {
    cat(sprintf("[%s] no marker genes present — skipping\n", grp))
    next
  }
  fp <- FeaturePlot(srt, feats, reduction = "pearsonbatchumap",
                    cols = c("white", "red"), order = TRUE) +
    plot_layout(ncol = 5)
  ggsave(file.path(png_dir, sprintf("marker_featureplot_%s.png", grp)),
         fp, width = 20, height = ceiling(length(feats) / 5) * 4,
         dpi = 300, limitsize = FALSE, bg = "white")
}

# ── 7. Save ───────────────────────────────────────────────────────────────────
qs_save(srt, file.path(out_dir, "normal_srt_refined_recluster.qs2"))
cat("\nDone. Outputs in", out_dir, "\n")
