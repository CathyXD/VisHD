#!/usr/bin/env Rscript
# 8.6.final_clear_normal_persample.R
# Per-sample propagation of the final normal annotation (array 1-8). Reads the
# refined object written by 8.5.4.refined_subtype_recluster.R, attaches its
# general-layer `final_annotation` AND its `refined_subtype_pearson` call to
# that sample's normal_srt.qs2 by bare-barcode match, drops the
# unmatched/"Removed" cells, re-embeds the cleaned subset (SpaNorm + Pearson
# PCA + BANKSY, as in 4.2.normal_split.R), writes normal_anno_srt.qs2 in the
# sample's normal/ dir, and saves an ImageDimPlot of each annotation in tissue
# space. If normal_srt.qs2 is missing (e.g. cleaned up from scratch after a
# prior successful run) but this sample's own normal_anno_srt.qs2 still
# exists, that object is updated in place instead — its embedding/clusters
# are reused as-is and only the two annotation columns are refreshed.
#
#   Rscript 8.6.final_clear_normal_persample.R <sample-index 1-8>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  library(pals)
  library(qs2)
})
source("~/VisHD/functions.R")   # do.spanorm(), do.pearson_pca()
options(future.globals.maxSize = 8 * 1024^3)

refined_srt_path <- path.expand("~/VisHD/8.5.4.refined_subtype_recluster/normal_srt_refined_recluster.qs2")

# ── DimPlot palette helper (removed/unknown in grey) ──────────────────────────
final_palette <- function(v, grey_level = "Removed") {
  others <- setdiff(sort(unique(v)), grey_level)
  pal <- setNames(as.vector(polychrome())[seq_along(others)], others)
  if (grey_level %in% v) pal <- c(pal, setNames("lightgrey", grey_level))
  pal
}

# ── CLI arg ───────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 8.6.final_clear_normal_persample.R <sample-index 1-8>")
sample_idx <- as.numeric(args[1])

if (!file.exists(refined_srt_path))
  stop("Refined subtype object not found: ", refined_srt_path,
       "\nRun the integrative stage first: Rscript 8.5.4.refined_subtype_recluster.R")

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[sample_idx]
i     <- basename(path)
cat("working at", path, "\n")

normal_dir <- file.path(path, "normal")
setwd(normal_dir)

# normal_srt.qs2 (4.2's output) is the normal re-embedding path; if it has
# since been cleaned up from scratch but this sample's normal_anno_srt.qs2
# (this script's own prior output) still exists, fall back to updating that
# object's meta.data in place — it already carries a valid SpaNorm/Pearson/
# BANKSY embedding, so there is nothing to re-embed.
needs_reembed <- file.exists("normal_srt.qs2")
if (needs_reembed) {
  normal_srt <- qs_read("normal_srt.qs2")
  cat("  loaded", ncol(normal_srt), "cells (from normal_srt.qs2)\n")
} else if (file.exists("normal_anno_srt.qs2")) {
  cat("  normal_srt.qs2 missing — updating existing normal_anno_srt.qs2 meta.data in place\n")
  normal_srt <- qs_read("normal_anno_srt.qs2")
  cat("  loaded", ncol(normal_srt), "cells (from normal_anno_srt.qs2, already embedded)\n")
} else {
  stop("Neither normal_srt.qs2 nor normal_anno_srt.qs2 found in ", normal_dir)
}

# Attach this slide's general-layer + refined-subtype annotation by
# bare-barcode match (cell_ID strips its "{slide}_" prefix, as in 8.4/8.5.4).
cat("Loading refined annotation from", refined_srt_path, "\n")
refined_srt <- qs_read(refined_srt_path)
anno <- data.frame(
  slide                   = as.character(refined_srt$slide),
  barcode                 = str_remove(as.character(refined_srt$cell_ID),
                                      paste0("^", refined_srt$slide, "_")),
  final_annotation        = as.character(refined_srt$final_annotation),
  refined_subtype_pearson = as.character(refined_srt$refined_subtype_pearson),
  stringsAsFactors = FALSE
)
rm(refined_srt); gc()

anno <- anno[anno$slide == i, ]
if (nrow(anno) == 0) stop("No final-annotation rows for slide ", i, " in ", refined_srt_path)
anno_map     <- setNames(anno$final_annotation,        anno$barcode)
subtype_map  <- setNames(anno$refined_subtype_pearson,  anno$barcode)
normal_srt$final_annotation        <- unname(anno_map[colnames(normal_srt)])
normal_srt$refined_subtype_pearson <- unname(subtype_map[colnames(normal_srt)])
cat("  matched", sum(!is.na(normal_srt$final_annotation)), "/", ncol(normal_srt),
    "cells to final annotation\n")

# Clear removed / unmatched cells. (9.4.1 already excludes Removed from the CSV,
# so unmatched cells here are the dropped clusters; both are cleared.)
keep <- colnames(normal_srt)[!is.na(normal_srt$final_annotation) &
                             normal_srt$final_annotation != "Removed"]
normal_srt <- subset(normal_srt, cells = keep)
cat("  kept", ncol(normal_srt), "cells after clearing Removed/unmatched\n")
print(table(normal_srt$final_annotation, useNA = "ifany"))

if (needs_reembed) {
  # Re-embed the cleaned subset (SpaNorm + Pearson PCA + BANKSY) — same as 4.2.
  normal_srt <- do.spanorm(normal_srt)
  normal_srt <- do.pearson_pca(normal_srt)
  normal_srt <- FindClusters(normal_srt, resolution = 0.8, algorithm = 4)
  cat("  SpaNorm + Pearson PCA done\n")

  normal_srt <- SeuratWrappers::RunBanksy(normal_srt, lambda = 0.2, verbose = TRUE,
                                          use_agf = TRUE, assay = "SpaNorm", slot = "data",
                                          k_geom = c(15), assay_name = "BANKSY_0.2")
  normal_srt <- RunPCA(normal_srt, npcs = 30, features = rownames(normal_srt),
                       reduction.name = "banksy0.2.pca")
  normal_srt <- RunUMAP(normal_srt, dims = 1:20, reduction = "banksy0.2.pca",
                        reduction.name = "banksy0.2.umap")
  normal_srt <- FindNeighbors(normal_srt, reduction = "banksy0.2.pca", dims = 1:20)
  normal_srt <- FindClusters(normal_srt, resolution = 1, algorithm = 4)
  cat("  BANKSY done\n")
} else {
  cat("  Skipping re-embedding — reusing normal_anno_srt.qs2's existing embedding\n")
}

# Re-assert both annotations as clean factors and set the general layer active.
normal_srt$final_annotation        <- factor(normal_srt$final_annotation)
normal_srt$refined_subtype_pearson <- factor(normal_srt$refined_subtype_pearson)
Idents(normal_srt) <- "final_annotation"

qs_save(normal_srt, "normal_anno_srt.qs2")

dir.create("png", showWarnings = FALSE, recursive = TRUE)
pal <- final_palette(as.character(normal_srt$final_annotation))
ggsave("png/final_annotation_ImageDimPlot.png",
       plot = ImageDimPlot(normal_srt, group.by = "final_annotation",
                           cols = pal, border.color = "#00000000", size = 0.3),
       width = 8, height = 7, dpi = 350)

subtype_v   <- as.character(normal_srt$refined_subtype_pearson)
subtype_v[is.na(subtype_v)] <- "Unknown"
subtype_pal <- final_palette(subtype_v, grey_level = "Unknown")
ggsave("png/refined_subtype_pearson_ImageDimPlot.png",
       plot = ImageDimPlot(normal_srt, group.by = "refined_subtype_pearson",
                           cols = subtype_pal, border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 350)

cat("===================== ", i, " final annotation done ====================\n")
