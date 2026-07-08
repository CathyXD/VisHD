#!/usr/bin/env Rscript
# 8.6.final_clear_normal_persample.R
# Per-sample propagation of the final normal annotation (array 1-8). Reads the
# final_annotation.csv written by 9.4.1, attaches it to that sample's
# normal_srt.qs2 by bare-barcode match, drops the "Removed"/unmatched cells,
# re-embeds the cleaned subset (SpaNorm + Pearson PCA + BANKSY, as in
# 4.2.normal_split.R), writes normal_anno_srt.qs2 in the sample's normal/ dir,
# and saves an ImageDimPlot of the final annotation in tissue space.
#
#   Rscript 8.6.final_clear_normal_persample.R <sample-index 1-8>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
})
source("~/VisHD/functions.R")   # do.spanorm(), do.pearson_pca()
options(future.globals.maxSize = 8 * 1024^3)

anno_csv <- path.expand("~/VisHD/8.4.final_clear_normal_integration/final_annotation.csv")

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

if (!file.exists(anno_csv))
  stop("Final annotation CSV not found: ", anno_csv,
       "\nRun the integrative stage first: Rscript 8.4.final_clear_normal_integration.R")

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[sample_idx]
i     <- basename(path)
cat("working at", path, "\n")

normal_dir <- file.path(path, "normal")
setwd(normal_dir)
normal_srt <- qs_read("normal_srt.qs2")
cat("  loaded", ncol(normal_srt), "cells\n")

# Attach this slide's final annotation by bare-barcode rowname match.
anno <- read.csv(anno_csv, check.names = FALSE, stringsAsFactors = FALSE)
anno <- anno[anno$slide == i, ]
if (nrow(anno) == 0) stop("No final-annotation rows for slide ", i, " in ", anno_csv)
anno_map <- setNames(anno$final_annotation, anno$barcode)
normal_srt$final_annotation <- unname(anno_map[colnames(normal_srt)])
cat("  matched", sum(!is.na(normal_srt$final_annotation)), "/", ncol(normal_srt),
    "cells to final annotation\n")

# Clear removed / unmatched cells. (9.4.1 already excludes Removed from the CSV,
# so unmatched cells here are the dropped clusters; both are cleared.)
keep <- colnames(normal_srt)[!is.na(normal_srt$final_annotation) &
                             normal_srt$final_annotation != "Removed"]
normal_srt <- subset(normal_srt, cells = keep)
cat("  kept", ncol(normal_srt), "cells after clearing Removed/unmatched\n")
print(table(normal_srt$final_annotation, useNA = "ifany"))

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

# Re-assert the final annotation as a clean factor and active identity.
normal_srt$final_annotation <- factor(normal_srt$final_annotation)
Idents(normal_srt) <- "final_annotation"

qs_save(normal_srt, "normal_anno_srt.qs2")

pal <- final_palette(as.character(normal_srt$final_annotation))
ggsave("png/final_annotation_ImageDimPlot.png",
       plot = ImageDimPlot(normal_srt, group.by = "final_annotation",
                           cols = pal, border.color = "#00000000", size = 0.3),
       width = 8, height = 7, dpi = 350)

cat("===================== ", i, " final annotation done ====================\n")
