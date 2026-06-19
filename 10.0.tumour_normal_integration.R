#!/usr/bin/env Rscript
# 10.0.tumour_normal_integration.R
# Combined tumour + clean-normal integration (run-once). Places ALL clean tumour
# cells and ALL clean normal cells into one batch-corrected Pearson PCA space.
#   - Normals come from the already-merged, already-cleaned 9.4.1 object
#     (Removed + tumour-module-positive already dropped; cell IDs slide-prefixed;
#     carries `final_annotation`).
#   - Tumour cells come from the 8 per-sample tumour/tumour_srt.qs2 objects and
#     are labelled `final_annotation = "Tumour"`.
# Keeps only the Spatial assay, runs do.pearson_pca(batch_variable = "slide"),
# stashes spatial centroids in @meta.data (FOV objects do not survive merge),
# saves the object, and rebuilds per-sample spatial DimPlots grouped by
# final_annotation (NA -> "Tumour") via ggplot on x_centroid/y_centroid.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
  library(ggdark, lib.loc = "~/R_Library/4.5")
})
source("~/VisHD/functions.R")   # do.pearson_pca()
options(future.globals.maxSize = 8 * 1024^3)

RES <- 1.0   # clustering resolution on the batch-corrected Pearson graph

# ── Paths ─────────────────────────────────────────────────────────────────────
normal_qs <- path.expand("~/VisHD/9.4.1.final_clear_normal_integration/normal_srt_final_anno.qs2")
out_dir   <- path.expand("~/VisHD/10.0.tumour_normal_integration")
png_dir   <- file.path(out_dir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

# ── Shared annotation palette (Tumour fixed to grey so colours match per panel) ─
final_palette <- function(v, grey_level = "Tumour") {
  others <- setdiff(sort(unique(v)), grey_level)
  pal <- setNames(as.vector(polychrome())[seq_along(others)], others)
  if (grey_level %in% v) pal <- c(pal, setNames("lightgrey", grey_level))
  pal
}

# ── 1a. Clean-normal mini object (from the merged 9.4.1 object) ───────────────
cat("Loading clean normals:", normal_qs, "\n")
normal_full <- qs_read(normal_qs)
normal_full <- JoinLayers(normal_full, assay = "Spatial")
if (!all(c("x_centroid", "y_centroid") %in% colnames(normal_full@meta.data)))
  stop("9.4.1 normal object lacks x_centroid/y_centroid in @meta.data — ",
       "cannot recover spatial coords (no FOV after merge).")
normal_meta <- normal_full@meta.data
normal_meta$final_annotation <- as.character(normal_meta$final_annotation)
normal_mini <- CreateSeuratObject(
  counts    = GetAssayData(normal_full, assay = "Spatial", layer = "counts"),
  meta.data = normal_meta, assay = "Spatial")
n_normal <- ncol(normal_mini)
cat("  normal cells:", n_normal, "\n")
rm(normal_full); gc()

# ── 1b. Per-sample tumour mini objects ────────────────────────────────────────
paths  <- system("realpath ~/VisHD/LUT-245-*/tumour/tumour_srt.qs2", intern = TRUE)
slides <- basename(dirname(dirname(paths)))
cat("\nLoading", length(paths), "tumour objects\n")

tumour_minis <- lapply(seq_along(paths), function(k) {
  cat("  Loading", slides[k], "\n")
  srt_full <- qs_read(paths[k])
  srt_full <- UpdateSeuratObject(srt_full)
  srt_full <- JoinLayers(srt_full, assay = "Spatial")
  counts <- GetAssayData(srt_full, assay = "Spatial", layer = "counts")
  meta   <- srt_full@meta.data
  meta$slide <- slides[k]
  coords <- GetTissueCoordinates(srt_full, which = "centroids")
  meta$x_centroid <- coords[rownames(meta), "x"]
  meta$y_centroid <- coords[rownames(meta), "y"]
  meta$final_annotation <- "Tumour"
  srt_mini <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
  # Prefix barcodes with the slide to match the 9.4.1 normal naming convention.
  srt_mini <- RenameCells(srt_mini, add.cell.id = slides[k])
  cat("    ", slides[k], ":", ncol(srt_mini), "tumour cells\n")
  rm(srt_full); gc()
  srt_mini
})
n_tumour <- sum(vapply(tumour_minis, ncol, numeric(1)))
cat("  total tumour cells:", n_tumour, "\n")

# ── 2. Merge (IDs already globally unique + slide-prefixed) ───────────────────
srt <- merge(normal_mini, y = tumour_minis)
DefaultAssay(srt) <- "Spatial"
srt <- JoinLayers(srt, assay = "Spatial")
cat("\nMerged:", ncol(srt), "cells (normal", n_normal, "+ tumour", n_tumour, ")\n")
stopifnot(ncol(srt) == n_normal + n_tumour)

# ── 3. Keep only the Spatial assay; clear stale reductions ────────────────────
srt <- DietSeurat(srt, assays = "Spatial")

# ── 4. Batch-corrected Pearson PCA (by slide) + UMAP + clusters ───────────────
srt <- do.pearson_pca(srt, batch_variable = "slide", assay = "Spatial",
                      resolution = RES)
Idents(srt) <- "pearson_clusters_batch"
n_clu <- nlevels(factor(srt$pearson_clusters_batch))
cat("Batch-corrected Pearson clustering (res", RES, "):", n_clu, "clusters\n")

# ── 5. Finalise final_annotation (residual NA -> "Tumour") ────────────────────
fa <- as.character(srt$final_annotation)
fa[is.na(fa)] <- "Tumour"
srt$final_annotation <- factor(fa)
cat("\nfinal_annotation:\n")
print(table(srt$final_annotation, useNA = "ifany"))

# ── 6. Save object ────────────────────────────────────────────────────────────
qs_save(srt, file.path(out_dir, "integrated_pearson_srt.qs2"))
cat("\nSaved integrated_pearson_srt.qs2\n")

# ── 7. Embedding DimPlots (pearsonbatchumap) ──────────────────────────────────
cols = c("#8DD3C7","#FFFFB3","#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
         "#FCCDE5", "#D9D9D9", "red")
anno_pal <- setNames(cols[seq_along(levels(srt$final_annotation))], levels(srt$final_annotation))

p_anno <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "final_annotation",
                  cols = anno_pal) +
  ggtitle("Tumour + normal integration — final annotation") +
  theme(legend.text = element_text(size = 7))
ggsave(file.path(png_dir, "1_DimPlot_final_annotation.png"), p_anno,
       width = 9, height = 7, dpi = 300)

p_clu <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "pearson_clusters_batch",
                 label = TRUE, label.size = 3, repel = TRUE,
                 cols = as.vector(polychrome())) +
  ggtitle(sprintf("Pearson batch clusters (res %.1f, n=%d)", RES, n_clu)) +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "2_DimPlot_clusters.png"), p_clu,
       width = 8, height = 7, dpi = 300)

p_slide <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "slide",
                   cols = as.vector(polychrome())) +
  ggtitle("Slide (batch-mixing check)") +
  theme(legend.text = element_text(size = 7))
ggsave(file.path(png_dir, "3_DimPlot_slide.png"), p_slide,
       width = 8, height = 7, dpi = 300)

p_category <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "category",
                   cols = as.vector(polychrome())) +
  ggtitle("Category (batch-mixing check)") +
  theme(legend.text = element_text(size = 7))
ggsave(file.path(png_dir, "3.2_DimPlot_category.png"), p_category,
       width = 8, height = 7, dpi = 300)
# ── 8. Per-sample spatial DimPlots (ggplot on centroids; no FOV after merge) ──
spatial_df <- srt@meta.data[, c("slide", "x_centroid", "y_centroid", "final_annotation")]
spatial_df$final_annotation <- factor(spatial_df$final_annotation,
                                      levels = names(anno_pal))

for (s in sort(unique(spatial_df$slide))) {
  d <- spatial_df[spatial_df$slide == s, ]
  p <- ggplot(d, aes(x_centroid, y_centroid, colour = final_annotation)) +
    geom_point(size = 0.3) +
    scale_y_reverse() +  
    scale_colour_manual(values = anno_pal, drop = FALSE, name = "Annotation") +
    coord_fixed() + ggtitle(s) + theme_classic()  +
    guides(colour = guide_legend(override.aes = list(size = 2)))
  ggsave(file.path(png_dir, paste0("spatial_", s, ".png")), p,
         width = 8, height = 7, dpi = 350)
}

p_all <- ggplot(spatial_df, aes(x_centroid, y_centroid, colour = final_annotation)) +
  geom_point(size = 0.15) +
  scale_y_reverse() +  
  scale_colour_manual(values = anno_pal, drop = FALSE, name = "Annotation") +
  facet_wrap(~ slide, scales = "free", ncol = 4) +
  theme_classic()  +theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(png_dir, "9_spatial_all_samples.png"), p_all,
       width = 20, height = 10, dpi = 300, limitsize = FALSE)

cat("\nDone. Outputs in", out_dir, "\n")
