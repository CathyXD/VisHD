#!/usr/bin/env Rscript
# 9.4.tumour_normal_integration.R
# Combined tumour + clean-normal integration (run-once). Places ALL clean tumour
# cells and ALL clean normal cells into one batch-corrected Pearson PCA space.
#   - Aggregates the 8 per-sample LUT-245-XX/tumour_normal_anno_srt.qs2 objects
#     written by 9.1.per_sample_tumour_normal.R. Each carries both compartments:
#     normal cells hold their `final_annotation`; tumour cells have NA (filled to
#     "Tumour" below). The `cell`/`compartment`/`cell_type` columns ride along.
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
  library(UCell, lib.loc = "~/R_Library/4.5")
})
source("~/VisHD/functions.R")   # do.pearson_pca()
options(future.globals.maxSize = 8 * 1024^3)

RES <- 1.0   # clustering resolution on the batch-corrected Pearson graph

# ── Paths ─────────────────────────────────────────────────────────────────────
out_dir   <- path.expand("~/VisHD/10.0.tumour_normal_integration")
png_dir   <- file.path(out_dir, "png")
integrated_file <- file.path(out_dir, "integrated_pearson_srt.qs2")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

# ── Shared annotation palette (Tumour fixed to grey so colours match per panel) ─
final_palette <- function(v, grey_level = "Tumour") {
  others <- setdiff(sort(unique(v)), grey_level)
  pal <- setNames(as.vector(polychrome())[seq_along(others)], others)
  if (grey_level %in% v) pal <- c(pal, setNames("lightgrey", grey_level))
  pal
}

if (!file.exists(integrated_file)) {
# ── 0. Reuse the existing integrated object (skip rebuild) ────────────────────
cat("Loading existing integrated object:", integrated_file, "\n")
srt <- qs_read(integrated_file)
} else {
# ── 1. Per-sample tumour+normal mini objects (from 10.0 per-sample output) ────
# Each LUT-245-XX/tumour_normal_anno_srt.qs2 carries both compartments with a
# single-sample FOV; keep only Spatial counts + the centroids/annotation we need.
paths  <- system("realpath ~/VisHD/LUT-245-*/tumour_normal_anno_srt.qs2", intern = TRUE)
slides <- basename(dirname(paths))
cat("\nLoading", length(paths), "per-sample tumour+normal objects\n")

minisres <- lapply(seq_along(paths), function(k) {
  cat("  Loading", slides[k], "\n")
  srt_full <- qs_read(paths[k])
  counts <- GetAssayData(srt_full, assay = "Spatial", layer = "counts")
  meta   <- srt_full@meta.data
  meta$slide <- slides[k]
  coords <- GetTissueCoordinates(srt_full, which = "centroids")
  idx <- match(rownames(meta), coords$cell)
  stopifnot(!anyNA(idx))
  meta$x_centroid <- coords$x[idx]
  meta$y_centroid <- coords$y[idx]

  # final_annotation: normal cell-types; NA on tumour cells (filled -> "Tumour" below)
  meta$final_annotation <- as.character(meta$final_annotation)
  srt_mini <- CreateSeuratObject(counts = counts,  assay = "Spatial")
  # Prefix barcodes with the slide so IDs are globally unique across samples.
  srt_mini <- RenameCells(srt_mini, add.cell.id = slides[k])
  cat("    ", slides[k], ":", ncol(srt_mini), "cells\n")
  res <- list(srt = srt_mini, meta = meta)
  rm(srt_full); gc()
  return(res)
})

minis <- lapply(minisres, `[[`, "srt")
metaList <- lapply(minisres, `[[`, "meta")
metaList <- lapply(metaList, function(x) x[, colnames(metaList[[1]])])  # ensure same column order
meta <- do.call(rbind, metaList)
rownames(meta) <- paste(meta$slide, meta$cell, sep = "_")
n_cells <- sum(vapply(minis, ncol, numeric(1)))
cat("  total cells:", n_cells, "\n")

# ── 2. Merge (IDs slide-prefixed) ─────────────────────────────────────────────
srt <- merge(minis[[1]], y = minis[-1])
DefaultAssay(srt) <- "Spatial"
srt <- JoinLayers(srt, assay = "Spatial")
cat("\nMerged:", ncol(srt), "cells\n")
stopifnot(ncol(srt) == n_cells)

# ── 3. Keep only the Spatial assay; clear stale reductions ────────────────────
srt <- DietSeurat(srt, assays = "Spatial")
srt <- AddMetaData(srt, meta)
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
qs_save(srt, integrated_file)
cat("\nSaved integrated_pearson_srt.qs2\n")
}

Idents(srt) <- "pearson_clusters_batch"
n_clu <- nlevels(factor(srt$pearson_clusters_batch))

# ── 7. Embedding DimPlots (pearsonbatchumap) ──────────────────────────────────
cols = c("red","gold", "#8DD3C7","#FFFFB3","#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
         "#FCCDE5", "#D9D9D9")
anno_pal <- setNames(cols[seq_along(unique(srt$cell_type))], unique(srt$cell_type))

p_anno <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "cell_type",
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
spatial_df <- srt@meta.data[, c("slide", "x_centroid", "y_centroid", "final_annotation", "cell_type")]
spatial_df$final_annotation <- factor(spatial_df$cell_type,
                                      levels = names(anno_pal))

for (s in sort(unique(spatial_df$slide))) {
  d <- spatial_df[spatial_df$slide == s, ]
  p <- ggplot(d, aes(x_centroid, y_centroid, colour = cell_type)) +
    geom_point(size = 0.3) +
    scale_y_reverse() +  
    scale_colour_manual(values = anno_pal, drop = FALSE, name = "Annotation") +
    coord_fixed() + ggtitle(s) + theme_classic()  +
    guides(colour = guide_legend(override.aes = list(size = 2)))
  ggsave(file.path(png_dir, paste0("spatial_", s, ".png")), p,
         width = 8, height = 7, dpi = 350)
}

p_all <- ggplot(spatial_df, aes(x_centroid, y_centroid, colour = cell_type)) +
  geom_point(size = 0.15) +
  scale_y_reverse() +  
  scale_colour_manual(values = anno_pal, drop = FALSE, name = "Annotation") +
  facet_wrap(~ slide, scales = "free", ncol = 4) +
  theme_classic()  +theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(png_dir, "9_spatial_all_samples.png"), p_all,
       width = 20, height = 10, dpi = 300, limitsize = FALSE)

# ── 9. Module-group (G1/G2/G3 signature) visualisation ────────────────────────
# Resolve the single "Tumour" bucket into its per-cell signature label using the
# 6.2.3 binarisation output (metas$Module_group), and score the three groupdeg
# signatures on the integrated object. Normals collapse to a single "Normal".
metas    <- readRDS("~/VisHD/6.2.3.signature_analysis/binarisation/metas.Rds")
groupdeg <- readRDS(paste0("~/VisHD/6.2archetype_downstream_tumour/archetype_module/",
                           "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))

# barcode key in the integrated object is slide_<cell>; metas carries slide + cell
mg_lookup <- setNames(as.character(metas$Module_group),
                      paste0(metas$slide, "_", metas$cell))
is_tum    <- srt$compartment == "Tumour"
module_anno          <- rep("Normal", ncol(srt))
module_anno[is_tum]  <- mg_lookup[colnames(srt)[is_tum]]
n_unmatched <- sum(is_tum & is.na(module_anno))
cat(sprintf("\nTumour module_anno match: %d/%d (%.1f%%); %d unmatched -> 'Neg'\n",
            sum(is_tum) - n_unmatched, sum(is_tum),
            100 * (sum(is_tum) - n_unmatched) / sum(is_tum), n_unmatched))
module_anno[is.na(module_anno)] <- "Neg"

# rebuild 6.2.3 group_levels/group_pal (derive levels from groupdeg so they match
# metas); make the palette order-invariant so combos like "G3/G1" are coloured
labs <- names(groupdeg)
group_combos <- unlist(lapply(seq_along(labs), function(k)
  combn(labs, k, FUN = function(x) paste(x, collapse = "/"))))
group_levels <- c("Neg", group_combos)
group_pal <- c("Neg"      = "lightblue",
               "G1"       = "red",
               "G2"       = "gold",
               "G3"       = "royalblue",
               "G1/G2"    = "orange",
               "G1/G3"    = "purple",
               "G2/G3"    = "green",
               "G1/G2/G3" = "grey")
canon   <- function(x) vapply(strsplit(x, "/"),
                              function(p) paste(sort(p), collapse = "/"), character(1))
present <- levels(droplevels(factor(module_anno)))
mg_pal  <- setNames(group_pal[canon(present)], present)
mg_pal["Normal"] <- "lightpink"

mg_levels <- c("Normal", group_levels)
srt$module_anno <- factor(module_anno, levels = mg_levels[mg_levels %in% module_anno])
cat("\nmodule_anno:\n"); print(table(srt$module_anno, useNA = "ifany"))

# UMAP DimPlot by module_anno
p_mg <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "module_anno",
                cols = mg_pal) +
  ggtitle("Tumour signature module group + Normal") +
  theme(legend.text = element_text(size = 7))
ggsave(file.path(png_dir, "10_DimPlot_module_group.png"), p_mg,
       width = 9, height = 7, dpi = 300)

# Spatial DimPlot by module_anno (ggplot on centroids; no FOV after merge)
mg_df <- srt@meta.data[, c("slide", "x_centroid", "y_centroid", "module_anno")]
for (s in sort(unique(mg_df$slide))) {
  d <- mg_df[mg_df$slide == s, ]
  p <- ggplot(d, aes(x_centroid, y_centroid, colour = module_anno)) +
    geom_point(size = 0.3) +
    scale_y_reverse() +
    scale_colour_manual(values = mg_pal, drop = FALSE, name = "Module group") +
    coord_fixed() + ggtitle(s) + theme_classic() +
    guides(colour = guide_legend(override.aes = list(size = 2)))
  ggsave(file.path(png_dir, paste0("spatial_modulegroup_", s, ".png")), p,
         width = 8, height = 7, dpi = 350)
}
p_mg_all <- ggplot(mg_df, aes(x_centroid, y_centroid, colour = module_anno)) +
  geom_point(size = 0.15) +
  scale_y_reverse() +
  scale_colour_manual(values = mg_pal, drop = FALSE, name = "Module group") +
  facet_wrap(~ slide, scales = "free", ncol = 4) +
  theme_classic() + theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(png_dir, "11_spatial_module_group_all.png"), p_mg_all,
       width = 20, height = 10, dpi = 300, limitsize = FALSE)

# Score the three groupdeg signatures on the integrated object (counts-only assay)
DefaultAssay(srt) <- "Spatial"
srt   <- NormalizeData(srt, verbose = FALSE)
feats <- lapply(groupdeg, function(g) intersect(g, rownames(srt)))
srt   <- AddModuleScore_UCell(srt, features = feats, name = "_GDmod")
mod_cols <- paste0(names(groupdeg), "_GDmod")            # module_G1/G2/G3
# FeaturePlot (UMAP) — steelblue/white/indianred gradient2 centred at 0
p_feat <- FeaturePlot(srt, reduction = "pearsonbatchumap", features = mod_cols,
                      ncol = length(mod_cols), order = TRUE) &
  scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                         midpoint = 0)
ggsave(file.path(png_dir, "12_FeaturePlot_UCellmodules.png"), p_feat,
       width = 5 * length(mod_cols), height = 5, dpi = 300)

# Spatial FeaturePlot per module (ggplot on centroids, faceted by slide)
for (m in mod_cols) {
  d <- srt@meta.data[, c("slide", "x_centroid", "y_centroid", m)]
  p <- ggplot(d, aes(x_centroid, y_centroid, colour = .data[[m]])) +
    geom_point(size = 0.15) +
    scale_y_reverse() +
    scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                           midpoint = 0, name = m) +
    facet_wrap(~ slide, scales = "free", ncol = 4) +
    ggtitle(m) + theme_classic() + theme(aspect.ratio = 1)
  ggsave(file.path(png_dir, paste0("13_spatial_UCell_", m, ".png")), p,
         width = 20, height = 10, dpi = 300, limitsize = FALSE)
}

srt <- AddModuleScore(srt, features = feats, name = "Module", assay = "Spatial")
colnames(srt@meta.data)[grep("Module", colnames(srt@meta.data))] <- paste0("Module_", names(groupdeg))  # Module_G1/G2/G3
mod_cols <- paste0("Module_", names(groupdeg)) 

# FeaturePlot (UMAP) — steelblue/white/indianred gradient2 centred at 0
p_feat <- FeaturePlot(srt, reduction = "pearsonbatchumap", features = mod_cols,
                      ncol = length(mod_cols), order = TRUE) &
  scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                         midpoint = 0)
ggsave(file.path(png_dir, "12_FeaturePlot_modules.png"), p_feat,
       width = 5 * length(mod_cols), height = 5, dpi = 300)

# Spatial FeaturePlot per module (ggplot on centroids, faceted by slide)
for (m in mod_cols) {
  d <- srt@meta.data[, c("slide", "x_centroid", "y_centroid", m)]
  p <- ggplot(d, aes(x_centroid, y_centroid, colour = .data[[m]])) +
    geom_point(size = 0.15) +
    scale_y_reverse() +
    scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                           midpoint = 0, name = m) +
    facet_wrap(~ slide, scales = "free", ncol = 4) +
    ggtitle(m) + theme_classic() + theme(aspect.ratio = 1)
  ggsave(file.path(png_dir, paste0("13_spatial_", m, ".png")), p,
         width = 20, height = 10, dpi = 300, limitsize = FALSE)
}

# Persist new fields (module_anno + module scores) onto the saved object
qs_save(srt, file.path(out_dir, "integrated_pearson_srt.qs2"))
cat("\nRe-saved integrated_pearson_srt.qs2 with module_anno + module scores\n")

cat("\nDone. Outputs in", out_dir, "\n")
