#!/usr/bin/env Rscript
# 8.5.2.general_layer_analysis.R   (per-general_layer-group, array 1-5)
# General-layer (coarse compartment) counterpart of 8.5.additional_analysis.R,
# split out from that script's former Sections 7-9. general_layer subsetting
# is independent of the fine cell-type pipeline, and running both in one
# script meant a crash partway through the fine cell-type loop (e.g. the
# IntegrateLayers "anchor cells < k.weight" RPCA error) prevented general_layer
# from ever being attempted in the same run.
#
# Split into an array job (one task per general_layer group) because looping
# over all 5 groups serially in a single job kept hitting the SLURM time
# limit partway through (Stromal alone is 80k+ cells) before reaching the
# smaller groups or the visualization section. Each array task now runs the
# full subset -> recluster/DE -> visualization pipeline for exactly one
# group, so a slow/large group (Stromal) no longer blocks the others.
#
# Requires 8.5.additional_analysis.R to have run at least once first — it
# reads `general_layer` from the SAME normal_srt_final_anno.qs2 object that
# script derives (from final_annotation) and saves back.
#
#   Rscript 8.5.2.general_layer_analysis.R <array-index 1-5>
#   (index maps to levels(srt$general_layer), i.e. Stromal/Myeloid/Lymphoid/
#   Epithelial/Neural in whatever order that factor was built with)
#
# Loads:
#   8.4.final_clear_normal_integration/normal_srt_final_anno.qs2   (must already
#     carry final_annotation + general_layer columns, written by 8.5)
# Writes (to ~/VisHD/8.5.normal_cell_subtypes/general_layer/):
#   <group>_srt.qs2                       (raw subset)
#   <group>_marker_dotplot.png            (curated markers by final_annotation)
#   <group>_marker_dotplot_cluster.png    (curated markers by pearson_clusters_batch)
#   <group>_recluster_srt.qs2             (cached per group)
#   <group>/<group>_deg_subcluster{,_rpca}.{Rds,csv}
#   <group>/<group>_DimPlot_{rpca,pearson}.png
#   <group>/<group>_nCount_vlnplot_{rpca,pearson}.png
#   <group>/<group>_marker_dotplot{,_rpca}.png
# Writes (to 8.4.../final_annotation_analysis/general_marker_scores/):
#   <group>_all_marker_scores.png         (all_marker purity check FeaturePlot)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
})
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")

# ── CLI arg ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 8.5.2.general_layer_analysis.R <array-index 1-5>")
arg <- as.numeric(args[1])

out_dir <- path.expand("~/VisHD/8.4.final_clear_normal_integration")
in_srt  <- file.path(out_dir, "normal_srt_final_anno.qs2")
viz_dir <- file.path(out_dir, "final_annotation_analysis")

cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
if (!"general_layer" %in% colnames(srt@meta.data))
  stop("`general_layer` column not found on ", in_srt,
       " — run 8.5.additional_analysis.R first")
cat("  ", ncol(srt), "cells\n")

general_levels <- levels(srt$general_layer)
if (arg < 1 || arg > length(general_levels))
  stop("array index ", arg, " out of range — general_layer has ",
       length(general_levels), " levels: ", paste(general_levels, collapse = ", "))
layer      <- general_levels[arg]
layer_safe <- gsub("[^A-Za-z0-9]+", "_", layer)
cat("Processing general_layer group:", layer, "\n")

# Same marker_list as 8.5.additional_analysis.R's Section 8 all_marker-purity
# check (kept in sync manually; not sourced from that script to stay standalone).
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

min_cells_subtype <- 50
options(future.globals.maxSize = 50 * 1024^3)

subtype_out_dir <- path.expand("~/VisHD/8.5.normal_cell_subtypes")

dedup_markers <- function(marker_groups) {
  seen <- character(0)
  lapply(marker_groups, function(genes) {
    genes <- genes[!genes %in% seen]
    seen <<- c(seen, genes)
    genes
  })
}

general_marker_map <- list(
  "Stromal"    = dedup_markers(c(fibro_marker, endo_marker, Stromal_subtype)),       # CAF + Smooth muscle + Endo/Pericyte
  "Myeloid"    = dedup_markers(c(macro_marker, Myeloid_subtype)),
  "Lymphoid"   = dedup_markers(c(Bcell_subtype, Tcell_subtype)),     # B/T cells + Plasma
  "Epithelial" = dedup_markers(Epithelial_subtype),                         # Epithelial + SVEC
  "Neural"     = dedup_markers(list(Neuron = Neuron_feature))
)

general_layer_dir      <- file.path(subtype_out_dir, "general_layer")
general_marker_viz_dir <- file.path(viz_dir, "general_marker_scores")
dir.create(general_layer_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(general_marker_viz_dir, showWarnings = FALSE, recursive = TRUE)

layer_path <- file.path(general_layer_dir, paste0(layer_safe, "_recluster_srt.qs2"))
layer_dir  <- file.path(general_layer_dir, layer_safe)
dir.create(layer_dir, showWarnings = FALSE, recursive = TRUE)

markers <- general_marker_map[[layer]]

# ── 1. General-layer subset ────────────────────────────────────────────────────
# Splits srt down to this one general_layer group and saves it for downstream
# reuse, alongside a curated-marker DotPlot (general_marker_map) faceted by
# the original final_annotation and by pearson_clusters_batch so the merged
# subset can be sanity-checked. Drops genes already claimed by an earlier
# sub-list so no gene is faceted under two groups at once (which duplicates
# DotPlot's feature.groups levels).
layer_srt <- subset(srt, subset = general_layer == layer)
cat(sprintf("[%s] %d cells\n", layer, ncol(layer_srt)))
qs_save(layer_srt, file.path(general_layer_dir, paste0(layer_safe, "_srt.qs2")))

if (!is.null(markers)) {
  feats <- lapply(markers, function(g) intersect(g, rownames(layer_srt)))
  feats <- feats[lengths(feats) > 0]
  if (length(feats)) {
    dp_layer <- DotPlot(layer_srt, features = feats, group.by = "final_annotation",
                      assay = "SCT") +
      ggtitle(sprintf("%s — curated markers by final_annotation", layer)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
    ggsave(file.path(general_layer_dir, sprintf("%s_marker_dotplot.png", layer_safe)),
           dp_layer, width = max(4, length(unlist(feats)) * 0.25 + 2), height = 5,
           dpi = 300, limitsize = FALSE, bg = "white")
    dp_layer <- DotPlot(layer_srt, features = feats, group.by = "pearson_clusters_batch",
                      assay = "SCT") +
      ggtitle(sprintf("%s — curated markers by pearson_clusters_batch", layer)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
    ggsave(file.path(general_layer_dir, sprintf("%s_marker_dotplot_cluster.png", layer_safe)),
           dp_layer, width = max(4, length(unlist(feats)) * 0.25 + 2), height = 5,
           dpi = 300, limitsize = FALSE, bg = "white")
  }
}
cat(sprintf("[%s] saved subset to %s\n", layer, general_layer_dir))

# ── 2. General-layer subtype processing ───────────────────────────────────────
# Subset + recluster (batch-corrected Pearson PCA) and score every all_marker
# signature, then RPCA-integrate and run MAST DE on both the RPCA and
# pearson_clusters_batch clusterings. Pure compute/cache — no plotting here;
# Section 3 draws everything from the cache.
if (ncol(layer_srt) < min_cells_subtype) {
  cat(sprintf("[%s] skipped — only %d cells (< %d)\n", layer, ncol(layer_srt), min_cells_subtype))
  quit(save = "no", status = 0)
}

rpca_ok <- FALSE

if (file.exists(layer_path)) {
  layer_srt <- qs_read(layer_path)
  rpca_ok <- "umap.rpca" %in% Reductions(layer_srt)
}

if (!all(names(all_marker) %in% colnames(layer_srt@meta.data))) {
  cat(sprintf("[%s] reclustering %d cells...\n", layer, ncol(layer_srt)))
  layer_srt <- JoinLayers(layer_srt, assay = "Spatial")
  layer_srt <- SCTransform(layer_srt, assay = "Spatial", verbose = FALSE)
  hvg_sct <- VariableFeatures(layer_srt, assay = "SCT")
  spatial_genes <- rownames(layer_srt[["Spatial"]])
  hvg <- intersect(hvg_sct, spatial_genes)
  VariableFeatures(layer_srt, assay = "Spatial") <- hvg
  layer_srt <- do.pearson_pca(layer_srt, batch_variable = "slide", assay = "Spatial",
                    find_hvgs = FALSE, reduction_prefix = "pearsonbatch",
                    clusters_col = "pearson_clusters_batch", resolution = 1)
  n_sct_genes <- nrow(layer_srt[["SCT"]])
  for (nb in c(24, 12, 6, 3, 1)) {
    scored <- tryCatch(
      AddModuleScore(layer_srt, features = marker_list, name = "Cluster",
                     assay = "SCT", nbin = nb, ctrl = min(100, n_sct_genes)),
      error = function(e) NULL
    )
    if (!is.null(scored)) { layer_srt <- scored; break }
  }
  added <- paste0("Cluster", seq_along(marker_list))
  colnames(layer_srt@meta.data)[match(added, colnames(layer_srt@meta.data))] <- names(marker_list)
  qs_save(layer_srt, layer_path)
} else {
  cat(sprintf("[%s] marker scores already present — skipping recluster\n", layer))
}

if (is.null(markers)) {
  cat(sprintf("[%s] no matching curated marker list — skipping RPCA/DE\n", layer))
  quit(save = "no", status = 0)
}

if (!rpca_ok) {
  cat(sprintf("[%s] RPCA integration...\n", layer))
  layer_srt[["Spatial"]] <- split(layer_srt[["Spatial"]], f = layer_srt$slide)
  DefaultAssay(layer_srt) <- "Spatial"
  layer_srt <- NormalizeData(layer_srt)
  layer_srt <- FindVariableFeatures(layer_srt)
  layer_srt <- ScaleData(layer_srt)
  layer_srt <- RunPCA(layer_srt)
  # IntegrateLayers can error if the actual anchor count (found internally,
  # governed by k.anchor / shared nearest neighbours) falls below k.weight
  # for a small slide; retry with progressively smaller k.weight.
  kw_max <- min(50, min(table(layer_srt$slide)))
  for (kw in unique(pmin(kw_max, c(50, 30, 20, 15, 10, 5)))) {
    integrated <- tryCatch(
      IntegrateLayers(object = layer_srt, method = RPCAIntegration, orig.reduction = "pca",
        new.reduction = "integrated.rpca", verbose = FALSE, k.weight = kw),
      error = function(e) {
        message(sprintf("  [%s] IntegrateLayers failed at k.weight=%d: %s", layer, kw, conditionMessage(e)))
        NULL
      })
    if (!is.null(integrated)) { layer_srt <- integrated; rpca_ok <- TRUE; break }
  }
  if (rpca_ok) {
    layer_srt <- FindNeighbors(layer_srt, reduction = "integrated.rpca", dims = 1:30)
    layer_srt <- FindClusters(layer_srt, resolution = 1, cluster.name = "rpca_clusters")
    layer_srt <- RunUMAP(layer_srt, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
    qs_save(layer_srt, layer_path)
  } else {
    cat(sprintf("[%s] RPCA integration failed at all k.weight values — skipping rpca DE\n", layer))
  }
}

de_rds_rpca <- file.path(layer_dir, sprintf("%s_deg_subcluster_rpca.Rds", layer_safe))
if (!rpca_ok) {
  cat(sprintf("[%s] no RPCA integration — skipping rpca DE\n", layer))
} else if (file.exists(de_rds_rpca)) {
  cat(sprintf("[%s] MAST DE (rpca) already cached\n", layer))
} else {
  Idents(layer_srt) <- "rpca_clusters"
  cat(sprintf("[%s] MAST DE between %d rpca subclusters...\n", layer, dplyr::n_distinct(Idents(layer_srt))))
  DE <- tryCatch(
    FindAllMarkers(layer_srt, assay = "SCT", test.use = "MAST",
                   only.pos = TRUE, verbose = FALSE),
    error = function(e) { message(sprintf("  FindAllMarkers failed for %s: %s", layer, conditionMessage(e))); NULL })
  if (!is.null(DE) && nrow(DE) > 0) {
    print(head(DE))
    saveRDS(DE, de_rds_rpca)
    write.csv(DE, file.path(layer_dir, sprintf("%s_deg_subcluster_rpca.csv", layer_safe)),
              row.names = FALSE)
  }
}

de_rds_pearson <- file.path(layer_dir, sprintf("%s_deg_subcluster.Rds", layer_safe))
if (dplyr::n_distinct(layer_srt$pearson_clusters_batch) < 2) {
  cat(sprintf("[%s] fewer than 2 pearson subclusters — skipping pearson DE\n", layer))
} else if (file.exists(de_rds_pearson)) {
  cat(sprintf("[%s] MAST DE (pearson) already cached\n", layer))
} else {
  Idents(layer_srt) <- "pearson_clusters_batch"
  cat(sprintf("[%s] MAST DE between %d pearson subclusters...\n", layer, dplyr::n_distinct(Idents(layer_srt))))
  DE <- tryCatch(
    FindAllMarkers(layer_srt, assay = "SCT", test.use = "MAST",
                   only.pos = TRUE, verbose = FALSE),
    error = function(e) { message(sprintf("  FindAllMarkers failed for %s: %s", layer, conditionMessage(e))); NULL })
  if (!is.null(DE) && nrow(DE) > 0) {
    print(head(DE))
    saveRDS(DE, de_rds_pearson)
    write.csv(DE, file.path(layer_dir, sprintf("%s_deg_subcluster.csv", layer_safe)),
              row.names = FALSE)
  }
}

# ── 3. General-layer subtype visualization ────────────────────────────────────
# Reads back the Section 2 cache: all-marker FeaturePlot (purity check) plus,
# since this group has a curated marker list, RPCA and pearson_clusters_batch
# DimPlot/VlnPlot/DotPlot against that group's markers.
fp <- FeaturePlot(layer_srt, names(marker_list), reduction = "pearsonbatchumap",
                  cols = c("white", "red"), order = TRUE) +
  plot_layout(ncol = 4)
ggsave(file.path(general_marker_viz_dir, paste0(layer_safe, "_all_marker_scores.png")),
       fp, width = 16, height = 12, dpi = 300, limitsize = FALSE)

if (!rpca_ok) {
  cat(sprintf("[%s] no RPCA integration cached — skipping subtype plots\n", layer))
  quit(save = "no", status = 0)
}

dp_sub <- DimPlot(layer_srt, reduction = "umap.rpca",
                  group.by = "rpca_clusters", label = TRUE,
                  repel = TRUE, cols = as.vector(polychrome())) +
  ggtitle(sprintf("%s — RPCA DimPlot", layer))
dp_sub2 <- DimPlot(layer_srt, reduction = "umap.rpca",
                  group.by = "slide", label = TRUE,
                  repel = TRUE, cols = as.vector(polychrome()))
ggsave(file.path(layer_dir, sprintf("%s_DimPlot_rpca.png", layer_safe)),
       dp_sub+dp_sub2, width = 10, height = 4, dpi = 300, bg = "white")

vln <- VlnPlot(layer_srt, features = "nCount_Spatial", group.by = "rpca_clusters") +
  ggtitle(sprintf("%s — nCount_Spatial per subcluster", layer)) +
  NoLegend()
ggsave(file.path(layer_dir, sprintf("%s_nCount_vlnplot_rpca.png", layer_safe)),
       vln, width = 6, height = 4, dpi = 300, bg = "white")

feats <- lapply(markers, function(g) intersect(g, rownames(layer_srt)))
feats <- feats[lengths(feats) > 0]
if (!length(feats)) {
  cat(sprintf("[%s] no marker genes present — skipping DotPlots\n", layer))
  quit(save = "no", status = 0)
}

Idents(layer_srt) <- "rpca_clusters"
dp <- DotPlot(layer_srt, features = feats, group.by = "rpca_clusters",
             assay = "SCT") +
  ggtitle(sprintf("%s — marker expression per RPCA subcluster", layer)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(file.path(layer_dir, sprintf("%s_marker_dotplot_rpca.png", layer_safe)),
       dp, height = 5, width = max(4, length(unlist(feats)) * 0.25 + 2),
       dpi = 300, limitsize = FALSE, bg = "white")

if (dplyr::n_distinct(layer_srt$pearson_clusters_batch) < 2) {
  cat(sprintf("[%s] fewer than 2 pearson subclusters — skipping pearson plots\n", layer))
  quit(save = "no", status = 0)
}

vln <- VlnPlot(layer_srt, features = "nCount_Spatial", group.by = "pearson_clusters_batch") +
  ggtitle(sprintf("%s — nCount_Spatial per subcluster", layer)) +
  NoLegend()
ggsave(file.path(layer_dir, sprintf("%s_nCount_vlnplot_pearson.png", layer_safe)),
       vln, width = 6, height = 4, dpi = 300, bg = "white")

dp_sub <- DimPlot(layer_srt, reduction = "pearsonbatchumap",
                  group.by = "pearson_clusters_batch", label = TRUE,
                  repel = TRUE, cols = as.vector(polychrome())) +
  ggtitle(sprintf("%s — pearson batch-corrected DimPlot", layer))
dp_sub2 <- DimPlot(layer_srt, reduction = "pearsonbatchumap",
                  group.by = "slide", label = TRUE,
                  repel = TRUE, cols = as.vector(polychrome()))
ggsave(file.path(layer_dir, sprintf("%s_DimPlot_pearson.png", layer_safe)),
       dp_sub+dp_sub2, width = 10, height = 4, dpi = 300, bg = "white")

Idents(layer_srt) <- "pearson_clusters_batch"
dp <- DotPlot(layer_srt, features = feats, group.by = "pearson_clusters_batch",
             assay = "SCT") +
  ggtitle(sprintf("%s — marker expression per subcluster", layer)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(file.path(layer_dir, sprintf("%s_marker_dotplot.png", layer_safe)),
       dp, height = 5, width = max(4, length(unlist(feats)) * 0.25 + 2),
       dpi = 300, limitsize = FALSE, bg = "white")

cat(sprintf("\n[%s] Done. Outputs in %s and %s\n", layer, general_marker_viz_dir, layer_dir))
