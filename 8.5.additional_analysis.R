#!/usr/bin/env Rscript
# 8.5.additional_analysis.R
# Extra analysis on the final normal annotation. Mirrors the descriptive /
# visual half of 8.1.scimilarity_check.R, but keyed on `final_annotation`
# (joined from 9.4.1's final_annotation.csv by cell_ID) instead of the
# SCimilarity hint, and on the 9.4.1 re-embedded object.
#
# Loads:
#   8.4.final_clear_normal_integration/normal_srt_final_anno.qs2
#   8.4.final_clear_normal_integration/final_annotation.csv   (cell_ID -> final_annotation)
# Writes (to 9.4.1 out_dir/final_annotation_analysis/):
#   0_all_marker_dotplot.png            DotPlot of all_marker by final_annotation
#   1_final_annotation_DimPlot.png      DimPlot on pearsonbatchumap
#   2_composition_bar.png               per-slide x category composition
#   3_final_annotation_boxplot.png      per-slide proportion of each annotation
#   3b_DT_vs_CB_boxplot.png             DT vs CB per cell type (paired Wilcoxon)
#   composition.csv
#   meta_program_scores/<celltype>/<program>.png        (AddModuleScore)
#   meta_program_scores/<celltype>_UCell/<program>.png  (UCell)
#   ../normal_subtype_reclustering/<celltype>_recluster_srt.qs2  (cached per subtype)
#   subtype_marker_scores/<celltype>_all_marker_scores.png       (all_marker on each subtype's own reclustering)
# Writes (to ~/VisHD/8.5.normal_cell_subtypes/cell_subtype/):
#   <celltype>_recluster_srt.qs2               (cached reclustering, flat)
#   <celltype>/<celltype>_deg_subcluster{,_rpca}.{Rds,csv}   (MAST DE, pearson + rpca subclusters)
#   <celltype>/<celltype>_DimPlot_{rpca,pearson}.png
#   <celltype>/<celltype>_nCount_vlnplot_{rpca,pearson}.png
#   <celltype>/<celltype>_marker_dotplot{,_rpca}.png
#
# The general_layer (coarse Stromal/Myeloid/Lymphoid/Epithelial/Neural) subset
# -> recluster -> RPCA-integrate -> MAST DE -> visualize pipeline previously
# run here as Sections 7-9 now lives in 8.5.2.general_layer_analysis.R, which
# reads the `general_layer` column this script derives and saves back onto
# normal_srt_final_anno.qs2 (must run after this script at least once).
#
# Skipped vs 9.2 (annotation is already final): TME annotation pipeline, the
# DT-vs-CB paired Wilcoxon, and the annotated qs2 / AnnData re-export.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(purrr)
  library(stringr)
  library(pals)
  library(ggpubr)
  library(qs2)
  library(UCell, lib.loc = "~/R_Library/4.5")
})
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")

out_dir <- path.expand("~/VisHD/8.4.final_clear_normal_integration")
in_srt  <- file.path(out_dir, "normal_srt_final_anno.qs2")
in_csv  <- file.path(out_dir, "final_annotation.csv")
viz_dir <- file.path(out_dir, "final_annotation_analysis")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
cat("  ", ncol(srt), "cells\n")

# ── Join final_annotation from the CSV by cell_ID (mirror of 9.2's hint join) ─
anno <- read.csv(in_csv, check.names = FALSE, stringsAsFactors = FALSE)
rownames(anno) <- anno$cell_ID
cat("Loaded", nrow(anno), "annotated cells from", in_csv, "\n")

srt$final_annotation <- anno[srt$cell_ID, "final_annotation"]
srt$source_cluster   <- anno[srt$cell_ID, "source_cluster"]
na_frac <- mean(is.na(srt$final_annotation))
cat(sprintf("NA fraction after join: %.4f (%d / %d)\n",
            na_frac, sum(is.na(srt$final_annotation)), ncol(srt)))

# ── 0. All-marker DotPlot across final_annotation ─────────────────────────────
# Cells grouped by final_annotation; markers grouped by all_marker's own list
# split (Hepatocyte, Neuron, Epithelial, ... — DotPlot facets by list name).
# Marker panels grouped by broad cell-type compartment
marker_list <- list(
  # --- Stromal / mesenchymal ---
  CAF           = c("COL1A1", "COL1A2", "DCN", "LUM", "SFRP4"),
  Smooth_muscle = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES"),
  Pericyte      = c("RGS5", "NOTCH3", "MCAM", "BCAM", "PDGFRB"),   # PDGFRB canonical, NOT in DEGs

  # --- Immune ---
  Macrophages   = c("CD68", "CD163", "CSF1R", "C1QA", "MS4A7"),
  T_cells       = c("TRAC", "TRBC2", "IL7R", "IL32", "IKZF1"),
  B_cells       = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),   # canonical, NOT in DEGs
  Plasma        = c("MZB1", "JCHAIN", "XBP1", "IGHG1", "TXNDC5"),

  # --- Epithelial ---
  Epithelial    = c("KRT5", "KRT15", "KRT17", "TACSTD2", "CLDN4"),
  SVEC          = c("SEMG1", "SEMG2", "MUC6", "PIP", "PATE1"),

  # --- Neural (Schwann cells/peripheral glia) ---
  Neurons       = c("S100B", "NRXN1", "SCN7A", "PMP22", "PTGDS") 
)

srt$final_annotation[srt$final_annotation == "Pericyte"] <- "Endo/Pericyte"
srt$final_annotation[srt$final_annotation == "B cells"] <- "B/T cells"
srt$final_annotation[srt$final_annotation == "Neurons"] <- "Glial cells"

celltype_levels <- c(
  "CAF", "Smooth muscle", "Endo/Pericyte",   # stromal
  "Macrophages", "B/T cells", "Plasma",   # immune
  "Epithelial", "SVEC",                 # epithelial
  "Glial cells"                             # neural
)
srt$final_annotation <- factor(srt$final_annotation, levels = celltype_levels)

# General-layer grouping: coarse compartment reference for subsetting
general_layer_map <- c(
  "B/T cells"     = "Lymphoid",
  "Plasma"        = "Lymphoid",
  "Macrophages"   = "Myeloid",
  "Epithelial"    = "Epithelial",
  "SVEC"          = "Epithelial",
  "Smooth muscle" = "Stromal",
  "CAF"           = "Stromal",
  "Endo/Pericyte" = "Stromal",
  "Glial cells" = "Neural"
)
srt$general_layer <- factor(unname(general_layer_map[as.character(srt$final_annotation)]),
                             levels = c("Stromal", "Myeloid", "Lymphoid", "Epithelial", "Neural"))

feats_all <- lapply(marker_list, function(g) intersect(g, rownames(srt)))
feats_all <- feats_all[lengths(feats_all) > 0]
dp_all <- DotPlot(srt, features = feats_all, group.by = "final_annotation",
                  assay = "SCT") +
  # coord_flip() +
  ggtitle("marker_list expression by final annotation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(file.path(viz_dir, "0_marker_list_dotplot.png"),
       dp_all, width = 20, height = 6,
       dpi = 300, limitsize = FALSE,bg = "white")


# ── 1. DimPlot of the final annotation on the batch-corrected UMAP ────────────
dp <- DimPlot(srt, reduction = "pearsonbatchumap",
              group.by = "final_annotation", label = TRUE, label.size = 6,
              repel = TRUE, cols = as.vector(polychrome())) +
  ggtitle("Final normal annotation (batch corrected)") +
  theme(legend.text = element_text(size = 15))
ggsave(file.path(viz_dir, "1_final_annotation_DimPlot.png"),
       dp, width = 9, height = 7, dpi = 400)

# ── 2. Per-slide x category composition (stacked bar faceted by slide) ────────
has_category <- "category" %in% colnames(srt@meta.data)
if (has_category) {
  comp <- srt@meta.data %>%
    mutate(cat2 = ifelse(startsWith(as.character(category), "CB"), "CB", "DT")) %>%
    count(slide, cat2, final_annotation) %>%
    group_by(slide, cat2) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(slide, cat2, final_annotation, fill = list(n = 0, prop = 0))
  bar <- ggplot(comp, aes(x = cat2, y = prop, fill = final_annotation)) +
    geom_col() +
    facet_grid(~ slide) +
    labs(x = "Category", y = "Proportion",
         title = "Final annotation composition per slide x category")
} else {
  comp <- srt@meta.data %>%
    count(slide, final_annotation) %>%
    group_by(slide) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(slide, final_annotation, fill = list(n = 0, prop = 0))
  bar <- ggplot(comp, aes(x = slide, y = prop, fill = final_annotation)) +
    geom_col() +
    labs(x = "Slide", y = "Proportion",
         title = "Final annotation composition per slide")
}
bar <- bar +
  scale_fill_manual(values = as.vector(polychrome())) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 9))
ggsave(file.path(viz_dir, "2_composition_bar.png"), bar, width = 8, height = 4, dpi = 400)
write.csv(comp, file.path(viz_dir, "composition.csv"), row.names = FALSE)

bar2 <- bar + facet_grid(final_annotation ~ slide, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 6))+
  theme_bw(base_size = 15) 
ggsave(file.path(viz_dir, "2_composition_bar_facet.png"), bar2, width = 11, height = 11, dpi = 400)

# ── 2c. DT - CB proportion difference per cell type x slide ───────────────────
if (has_category) {
  diff_df <- comp %>%
    select(slide, cat2, final_annotation, prop) %>%
    pivot_wider(names_from = cat2, values_from = prop) %>%
    mutate(diff = DT - CB)

  diff_bar <- ggplot(diff_df, aes(x = "", y = diff, fill = diff > 0)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "blue"), guide = "none") +
    facet_grid(final_annotation ~ slide, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Proportion difference (DT - CB)",
         title = "Final annotation — DT minus CB proportion per cell type x slide") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 6))
  ggsave(file.path(viz_dir, "2c_DT_minus_CB_diff_bar.png"),
         diff_bar, width = 8, height = 8, dpi = 400)
}

# ── 3. Boxplot: each annotation's proportion across the 8 slides ──────────────
comp_slide <- srt@meta.data %>%
  count(slide, final_annotation) %>%
  group_by(slide) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(slide, final_annotation, fill = list(n = 0, prop = 0))

box <- ggplot(comp_slide,
              aes(x = reorder(final_annotation, -prop, FUN = median), y = prop)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(aes(color = slide), width = 0.15, alpha = 0.6, size = 1) +
  labs(x = NULL, y = "Proportion (per slide)",
       title = "Per-slide variability of each final annotation") +
  theme_minimal(base_size = 12) +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(viz_dir, "3_final_annotation_boxplot.png"),
       box, width = 10, height = 2, dpi = 400)

# ── 3b. Paired DT vs CB boxplot per cell type (slide is the pairing unit) ─────
# Uses the per-slide x cat2 proportions from section 2 (CB = CB 0 + CB 1).
if (has_category) {
  comp_paired <- comp %>%
    group_by(final_annotation, slide) %>%
    filter(n() == 2) %>%                       # keep slides with both DT and CB
    ungroup() %>%
    mutate(group = factor(cat2, levels = c("DT", "CB")))

  box_pair <- ggpaired(comp_paired,
                       x = "group", y = "prop", id = "slide",
                       facet.by = "final_annotation",
                       fill = "group", palette = c("#E69F00", "#56B4E9"),
                       line.color = "grey60", line.size = 0.3,
                       point.size = 1, short.panel.labs = TRUE) +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       label = "p.format", size = 2.6, label.y.npc = 0.7) +
    labs(x = NULL, y = "Proportion (per slide x category)",
         title = "Final annotation — DT vs CB per cell type (paired Wilcoxon)") +
    theme(axis.text.x = element_text(size = 8),
          strip.text = element_text(size = 7))
  ggsave(file.path(viz_dir, "3b_DT_vs_CB_boxplot.png"),
         box_pair, width = 6, height = 5, dpi = 400)
}

# ── 4. Meta-program module scores (FeaturePlot per program, per cell-type) ────
meta_xlsx <- path.expand("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx")
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
meta_programs <- set_names(sheetname, sheetname) |>
  purrr::map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  purrr::map(~ purrr::map(.x, ~ as.character(na.omit(.x))))

mp_root <- file.path(viz_dir, "meta_program_scores")
dir.create(mp_root, showWarnings = FALSE, recursive = TRUE)
genes_in_srt <- rownames(srt)

# AddModuleScore (one program at a time, blue gradient)
for (ct in names(meta_programs)) {
  ct_dir <- file.path(mp_root, gsub("[^A-Za-z0-9]+", "_", ct))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  for (prog in names(meta_programs[[ct]])) {
    genes <- intersect(meta_programs[[ct]][[prog]], genes_in_srt)
    if (length(genes) < 3) next
    srt <- AddModuleScore(srt, features = list(genes), name = "mp_score", assay = "SCT")
    fp <- FeaturePlot(srt, features = "mp_score1",
                      reduction = "pearsonbatchumap", order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, NA)) +
      ggtitle(sprintf("%s - %s (%d genes)", ct, prog, length(genes)))
    ggsave(file.path(ct_dir, sprintf("%s.png", gsub("[^A-Za-z0-9]+", "_", prog))),
           fp, width = 6, height = 5, dpi = 300)
    srt[["mp_score1"]] <- NULL
  }
}

# UCell (whole sheet at once, red gradient)
for (ct in names(meta_programs)) {
  ct_dir <- file.path(mp_root, paste(gsub("[^A-Za-z0-9]+", "_", ct), "UCell", sep = "_"))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  srt <- AddModuleScore_UCell(srt, features = meta_programs[[ct]], missing_genes = "skip")
  for (prog in names(meta_programs[[ct]])) {
    genes <- intersect(meta_programs[[ct]][[prog]], genes_in_srt)
    fp <- FeaturePlot(srt, features = paste0(prog, "_UCell"),
                      reduction = "pearsonbatchumap", order = TRUE) +
      scale_color_gradient(low = "snow", high = "red", limits = c(0, NA)) +
      ggtitle(sprintf("%s - %s (%d genes)", ct, prog, length(genes)))
    ggsave(file.path(ct_dir, sprintf("%s.png", gsub("[^A-Za-z0-9]+", "_", prog))),
           fp, width = 6, height = 5, dpi = 300)
  }
}
# ── 5. Per-subtype processing ─────────────────────────────────────────────────
# For each final_annotation cell type: subset + recluster (batch-corrected
# Pearson PCA) and score every all_marker signature, then RPCA-integrate and
# run MAST DE on both the RPCA and pearson_clusters_batch clusterings. Pure
# compute/cache — no plotting here; Section 6 draws everything from the cache.
marker_viz_dir  <- file.path(viz_dir, "subtype_marker_scores")
subtype_out_dir <- path.expand("~/VisHD/8.5.normal_cell_subtypes")
subtype_dir     <- file.path(subtype_out_dir, "cell_subtype")
dir.create(subtype_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(marker_viz_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(subtype_out_dir, showWarnings = FALSE, recursive = TRUE)
options(future.globals.maxSize = 50 * 1024^3)

min_cells_subtype <- 50
ct_marker_map <- list(
  "B/T cells"     = c(bcell_marker, Tcell_subtype),
  "Plasma"        = bcell_marker,
  "Macrophages"   = macro_marker,
  "Epithelial"    = epi_marker,
  "SVEC"          = epi_marker,        # Seminal Vesicle Epithelial Cell
  "CAF"           = fibro_marker,
  "Smooth muscle" = fibro_marker,
  "Endo/Pericyte" = endo_marker,
  "Glial cells" = list(Neuron = Neuron_feature)
)


celltypes <- sort(unique(na.omit(srt$final_annotation)))

cat("Saving updated final_annotation back to", in_srt, "\n")
qs_save(srt, in_srt)

for (ct in celltypes) {
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  ct_path <- file.path(subtype_dir, paste0(ct_safe, "_recluster_srt.qs2"))

  if (file.exists(ct_path)) {
    sub_srt <- qs_read(ct_path)
  } else {
    sub_srt <- subset(srt, subset = final_annotation == ct)
  }

  if (ncol(sub_srt) < min_cells_subtype) {
    cat(sprintf("[%s] skipped — only %d cells (< %d)\n", ct, ncol(sub_srt), min_cells_subtype))
    next
  }

  if (!all(names(all_marker) %in% colnames(sub_srt@meta.data))) {
    cat(sprintf("[%s] reclustering %d cells...\n", ct, ncol(sub_srt)))
    sub_srt <- SCTransform(sub_srt, assay = "Spatial", verbose = FALSE)
    hvg_sct <- VariableFeatures(sub_srt, assay = "SCT")
    spatial_genes <- rownames(GetAssayData(sub_srt, assay = "Spatial", layer = "counts"))
    hvg <- intersect(hvg_sct, spatial_genes)
    VariableFeatures(sub_srt, assay = "Spatial") <- hvg
    sub_srt <- do.pearson_pca(sub_srt,batch_variable = "slide", assay = "Spatial",
                      find_hvgs = FALSE, reduction_prefix = "pearsonbatch",
                      clusters_col = "pearson_clusters_batch", resolution = 1)
    n_sct_genes <- nrow(GetAssayData(sub_srt, assay = "SCT"))
    for (nb in c(24, 12, 6, 3, 1)) {
      scored <- tryCatch(
        AddModuleScore(sub_srt, features = marker_list, name = "Cluster",
                       assay = "SCT", nbin = nb, ctrl = min(100, n_sct_genes)),
        error = function(e) NULL
      )
      if (!is.null(scored)) { sub_srt <- scored; break }
    }
    added <- paste0("Cluster", seq_along(marker_list))
    colnames(sub_srt@meta.data)[match(added, colnames(sub_srt@meta.data))] <- names(marker_list)
    qs_save(sub_srt, ct_path)
  } else {
    cat(sprintf("[%s] marker scores already present — skipping recluster\n", ct))
  }

  markers <- ct_marker_map[[ct]]
  if (is.null(markers)) {
    cat(sprintf("[%s] no matching curated marker list — skipping RPCA/DE\n", ct))
    next
  }

  rpca_ok <- "umap.rpca" %in% Reductions(sub_srt)
  if (!rpca_ok) {
    cat(sprintf("[%s] RPCA integration...\n", ct))
    sub_srt[["Spatial"]] <- split(sub_srt[["Spatial"]], f = sub_srt$slide)
    DefaultAssay(sub_srt) <- "Spatial"
    sub_srt <- NormalizeData(sub_srt)
    sub_srt <- FindVariableFeatures(sub_srt)
    sub_srt <- ScaleData(sub_srt)
    sub_srt <- RunPCA(sub_srt)
    # IntegrateLayers can error if the actual anchor count (found internally,
    # governed by k.anchor / shared nearest neighbours) falls below k.weight
    # for a small slide; retry with progressively smaller k.weight.
    kw_max <- min(50, min(table(sub_srt$slide)))
    for (kw in unique(pmin(kw_max, c(50, 30, 20, 15, 10, 5)))) {
      integrated <- tryCatch(
        IntegrateLayers(object = sub_srt, method = RPCAIntegration, orig.reduction = "pca",
          new.reduction = "integrated.rpca", verbose = FALSE, k.weight = kw),
        error = function(e) {
          message(sprintf("  [%s] IntegrateLayers failed at k.weight=%d: %s", ct, kw, conditionMessage(e)))
          NULL
        })
      if (!is.null(integrated)) { sub_srt <- integrated; rpca_ok <- TRUE; break }
    }
    if (rpca_ok) {
      sub_srt <- FindNeighbors(sub_srt, reduction = "integrated.rpca", dims = 1:30)
      sub_srt <- FindClusters(sub_srt, resolution = 1, cluster.name = "rpca_clusters")
      sub_srt <- RunUMAP(sub_srt, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
      qs_save(sub_srt, ct_path)
    } else {
      cat(sprintf("[%s] RPCA integration failed at all k.weight values — skipping rpca DE\n", ct))
    }
  }

  ct_dir <- file.path(subtype_dir, ct_safe)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  de_rds_rpca <- file.path(ct_dir, sprintf("%s_deg_subcluster_rpca.Rds", ct_safe))
  if (!rpca_ok) {
    cat(sprintf("[%s] no RPCA integration — skipping rpca DE\n", ct))
  } else if (file.exists(de_rds_rpca)) {
    cat(sprintf("[%s] MAST DE (rpca) already cached\n", ct))
  } else {
    Idents(sub_srt) <- "rpca_clusters"
    cat(sprintf("[%s] MAST DE between %d rpca subclusters...\n", ct, dplyr::n_distinct(Idents(sub_srt))))
    DE <- tryCatch(
      FindAllMarkers(sub_srt, assay = "SCT", test.use = "MAST",
                     only.pos = TRUE, verbose = FALSE),
      error = function(e) { message(sprintf("  FindAllMarkers failed for %s: %s", ct, conditionMessage(e))); NULL })
    if (!is.null(DE) && nrow(DE) > 0) {
      print(head(DE))
      saveRDS(DE, de_rds_rpca)
      write.csv(DE, file.path(ct_dir, sprintf("%s_deg_subcluster_rpca.csv", ct_safe)),
                row.names = FALSE)
    }
  }

  if (dplyr::n_distinct(sub_srt$pearson_clusters_batch) < 2) {
    cat(sprintf("[%s] fewer than 2 pearson subclusters — skipping pearson DE\n", ct))
    next
  }

  de_rds_pearson <- file.path(ct_dir, sprintf("%s_deg_subcluster.Rds", ct_safe))
  if (file.exists(de_rds_pearson)) {
    cat(sprintf("[%s] MAST DE (pearson) already cached\n", ct))
  } else {
    Idents(sub_srt) <- "pearson_clusters_batch"
    cat(sprintf("[%s] MAST DE between %d pearson subclusters...\n", ct, dplyr::n_distinct(Idents(sub_srt))))
    DE <- tryCatch(
      FindAllMarkers(sub_srt, assay = "SCT", test.use = "MAST",
                     only.pos = TRUE, verbose = FALSE),
      error = function(e) { message(sprintf("  FindAllMarkers failed for %s: %s", ct, conditionMessage(e))); NULL })
    if (!is.null(DE) && nrow(DE) > 0) {
      print(head(DE))
      saveRDS(DE, de_rds_pearson)
      write.csv(DE, file.path(ct_dir, sprintf("%s_deg_subcluster.csv", ct_safe)),
                row.names = FALSE)
    }
  }
}

# ── 6. Per-subtype visualization ──────────────────────────────────────────────
# Reads back the cache Section 5 built: all-marker FeaturePlot (purity check)
# plus, for cell types with a curated marker list, RPCA and pearson_clusters_batch
# DimPlot/VlnPlot/DotPlot against that cell type's markers.
for (ct in celltypes) {
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  ct_path <- file.path(subtype_dir, paste0(ct_safe, "_recluster_srt.qs2"))
  if (!file.exists(ct_path)) {
    cat(sprintf("[%s] skipped — no cached reclustering (too few cells)\n", ct))
    next
  }
  sub_srt <- qs_read(ct_path)

  fp <- FeaturePlot(sub_srt, names(marker_list), reduction = "pearsonbatchumap",
                    cols = c("white", "red"), order = TRUE) +
    plot_layout(ncol = 4)
  ggsave(file.path(marker_viz_dir, paste0(ct_safe, "_all_marker_scores.png")),
         fp, width = 16, height = 12, dpi = 300, limitsize = FALSE)

  markers <- ct_marker_map[[ct]]
  if (is.null(markers) || !"umap.rpca" %in% Reductions(sub_srt)) {
    cat(sprintf("[%s] no RPCA integration cached — skipping subtype plots\n", ct))
    next
  }
  ct_dir <- file.path(subtype_dir, ct_safe)

  dp_sub <- DimPlot(sub_srt, reduction = "umap.rpca",
                    group.by = "rpca_clusters", label = TRUE,
                    repel = TRUE, cols = as.vector(polychrome())) +
    ggtitle(sprintf("%s — RPCA DimPlot", ct))
  dp_sub2 <- DimPlot(sub_srt, reduction = "umap.rpca",
                    group.by = "slide", label = TRUE,
                    repel = TRUE, cols = as.vector(polychrome()))
  ggsave(file.path(ct_dir, sprintf("%s_DimPlot_rpca.png", ct_safe)),
         dp_sub+dp_sub2, width = 10, height = 4, dpi = 300, bg = "white")

  vln <- VlnPlot(sub_srt, features = "nCount_Spatial", group.by = "rpca_clusters") +
    ggtitle(sprintf("%s — nCount_Spatial per subcluster", ct)) +
    NoLegend()
  ggsave(file.path(ct_dir, sprintf("%s_nCount_vlnplot_rpca.png", ct_safe)),
         vln, width = 6, height = 4, dpi = 300, bg = "white")

  feats <- lapply(markers, function(g) intersect(g, rownames(sub_srt)))
  feats <- feats[lengths(feats) > 0]
  if (!length(feats)) {
    cat(sprintf("[%s] no marker genes present — skipping DotPlots\n", ct))
    next
  }

  Idents(sub_srt) <- "rpca_clusters"
  dp <- DotPlot(sub_srt, features = feats, group.by = "rpca_clusters",
               assay = "SCT") +
    ggtitle(sprintf("%s — marker expression per RPCA subcluster", ct)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  ggsave(file.path(ct_dir, sprintf("%s_marker_dotplot_rpca.png", ct_safe)),
         dp, height = 5, width = max(4, length(unlist(feats)) * 0.25 + 2),
         dpi = 300, limitsize = FALSE, bg = "white")

  if (dplyr::n_distinct(sub_srt$pearson_clusters_batch) < 2) {
    cat(sprintf("[%s] fewer than 2 pearson subclusters — skipping pearson plots\n", ct))
    next
  }

  vln <- VlnPlot(sub_srt, features = "nCount_Spatial", group.by = "pearson_clusters_batch") +
    ggtitle(sprintf("%s — nCount_Spatial per subcluster", ct)) +
    NoLegend()
  ggsave(file.path(ct_dir, sprintf("%s_nCount_vlnplot_pearson.png", ct_safe)),
         vln, width = 6, height = 4, dpi = 300, bg = "white")

  dp_sub <- DimPlot(sub_srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE,
                    repel = TRUE, cols = as.vector(polychrome())) +
    ggtitle(sprintf("%s — pearson batch-corrected DimPlot", ct))
  dp_sub2 <- DimPlot(sub_srt, reduction = "pearsonbatchumap",
                    group.by = "slide", label = TRUE,
                    repel = TRUE, cols = as.vector(polychrome()))
  ggsave(file.path(ct_dir, sprintf("%s_DimPlot_pearson.png", ct_safe)),
         dp_sub+dp_sub2, width = 10, height = 4, dpi = 300, bg = "white")

  Idents(sub_srt) <- "pearson_clusters_batch"
  dp <- DotPlot(sub_srt, features = feats, group.by = "pearson_clusters_batch",
               assay = "SCT") +
    ggtitle(sprintf("%s — marker expression per subcluster", ct)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  ggsave(file.path(ct_dir, sprintf("%s_marker_dotplot.png", ct_safe)),
         dp, height = 5, width = max(4, length(unlist(feats)) * 0.25 + 2),
         dpi = 300, limitsize = FALSE, bg = "white")
}

cat("\nDone. Outputs in", viz_dir, "and", subtype_out_dir, "\n")
cat("Run 8.5.2.general_layer_analysis.R next for the general_layer pipeline.\n")
