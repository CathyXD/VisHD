#!/usr/bin/env Rscript
# 8.5.3.align_refined_subtypes.R   (run-once, loops over all 5 general_layer groups)
#
# Joins the manually-refined subcluster identity calls in
# 8.5.normal_cell_subtypes/refined_tables/ALL_subtypes_refined.csv back onto
# the per-group Seurat objects produced by 8.5.2.general_layer_analysis.R,
# builds one combined per-cell annotation table across all 5 groups, and
# visualizes the refined subtypes on both the pearson and rpca reductions.
#
# `method` in ALL_subtypes_refined.csv: "default" -> pearson_clusters_batch /
# pearsonbatchumap reduction; "rpca" -> rpca_clusters / umap.rpca reduction.
#
# Loads:
#   8.5.normal_cell_subtypes/refined_tables/ALL_subtypes_refined.csv
#   8.5.normal_cell_subtypes/general_layer/<group>_recluster_srt.qs2   (x5)
# Writes:
#   8.5.normal_cell_subtypes/general_layer/<group>_recluster_srt.qs2   (overwritten,
#     gains refined_subtype_{pearson,rpca}/status_{pearson,rpca}/quality_flag_{pearson,rpca})
#   8.5.normal_cell_subtypes/refined_tables/integrated_refined_annotation.{csv,qs2}
#   8.5.normal_cell_subtypes/general_layer/<group>/<group>_DimPlot_refined_{pearson,rpca}.png
#   8.5.normal_cell_subtypes/general_layer/<group>/<group>_refined_subtype_markers_{pearson,rpca}.png

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
})

general_layer_dir <- path.expand("~/VisHD/8.5.normal_cell_subtypes/general_layer")
refined_dir       <- path.expand("~/VisHD/8.5.normal_cell_subtypes/refined_tables")
general_levels    <- c("Stromal", "Myeloid", "Lymphoid", "Epithelial", "Neural")

refined <- read.csv(file.path(refined_dir, "ALL_subtypes_refined.csv"), stringsAsFactors = FALSE)

# ── helpers ──────────────────────────────────────────────────────────────────
join_refined_annotation <- function(srt, layer, method, cluster_col, suffix) {
  df <- dplyr::filter(refined, cell_type == layer, method == !!method)
  if (nrow(df) == 0) return(srt)
  key         <- as.character(df$cluster)
  cluster_chr <- as.character(srt@meta.data[[cluster_col]])
  lut_subtype <- setNames(df$refined_subtype, key)
  lut_status  <- setNames(df$status, key)
  lut_quality <- setNames(df$quality_flag, key)
  srt@meta.data[[paste0("refined_subtype_", suffix)]] <- unname(lut_subtype[cluster_chr])
  srt@meta.data[[paste0("status_", suffix)]]          <- unname(lut_status[cluster_chr])
  srt@meta.data[[paste0("quality_flag_", suffix)]]    <- unname(lut_quality[cluster_chr])
  srt
}

# One representative gene per on-lineage refined_subtype label (deduped —
# shared identity, not per subcluster instance), walking the ranked
# top_markers list until a gene present in the object is found.
get_representative_genes <- function(layer, method, gene_pool) {
  df <- dplyr::filter(refined, cell_type == layer, method == !!method, status == "on-lineage")
  df <- df[!duplicated(df$refined_subtype), , drop = FALSE]
  genes <- setNames(rep(NA_character_, nrow(df)), df$refined_subtype)
  for (i in seq_len(nrow(df))) {
    candidates <- trimws(strsplit(df$top_markers[i], ";")[[1]])
    hit <- candidates[candidates %in% gene_pool]
    if (length(hit)) genes[i] <- hit[1]
  }
  genes[!is.na(genes)]
}

combined_meta <- list()

for (layer in general_levels) {
  layer_path <- file.path(general_layer_dir, paste0(layer, "_recluster_srt.qs2"))
  if (!file.exists(layer_path)) {
    cat(sprintf("[%s] %s not found — skipping\n", layer, layer_path))
    next
  }
  cat("Loading", layer_path, "\n")
  srt <- qs_read(layer_path)

  already_done <- all(c("refined_subtype_pearson", "status_pearson") %in% colnames(srt@meta.data))
  rpca_ok <- "umap.rpca" %in% Reductions(srt) && "rpca_clusters" %in% colnames(srt@meta.data)

  if (!already_done) {
    cat(sprintf("[%s] joining refined annotation (pearson)...\n", layer))
    srt <- join_refined_annotation(srt, layer, "default", "pearson_clusters_batch", "pearson")
    if (rpca_ok) {
      cat(sprintf("[%s] joining refined annotation (rpca)...\n", layer))
      srt <- join_refined_annotation(srt, layer, "rpca", "rpca_clusters", "rpca")
    }
    qs_save(srt, layer_path)
  } else {
    cat(sprintf("[%s] refined annotation already present — skipping join\n", layer))
  }

  # ── combined per-cell meta ──────────────────────────────────────────────
  combined_meta[[layer]] <- data.frame(
    barcode                  = colnames(srt),
    slide                    = srt$slide,
    general_layer            = layer,
    pearson_clusters_batch   = as.character(srt$pearson_clusters_batch),
    rpca_clusters            = if (rpca_ok) as.character(srt$rpca_clusters) else NA_character_,
    refined_subtype_pearson  = srt$refined_subtype_pearson,
    status_pearson           = srt$status_pearson,
    refined_subtype_rpca     = if ("refined_subtype_rpca" %in% colnames(srt@meta.data)) srt$refined_subtype_rpca else NA_character_,
    status_rpca              = if ("status_rpca" %in% colnames(srt@meta.data)) srt$status_rpca else NA_character_,
    stringsAsFactors = FALSE
  )

  # ── visualization ────────────────────────────────────────────────────────
  layer_dir <- file.path(general_layer_dir, layer)
  dir.create(layer_dir, showWarnings = FALSE, recursive = TRUE)

  dp_pearson <- DimPlot(srt, reduction = "pearsonbatchumap",
                        group.by = "refined_subtype_pearson", label = TRUE,
                        repel = TRUE, cols = as.vector(polychrome())) +
    ggtitle(sprintf("%s — refined subtype (pearson)", layer))
  ggsave(file.path(layer_dir, sprintf("%s_DimPlot_refined_pearson.png", layer)),
         dp_pearson, width = 8, height = 6, dpi = 300, bg = "white")

  fp_genes_pearson <- get_representative_genes(layer, "default", rownames(srt))
  if (length(fp_genes_pearson)) {
    fp <- FeaturePlot(srt, unname(fp_genes_pearson), reduction = "pearsonbatchumap",
                      cols = c("white", "red"), order = TRUE) +
      plot_layout(ncol = 4)
    ggsave(file.path(layer_dir, sprintf("%s_refined_subtype_markers_pearson.png", layer)),
           fp, width = 16, height = ceiling(length(fp_genes_pearson) / 4) * 4,
           dpi = 300, limitsize = FALSE, bg = "white")
  } else {
    cat(sprintf("[%s] no representative genes found for pearson subtypes — skipping FeaturePlot\n", layer))
  }

  if (!rpca_ok) {
    cat(sprintf("[%s] no RPCA reduction cached — skipping rpca plots\n", layer))
    next
  }

  dp_rpca <- DimPlot(srt, reduction = "umap.rpca",
                     group.by = "refined_subtype_rpca", label = TRUE,
                     repel = TRUE, cols = as.vector(polychrome())) +
    ggtitle(sprintf("%s — refined subtype (rpca)", layer))
  ggsave(file.path(layer_dir, sprintf("%s_DimPlot_refined_rpca.png", layer)),
         dp_rpca, width = 8, height = 6, dpi = 300, bg = "white")

  fp_genes_rpca <- get_representative_genes(layer, "rpca", rownames(srt))
  if (length(fp_genes_rpca)) {
    fp <- FeaturePlot(srt, unname(fp_genes_rpca), reduction = "umap.rpca",
                      cols = c("white", "red"), order = TRUE) +
      plot_layout(ncol = 4)
    ggsave(file.path(layer_dir, sprintf("%s_refined_subtype_markers_rpca.png", layer)),
           fp, width = 16, height = ceiling(length(fp_genes_rpca) / 4) * 4,
           dpi = 300, limitsize = FALSE, bg = "white")
  } else {
    cat(sprintf("[%s] no representative genes found for rpca subtypes — skipping FeaturePlot\n", layer))
  }
}

# ── combined integrated meta across all 5 groups ────────────────────────────
integrated_meta <- dplyr::bind_rows(combined_meta)
write.csv(integrated_meta, file.path(refined_dir, "integrated_refined_annotation.csv"), row.names = FALSE)
qs_save(integrated_meta, file.path(refined_dir, "integrated_refined_annotation.qs2"))
cat(sprintf("\nSaved integrated meta for %d cells across %d groups to %s\n",
            nrow(integrated_meta), length(combined_meta), refined_dir))
