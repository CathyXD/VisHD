#!/usr/bin/env Rscript
# 8.3.additional_visual.R
# Extra visualisations built from 8.2.DEG_cluster_annotation.R outputs.
#
#  (A) Per-cluster bar plots of cluster-DEG overlap with the meta-programs:
#      one plot per cluster, every overlapping program ranked most -> least,
#      bars coloured by broad cell type (the meta-program first layer / sheet).
#  (B) The same overlap check + per-cluster bar plots against the
#      normal_markers.R `all_marker` normal-cell-type signatures.
#
# Reads (from 9.3 out_dir):
#   cluster_metaprogram_overlap.csv   (cluster, first_layer, program, n_overlap, ...)
#   significant_DEGs.csv              (gene, cluster, ...)
# Writes:
#   png/per_cluster/metaprogram_cluster_<k>.png
#   png/per_cluster/normalmarker_cluster_<k>.png
#   cluster_normalmarker_overlap.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(pals)
})
source("~/VisHD/normal_markers.R")   # provides `all_marker` (named list of signatures)

out_dir        <- path.expand("~/VisHD/8.2.DEG_cluster_annotation")
percluster_dir <- file.path(out_dir, "png", "per_cluster")
dir.create(percluster_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: one cluster's overlaps as a ranked bar (most -> least) ────────────
# Only programs/signatures with >0 overlap are drawn; x ordered by overlap desc.
ranked_bar <- function(d, label_col, fill_col, fill_pal, fill_name,
                       title, xlab, file) {
  d <- d[d$n_overlap > 0, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  d <- d[order(-d$n_overlap), , drop = FALSE]
  d$.pos <- factor(sprintf("%03d", seq_len(nrow(d))))   # ranked, lexical-safe order
  p <- ggplot(d, aes(x = .pos, y = n_overlap, fill = .data[[fill_col]])) +
    geom_col() +
    scale_x_discrete(labels = d[[label_col]]) +
    scale_fill_manual(values = fill_pal, name = fill_name) +
    labs(title = title, x = xlab, y = "# overlapping DEGs") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.key.size = unit(0.3, "cm"))
  ggsave(file, p, dpi = 200, limitsize = FALSE,
         width = max(6, 0.16 * nrow(d) + 3), height = 6)
}

# ── A. Meta-program overlap, one plot per cluster ─────────────────────────────
overlap_df <- read.csv(file.path(out_dir, "cluster_metaprogram_overlap.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
overlap_df$cluster <- as.character(overlap_df$cluster)

# Stable broad-cell-type palette shared across all per-cluster plots.
layers     <- sort(unique(overlap_df$first_layer))
layer_cols <- setNames(as.vector(polychrome())[seq_along(layers)], layers)

for (cl in sort(unique(overlap_df$cluster))) {
  ranked_bar(overlap_df[overlap_df$cluster == cl, ],
             label_col = "program", fill_col = "first_layer",
             fill_pal = layer_cols, fill_name = "Broad cell type",
             title = sprintf("Cluster %s - DEG overlap with meta-programs (ranked)", cl),
             xlab  = "Meta-program (most -> least overlap)",
             file  = file.path(percluster_dir, sprintf("metaprogram_cluster_%s.png", cl)))
}
cat("A. Meta-program per-cluster bars:", length(unique(overlap_df$cluster)), "clusters\n")

# ── B. normal_markers all_marker overlap, one plot per cluster ────────────────
sig_deg <- read.csv(file.path(out_dir, "significant_DEGs.csv"), stringsAsFactors = FALSE)
sig_deg$cluster <- as.character(sig_deg$cluster)
deg_by_cluster  <- split(sig_deg$gene, sig_deg$cluster)

marker_overlap <- map_dfr(names(deg_by_cluster), function(cl) {
  degs <- unique(deg_by_cluster[[cl]])
  map_dfr(names(all_marker), function(ct) {
    mg <- unique(all_marker[[ct]])
    ov <- intersect(degs, mg)
    data.frame(cluster = cl, celltype = ct,
               n_deg = length(degs), n_marker_genes = length(mg),
               n_overlap = length(ov),
               overlap_genes = paste(ov, collapse = ";"),
               stringsAsFactors = FALSE)
  })
})
write.csv(marker_overlap, file.path(out_dir, "cluster_normalmarker_overlap.csv"),
          row.names = FALSE)

celltypes <- sort(unique(marker_overlap$celltype))
ct_cols   <- setNames(as.vector(polychrome())[seq_along(celltypes)], celltypes)

for (cl in sort(unique(marker_overlap$cluster))) {
  ranked_bar(marker_overlap[marker_overlap$cluster == cl, ],
             label_col = "celltype", fill_col = "celltype",
             fill_pal = ct_cols, fill_name = "Cell type",
             title = sprintf("Cluster %s - DEG overlap with normal_markers (ranked)", cl),
             xlab  = "Normal-cell signature (most -> least overlap)",
             file  = file.path(percluster_dir, sprintf("normalmarker_cluster_%s.png", cl)))
}
cat("B. normal_markers per-cluster bars:", length(unique(marker_overlap$cluster)), "clusters\n")

cat("\nDone. Per-cluster plots in", percluster_dir, "\n")
