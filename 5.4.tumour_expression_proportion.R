#!/usr/bin/env Rscript
# 5.4.tumour_expression_proportion.R
# Per-gene proportion of cells with positive expression (counts > 0):
# DT vs CB across all tumour samples.
# Output: bar plots (positive proportion per sample x category) and
# DT-vs-CB boxplots with paired Wilcoxon p-values.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
source("~/VisHD/functions.R")

# ── Config ────────────────────────────────────────────────────────────────────
samples  <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
              "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")
base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "5.4.expression_proportion")
png_dir  <- file.path(outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

# ── Gene set: AR + FOLH1 + all DEGs (padj<0.05, |lfc|>1.25) from 5.2 ─────────
deg_path <- file.path(base_dir, "5.3.DESeq2_results", "paired_sample", "deseq2_res_DT_vs_CB.Rds")
top_degs <- character(0)
if (file.exists(deg_path)) {
  res_df <- readRDS(deg_path) %>% filter(!is.na(pvalue))
  top_degs <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1.25) %>%
    pull(gene)
  cat("Loaded", length(top_degs), "DEGs (padj<0.05, |lfc|>1.25) from", deg_path, "\n")
} else {
  message("DEG file not found - proceeding with AR + FOLH1 only.")
}
# Curated marker panels for individual per-gene comparison (SLC18A1 = VMAT1)
ar_genes <- c("AR", "FOLH1", "KLK3", "NKX3-1", "TMPRSS2", "KLK4", "STEAP2", "STEAP1")
ne_genes <- c("CHGA", "CHGB", "SCG2", "SLC18A1", "SYNGR4", "NPB", "PTPN5", "HDAC9", "SYP", "ENO2", "NCAM1", "ELAVL4", "GABRG2", "GABRA1", "GABRB2")

genes <- unique(c("AR", "FOLH1", top_degs, ar_genes, ne_genes))
cat("Testing", length(genes), "genes total\n")

# ── Spatial dark_feature_plots (per gene, per sample) + SpaNorm expression ───
dark_dir <- file.path(png_dir, "dark_feature_plots")
dir.create(dark_dir, showWarnings = FALSE, recursive = TRUE)
dark_plots   <- setNames(vector("list", length(genes)), genes)
spanorm_list <- list()

# Tissue coordinates for a per-sample object: FOV centroids, falling back to
# x_centroid/y_centroid meta.data (mirrors dark_feature_plot()'s own fallback).
.get_srt_coords <- function(srt) {
  coords <- tryCatch(Seurat::GetTissueCoordinates(srt, which = "centroids"),
                      error = function(e) NULL)
  if (!is.null(coords)) {
    id_col <- intersect(c("cell", "barcode", "id"), names(coords))
    idx <- match(colnames(srt), coords[[id_col[1]]])
    return(data.frame(x = coords$x[idx], y = coords$y[idx]))
  }
  if (all(c("x_centroid", "y_centroid") %in% colnames(srt@meta.data)))
    return(data.frame(x = srt$x_centroid, y = srt$y_centroid))
  stop(".get_srt_coords: no FOV and no x_centroid/y_centroid meta.data.")
}

# ── Per-sample proportion calculation ─────────────────────────────────────────
prop_list <- list()
for (s in samples) {
  f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(f)) { message("Skipping ", s, " - missing tumour_srt.qs2"); next }
  message("Loading ", s)
  obj <- qs_read(f)
  obj$category_bin <- ifelse(obj$category == "DT", "DT",
                       ifelse(grepl("^CB", obj$category), "CB", NA_character_))
  DefaultAssay(obj) <- "Spatial"
  available <- intersect(genes, rownames(obj))
  if (length(available) == 0) next

  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")[available, , drop = FALSE]

  # SpaNorm expression + spatial coords for this sample's available genes
  spanorm_mat <- GetAssayData(obj, assay = "SpaNorm", layer = "data")[available, , drop = FALSE]
  xy <- .get_srt_coords(obj)
  spanorm_list[[length(spanorm_list) + 1]] <- data.frame(
    sample     = s,
    cell       = rep(colnames(obj), length(available)),
    gene       = rep(available, each = ncol(obj)),
    expression = as.vector(as.matrix(t(spanorm_mat))),
    category   = rep(obj$category_bin, length(available)),
    pos        = as.vector(as.matrix(t(counts))) > 0,
    x          = rep(xy$x, length(available)),
    y          = rep(xy$y, length(available)),
    stringsAsFactors = FALSE
  )

  # Per-gene dark_feature_plot for this sample, titled with the sample name
  # (samples are combined side-by-side per gene once the loop finishes)
  for (g in available) {
    dark_plots[[g]][[s]] <- dark_feature_plot(obj, features = g, pt_size_range = c(0.01, 1.2),
                                               assay = "SpaNorm", layer = "data") &
      ggtitle(s)
  }

  for (cat_v in c("DT", "CB")) {
    cells <- colnames(obj)[!is.na(obj$category_bin) & obj$category_bin == cat_v]
    if (length(cells) == 0) next
    pos_count <- Matrix::rowSums(counts[, cells, drop = FALSE] > 0)
    prop_list[[length(prop_list) + 1]] <- data.frame(
      sample   = s,
      gene     = available,
      category = cat_v,
      n_cells  = length(cells),
      n_pos    = pos_count[available],
      prop     = pos_count[available] / length(cells),
      stringsAsFactors = FALSE
    )
  }
}
prop_df <- bind_rows(prop_list)
prop_df$gene <- factor(prop_df$gene, levels = genes[genes %in% prop_df$gene])
saveRDS(prop_df,  file.path(outdir, "expression_proportion.Rds"))
write.csv(prop_df, file.path(outdir, "expression_proportion.csv"), row.names = FALSE)

# ── Aggregated per-gene spatial dark_feature_plots (samples 1-8 side by side) ─
for (g in genes) {
  plist <- dark_plots[[g]]
  if (length(plist) == 0) next
  n_col <- min(4, length(plist))
  n_row <- ceiling(length(plist) / n_col)
  agg <- patchwork::wrap_plots(plist, ncol = n_col) +
    patchwork::plot_annotation(
      title = g,
      theme = theme(plot.background = element_rect(fill = "black", colour = "black"),
                    plot.title = element_text(colour = "white", hjust = 0.5, face = "bold"))
    ) &
    theme(plot.background = element_rect(fill = "black", colour = "black"))
  ggsave(file.path(dark_dir, paste0(g, "_spatial_dark.png")),
         agg, width = n_col * 3, height = n_row * 3, dpi = 200,
         bg = "black", limitsize = FALSE)
}
message("Saved aggregated dark_feature_plots to: ", dark_dir)

# ── SpaNorm expression per cell + spatial coordinates (all genes, all samples) ─
spanorm_df <- bind_rows(spanorm_list)
saveRDS(spanorm_df, file.path(outdir, "spanorm_expression_per_cell.Rds"))
write.csv(spanorm_df, file.path(outdir, "spanorm_expression_per_cell.csv"), row.names = FALSE)
cat("Saved SpaNorm per-cell expression + coordinates:", nrow(spanorm_df), "rows\n")

# ── Paired Wilcoxon test per gene (DT vs CB across samples) ──────────────────
wilcox_df <- prop_df %>%
  select(sample, gene, category, prop) %>%
  pivot_wider(names_from = category, values_from = prop) %>%
  group_by(gene) %>%
  summarise(
    n_pairs = sum(!is.na(DT) & !is.na(CB)),
    p = if (sum(!is.na(DT) & !is.na(CB)) >= 3) {
      tryCatch(wilcox.test(DT, CB, paired = TRUE)$p.value,
               error = function(e) NA_real_)
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(p_lab = ifelse(is.na(p), "p = NA", sprintf("p = %.3g", p)))
write.csv(wilcox_df, file.path(outdir, "wilcoxon_DT_vs_CB.csv"), row.names = FALSE)

# ── Plots: AR + FOLH1 and other DEGs plotted in separate figures ─────────────
# Bar plot (per-sample positive proportion) + box plot (DT vs CB, paired
# Wilcoxon p) for a subset of genes, saved with a group tag in the filename.
make_prop_plots <- function(df, wilcox_df, tag, plots = c("bar", "box")) {
  n_genes <- length(unique(df$gene))
  if (n_genes == 0) { message("No genes for ", tag, " - skipping plots"); return(invisible(NULL)) }
  n_col <- min(4, n_genes)
  n_row <- ceiling(n_genes / n_col)

  if ("bar" %in% plots) {
    p_bar <- ggplot(df, aes(sample, prop, fill = category)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      facet_wrap(~ gene, scales = "free_y", ncol = n_col) +
      scale_fill_manual(values = c(CB = "steelblue", DT = "firebrick")) +
      labs(title = paste0("Positive expression proportion - DT vs CB (", tag, ")"),
           x = NULL, y = "Positive proportion") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text  = element_text(size = 9, face = "bold"))
    ggsave(file.path(png_dir, paste0("1_barplot_positive_proportion_", tag, ".png")),
           p_bar, width = n_col * 3.5, height = n_row * 2.5, dpi = 200, limitsize = FALSE)
  }

  if ("box" %in% plots) {
    ann_y  <- df %>% group_by(gene) %>%
      summarise(y = max(prop, na.rm = TRUE) * 1.08, .groups = "drop")
    ann_df <- left_join(filter(wilcox_df, gene %in% unique(df$gene)), ann_y, by = "gene")
    p_box <- ggplot(df, aes(category, prop, fill = category)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
      geom_text(data = ann_df, aes(x = 1.5, y = y, label = p_lab),
                inherit.aes = FALSE, size = 2.8) +
      facet_wrap(~ gene, scales = "free_y", ncol = n_col) +
      scale_fill_manual(values = c(CB = "steelblue", DT = "firebrick")) +
      labs(title = paste0("DT vs CB positive expression proportion - paired Wilcoxon (", tag, ")"),
           x = NULL, y = "Positive proportion") +
      theme_classic() +
      theme(legend.position = "none",
            strip.text = element_text(size = 9, face = "bold"))
    ggsave(file.path(png_dir, paste0("2_boxplot_DT_vs_CB_", tag, ".png")),
           p_box, width = n_col * 2, height = n_row * 2.8, dpi = 200, limitsize = FALSE)
  }
}

key_genes <- c("AR", "FOLH1")
make_prop_plots(filter(prop_df,  gene %in% key_genes), wilcox_df, "AR_FOLH1")
make_prop_plots(filter(prop_df, !gene %in% key_genes), wilcox_df, "DEGs")

# ── Boxplot: SpaNorm expression per sample, positive-count cells only ────────
make_expr_boxplot <- function(df, tag) {
  n_genes <- length(unique(df$gene))
  if (n_genes == 0) { message("No genes for ", tag, " - skipping expression boxplot"); return(invisible(NULL)) }
  n_col <- min(4, n_genes)
  n_row <- ceiling(n_genes / n_col)
  p <- ggplot(df, aes(sample, expression, fill = category)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ gene, scales = "free_y", ncol = n_col) +
    scale_fill_manual(values = c(CB = "steelblue", DT = "firebrick")) +
    labs(title = paste0("SpaNorm expression, positive-count cells only (", tag, ")"),
         x = NULL, y = "SpaNorm expression") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text  = element_text(size = 9, face = "bold"))
  ggsave(file.path(png_dir, paste0("3_boxplot_spanorm_expression_", tag, ".png")),
         p, width = n_col * 3.5, height = n_row * 2.8, dpi = 200, limitsize = FALSE)
}

spanorm_pos <- spanorm_df %>% filter(pos)
spanorm_pos$gene <- factor(spanorm_pos$gene, levels = genes[genes %in% spanorm_pos$gene])
make_expr_boxplot(filter(spanorm_pos,  gene %in% key_genes), "AR_FOLH1")
make_expr_boxplot(filter(spanorm_pos, !gene %in% key_genes), "DEGs")

# ── Individual per-gene comparison for the AR and NE marker panels ────────────
# Box-only (DT vs CB), one figure per panel, facets ordered as the panel is listed.
plot_panel <- function(panel, tag) {
  df <- prop_df %>% filter(gene %in% panel) %>%
    mutate(gene = factor(gene, levels = intersect(panel, unique(as.character(gene)))))
  make_prop_plots(df, wilcox_df, tag, plots = "box")
}
plot_panel(ar_genes, "AR_genes")
plot_panel(ne_genes, "NE_genes")

# ── Genes with significant Wilcoxon (p < 0.05) ───────────────────────────────
sig_wilcox <- wilcox_df %>% filter(!is.na(p), p < 0.05) %>% arrange(p)
writeLines(as.character(sig_wilcox$gene),
           file.path(outdir, "significant_wilcoxon_genes.txt"))
write.csv(sig_wilcox, file.path(outdir, "significant_wilcoxon_genes.csv"),
          row.names = FALSE)
cat("\nSignificant Wilcoxon genes (p<0.05):", nrow(sig_wilcox), "\n")
print(sig_wilcox)

message("\nDone. Results in: ", outdir)
