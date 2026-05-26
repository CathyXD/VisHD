#!/usr/bin/env Rscript
# 5.1.aggregate_DT_deg.R
# Aggregate per-sample 5.1.DT_deg.R results across all 8 samples.
# Outputs:
#   1_DEG_recurrent_heatmap.png      — recurrent significant genes, FC coloured, non-sig masked
#   2_Hallmark_enrichment_heatmap.png
#   2_C6_enrichment_heatmap.png
#   2_C5_enrichment_heatmap.png      — recurrent pathways, NES coloured, non-sig masked

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# ── Config ─────────────────────────────────────────────────────────────────────
samples  <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
              "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")
base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "5.1.aggregate_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

min_recur  <- 2    # must be significant in >= this many samples
p_thr      <- 0.05
top_n_degs <- 100   # cap on rows in DEG heatmap

strip_prefix <- function(x) sub("^[^_]+_", "", x)

# ══════════════════════════════════════════════════════════════════════════════
# 1.  DEG heatmap
# ══════════════════════════════════════════════════════════════════════════════
cat("── Loading per-sample DEG results ──\n")
deg_raw <- list()
for (s in samples) {
  f <- file.path(base_dir, s, "tumour", "deg_DTvsCB_overall.Rds")
  if (!file.exists(f)) { message("Missing: ", f); next }
  df <- readRDS(f)
  df$gene <- rownames(df)
  deg_raw[[s]] <- df[, c("gene", "avg_log2FC", "p_val_adj")]
}
cat(length(deg_raw), "samples with DEG results\n")
if (length(deg_raw) < 2) stop("Need >= 2 samples.")

all_genes <- Reduce(union, lapply(deg_raw, `[[`, "gene"))

fc_mat  <- matrix(NA_real_, nrow = length(all_genes), ncol = length(deg_raw),
                  dimnames = list(all_genes, names(deg_raw)))
sig_mat <- matrix(FALSE,    nrow = length(all_genes), ncol = length(deg_raw),
                  dimnames = list(all_genes, names(deg_raw)))

for (s in names(deg_raw)) {
  df <- deg_raw[[s]]
  fc_mat[df$gene, s]  <- df$avg_log2FC
  sig_mat[df$gene, s] <- !is.na(df$p_val_adj) & df$p_val_adj < p_thr
}

# Recurrent significant genes
n_sig_gene  <- rowSums(sig_mat, na.rm = TRUE)
recur_genes <- names(which(n_sig_gene >= min_recur))
cat(length(recur_genes), "genes significant in >=", min_recur, "samples\n")

if (length(recur_genes) == 0) {
  cat("No recurrent DEGs — skipping DEG heatmap\n")
} else {
  # Rank by recurrence desc, then mean |FC| across significant samples
  mean_abs_fc <- sapply(recur_genes, function(g) {
    mean(abs(fc_mat[g, sig_mat[g, ], drop = TRUE]), na.rm = TRUE)
  })
  row_ord <- order(-n_sig_gene[recur_genes], -mean_abs_fc)
  if (length(recur_genes) > top_n_degs) row_ord <- row_ord[seq_len(top_n_degs)]

  sel_genes <- recur_genes[row_ord]
  fc_sub    <- fc_mat[sel_genes, , drop = FALSE]
  sig_sub   <- sig_mat[sel_genes, , drop = FALSE]
  n_sig_ord <- n_sig_gene[sel_genes]

  # Mask non-significant cells
  fc_display <- fc_sub
  fc_display[!sig_sub] <- NA

  # Drop columns where >30% of display values are NA
  na_frac   <- colMeans(is.na(fc_display))
  drop_cols <- names(na_frac)[na_frac > 0.30]
  if (length(drop_cols) > 0) {
    message("DEG heatmap: removing samples with >30% NA — ",
            paste(sprintf("%s (%.0f%%)", drop_cols, na_frac[drop_cols] * 100), collapse = ", "))
    keep_cols  <- setdiff(colnames(fc_display), drop_cols)
    fc_display <- fc_display[, keep_cols, drop = FALSE]
    sig_sub    <- sig_sub[,    keep_cols, drop = FALSE]
  }
  # Drop rows that became all-NA after column removal (would break hclust)
  keep_rows  <- rowSums(!is.na(fc_display)) > 0
  fc_display <- fc_display[keep_rows, , drop = FALSE]
  sig_sub    <- sig_sub[keep_rows,    , drop = FALSE]
  n_sig_ord  <- rowSums(sig_sub, na.rm = TRUE)

  do_clust <- nrow(fc_display) > 1 && ncol(fc_display) > 1

  fc_lim  <- min(max(abs(fc_sub), na.rm = TRUE), 3)
  col_fc  <- colorRamp2(c(-fc_lim, 0, fc_lim), c("#2166AC", "white", "#B2182B"))

  ra <- rowAnnotation(
    `N sig` = anno_barplot(
      n_sig_ord, gp = gpar(fill = "grey40"),
      width = unit(1.5, "cm"),
      axis_param = list(side = "top", gp = gpar(fontsize = 7))
    ),
    annotation_name_side = "top",
    annotation_name_rot  = 0,
    annotation_name_gp   = gpar(fontsize = 8)
  )

  n_rows <- nrow(fc_display)
  png(file.path(outdir, "1_DEG_recurrent_heatmap.png"),
      width = 2000, height = max(800, n_rows * 28 + 300), res = 200)
  draw(Heatmap(
    fc_display,
    name              = "avg log2FC",
    col               = col_fc,
    na_col            = "grey88",
    cluster_rows      = do_clust,
    cluster_columns   = do_clust,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 8),
    column_names_gp   = gpar(fontsize = 9),
    column_names_rot  = 45,
    column_title      = paste0("Recurrent DT vs CB DEGs  (sig in ≥", min_recur,
                               " of ", length(deg_raw), " samples;  grey = non-sig / not tested)"),
    column_title_gp   = gpar(fontsize = 10),
    right_annotation  = ra,
    rect_gp           = gpar(col = "white", lwd = 0.5)
  ))
  dev.off()
  cat("DEG heatmap saved\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2.  Enrichment heatmaps (one per gene set)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Loading per-sample enrichment results ──\n")
enrich_raw <- list()
for (s in samples) {
  f <- file.path(base_dir, s, "tumour", "enrich_DTvsCB_overall.Rds")
  if (!file.exists(f)) { message("Missing enrichment: ", f); next }
  enrich_raw[[s]] <- readRDS(f)
}
cat(length(enrich_raw), "samples with enrichment results\n")

gene_sets <- c("Hallmark", "C6", "C5")

for (gs in gene_sets) {
  cat("\n  Gene set:", gs, "\n")

  # Extract NES + p.adjust per pathway per sample
  nes_frames <- list()
  for (s in names(enrich_raw)) {
    obj <- enrich_raw[[s]][[gs]]
    if (is.null(obj) || nrow(obj@result) == 0) next
    r <- obj@result[, c("Description", "NES", "p.adjust")]
    r <- r[!is.na(r$NES), ]
    if (nrow(r) > 0) nes_frames[[s]] <- r
  }

  if (length(nes_frames) < 2) { cat("  Insufficient samples — skipping\n"); next }

  all_paths <- Reduce(union, lapply(nes_frames, `[[`, "Description"))

  nes_mat  <- matrix(NA_real_, nrow = length(all_paths), ncol = length(nes_frames),
                     dimnames = list(all_paths, names(nes_frames)))
  esig_mat <- matrix(FALSE,    nrow = length(all_paths), ncol = length(nes_frames),
                     dimnames = list(all_paths, names(nes_frames)))

  for (s in names(nes_frames)) {
    df <- nes_frames[[s]]
    nes_mat[df$Description, s]  <- df$NES
    esig_mat[df$Description, s] <- !is.na(df$p.adjust) & df$p.adjust < p_thr
  }

  n_sig_path  <- rowSums(esig_mat, na.rm = TRUE)
  recur_paths <- names(which(n_sig_path >= min_recur))
  cat("  Recurrent pathways:", length(recur_paths), "\n")
  if (length(recur_paths) == 0) { cat("  None — skipping\n"); next }

  # Rank by recurrence desc, then mean |NES|
  mean_abs_nes <- rowMeans(abs(nes_mat[recur_paths, , drop = FALSE]), na.rm = TRUE)
  row_ord      <- order(-n_sig_path[recur_paths], -mean_abs_nes)

  sel_paths <- recur_paths[row_ord]
  nes_sub   <- nes_mat[sel_paths, , drop = FALSE]
  esig_sub  <- esig_mat[sel_paths, , drop = FALSE]
  n_sig_ord <- n_sig_path[sel_paths]

  # Clean row labels
  rownames(nes_sub)  <- strip_prefix(rownames(nes_sub))
  rownames(esig_sub) <- rownames(nes_sub)

  # Mask non-significant cells
  nes_display <- nes_sub
  nes_display[!esig_sub] <- NA

  # Drop columns where >30% of display values are NA
  na_frac   <- colMeans(is.na(nes_display))
  drop_cols <- names(na_frac)[na_frac > 0.30]
  if (length(drop_cols) > 0) {
    message(gs, " heatmap: removing samples with >30% NA — ",
            paste(sprintf("%s (%.0f%%)", drop_cols, na_frac[drop_cols] * 100), collapse = ", "))
    keep_cols   <- setdiff(colnames(nes_display), drop_cols)
    nes_display <- nes_display[, keep_cols, drop = FALSE]
    esig_sub    <- esig_sub[,   keep_cols, drop = FALSE]
  }
  # Drop rows that became all-NA after column removal (would break hclust)
  keep_rows   <- rowSums(!is.na(nes_display)) > 0
  nes_display <- nes_display[keep_rows, , drop = FALSE]
  esig_sub    <- esig_sub[keep_rows,    , drop = FALSE]
  n_sig_ord   <- rowSums(esig_sub, na.rm = TRUE)

  do_clust <- nrow(nes_display) > 1 && ncol(nes_display) > 1

  nes_lim <- min(max(abs(nes_sub), na.rm = TRUE), 3.5)
  col_nes <- colorRamp2(c(-nes_lim, 0, nes_lim), c("#4575B4", "white", "#D73027"))

  ra_e <- rowAnnotation(
    `N sig` = anno_barplot(
      n_sig_ord, gp = gpar(fill = "grey40"),
      width = unit(1.5, "cm"),
      axis_param = list(side = "top", gp = gpar(fontsize = 7))
    ),
    annotation_name_side = "top",
    annotation_name_rot  = 0,
    annotation_name_gp   = gpar(fontsize = 8)
  )

  n_rows <- nrow(nes_display)
  png(file.path(outdir, paste0("2_", gs, "_enrichment_heatmap.png")),
      width = 2400, height = max(600, n_rows * 30 + 350), res = 200)
  draw(Heatmap(
    nes_display,
    name              = "NES",
    col               = col_nes,
    na_col            = "grey88",
    cluster_rows      = do_clust,
    cluster_columns   = do_clust,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 8),
    column_names_gp   = gpar(fontsize = 9),
    column_names_rot  = 45,
    column_title      = paste0(gs, "  —  Recurrent DT vs CB pathways  (sig in ≥", min_recur,
                               " of ", length(nes_frames), " samples;  grey = non-sig / absent)"),
    column_title_gp   = gpar(fontsize = 10),
    right_annotation  = ra_e,
    rect_gp           = gpar(col = "white", lwd = 0.5)
  ))
  dev.off()
  cat("  Saved", gs, "heatmap\n")
}

message("\nDone. Results in: ", outdir)
