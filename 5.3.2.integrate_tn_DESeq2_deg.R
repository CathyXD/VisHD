#!/usr/bin/env Rscript
# 5.3.2.integrate_tn_DESeq2_deg.R
# Unpaired cross-sample pseudobulk DESeq2: DT vs CB across all tumour samples,
# sourced from the per-sample tumour+normal object (tumour_normal_anno_srt.qs2)
# instead of the tumour-only tumour_srt.qs2 used by 5.3.integrate_DESeq2_deg.R.
# Only the unpaired_DT_vs_CB comparison from that script is reproduced here
# (design ~category_bin, no subclone-pairing requirement) — the paired,
# paired_sample, and cohort A/B comparisons are out of scope.
#
# IMPORTANT: on tumour_normal_anno_srt.qs2, `category` is NOT NA for Normal
# cells — it carries a stale CB1/CB2/DT label inherited from an earlier
# pipeline stage (verified: ~36k Normal cells still tagged CB1/CB2/DT).
# Cells are therefore explicitly filtered to `tumour_anno == "Tumour"`
# immediately after loading, before category_bin/sample_subclone are derived,
# so Normal cells cannot leak into the DT-vs-CB pseudobulk groups.
# Pseudobulk unit = sample x subclone x category (subclone is "1"/"2" for
# tumour cells post-filter; the "Normal" subclone string never appears here).

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tibble)
  library(DESeq2, lib.loc = "~/R_Library/4.5")
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(readxl)
  library(purrr)
})
source("~/VisHD/functions.R")

# ── Config ────────────────────────────────────────────────────────────────────
samples  <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
              "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")
base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "5.3.2.tn_DESeq2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))

# Public signature gene sets
archetype_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
# Load only the Malignant sheet; each column is one meta-program, NA-padded.
meta_programs <- read_excel(
  "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx",
  sheet = "Malignant", col_names = TRUE
) |>
  map(~ as.character(na.omit(.x)))

# Convert named gene lists -> clusterProfiler TERM2GENE (term, gene)
list_to_term2gene <- function(lst) {
  data.frame(
    term = rep(names(lst), lengths(lst)),
    gene = unlist(lst, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}
Archetype <- list_to_term2gene(archetype_module)
MetaProgM <- list_to_term2gene(meta_programs)

min_cells <- 10   # minimum cells to keep a pseudobulk sample

# ── 1. Load all tumour+normal objects, keep Tumour cells only ────────────────
srt_list <- list()
for (s in samples) {
  f <- file.path(base_dir, s, "tumour_normal_anno_srt.qs2")
  if (!file.exists(f)) { message("Skipping ", s, " — missing tumour_normal_anno_srt.qs2"); next }
  message("Loading ", s, " ...")
  obj <- qs_read(f)
  obj <- subset(obj, subset = tumour_anno == "Tumour")
  obj$sample_id       <- s
  obj$category_bin    <- ifelse(obj$category == "DT", "DT",
                          ifelse(grepl("^CB", obj$category), "CB", NA_character_))
  obj$sample_subclone <- paste0(s, "_sc", obj$subclone)
  srt_list[[s]] <- obj
}
cat(length(srt_list), "samples loaded\n")

# ── 2. Merge and save integrated object ───────────────────────────────────────
merged_f <- file.path(outdir, "merged_tumour_srt.qs2")
if (file.exists(merged_f)) {
  cat("Loading existing merged object ...\n")
  srt_merged <- qs_read(merged_f)
} else {
  cat("Merging ...\n")
  srt_merged <- merge(
    srt_list[[1]],
    y            = srt_list[-1],
    add.cell.ids = names(srt_list),
    merge.data   = FALSE
  )
  srt_merged <- JoinLayers(srt_merged, assay = "Spatial")
  qs_save(srt_merged, merged_f)
  cat("Integrated object saved.\n")
}
cat("Merged:", ncol(srt_merged), "cells,", nrow(srt_merged), "genes\n")

# ── 3. Pseudobulk aggregation (unpaired universe — no subclone-pairing filter) ─
meta <- srt_merged@meta.data %>%
  filter(!is.na(category_bin), !is.na(subclone))
valid_cells <- rownames(meta)

counts_raw <- GetAssayData(srt_merged, assay = "Spatial", layer = "counts")[, valid_cells]

group_ids <- setNames(
  paste0(meta$sample_subclone, "__", meta$category_bin),
  valid_cells
)
unique_groups <- unique(group_ids)
cat(length(unique_groups), "pseudobulk groups before filtering\n")

# Sum counts per group (sparse-friendly)
pb_counts <- vapply(unique_groups, function(g) {
  cells_g <- names(group_ids)[group_ids == g]
  Matrix::rowSums(counts_raw[, cells_g, drop = FALSE])
}, numeric(nrow(counts_raw)))

pb_meta_all <- data.frame(
  group           = unique_groups,
  sample_subclone = sub("__.*", "", unique_groups),
  category_bin    = sub(".*__", "", unique_groups),
  n_cells         = as.integer(table(group_ids)[unique_groups]),
  row.names       = unique_groups,
  stringsAsFactors = FALSE
) %>% filter(n_cells >= min_cells)
pb_meta_all$sample_id <- sub("_sc.*", "", pb_meta_all$sample_subclone)
pb_counts_all         <- pb_counts[, rownames(pb_meta_all), drop = FALSE]
cat(nrow(pb_meta_all), "pseudobulks for the unpaired DT-vs-CB analysis\n")

if (nrow(pb_meta_all) < 4) stop("Too few pseudobulk samples for DESeq2 — check cell numbers.")

# ── 4. Unpaired DESeq2 + full visualization/GSEA suite ────────────────────────
# Reusable block: DESeq2 (single-factor design) + standard plots + GSEA.
gsea_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5,
                  Archetype = Archetype, MetaProgMalignant = MetaProgM)

run_unpaired_block <- function(counts, meta, group_var, contrast,
                               out_dir, group_cols, title_tag) {
  pdir <- file.path(out_dir, "png")
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)

  if (ncol(counts) < 4 || length(unique(meta[[group_var]])) < 2) {
    cat("[", title_tag, "] skipped — too few pseudobulks / groups\n"); return(invisible(NULL))
  }

  # Lowly-expressed gene filter
  keep   <- rowSums(counts >= 10) >= max(2, floor(ncol(counts) / 5))
  counts <- counts[keep, , drop = FALSE]

  meta[[group_var]] <- factor(meta[[group_var]], levels = c(contrast[3], contrast[2]))
  if (!"sample_id" %in% colnames(meta))
    meta$sample_id <- sub("_sc.*", "", meta$sample_subclone)
  cat("\n[", title_tag, "] ", ncol(counts), " pseudobulks, ", sum(keep), " genes\n", sep = "")

  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = meta,
                                design = stats::as.formula(paste("~", group_var)))
  dds <- DESeq(dds, parallel = FALSE)
  res_df <- as.data.frame(results(dds, contrast = contrast, alpha = 0.05)) %>%
    rownames_to_column("gene") %>% arrange(padj)

  saveRDS(dds,    file.path(out_dir, "deseq2_dds.Rds"))
  saveRDS(res_df, file.path(out_dir, "deseq2_res.Rds"))
  write.csv(res_df, file.path(out_dir, "deseq2_res.csv"), row.names = FALSE)

  n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5)
  cat("  significant DEGs (padj<0.05, |lfc|>0.5):", n_sig, "\n")

  vsd <- vst(dds, blind = TRUE)

  # PCA
  pca <- plotPCA(vsd, intgroup = c(group_var, "sample_id"), returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  p_pca <- ggplot(pca, aes(PC1, PC2, colour = .data[[group_var]], shape = sample_id)) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = group_cols) +
    scale_shape_manual(values = seq_len(length(unique(pca$sample_id)))) +
    labs(title = paste0("PCA — ", title_tag),
         x = paste0("PC1: ", pv[1], "% variance"),
         y = paste0("PC2: ", pv[2], "% variance"),
         colour = group_var, shape = "Sample") +
    theme_classic()
  ggsave(file.path(pdir, "1_PCA.png"), p_pca, width = 9, height = 6, dpi = 200)

  # Volcano
  p_vol <- res_df %>% filter(!is.na(pvalue)) %>%
    mutate(sig       = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.25,
           direction = case_when(
             sig & log2FoldChange > 0 ~ "Up",
             sig & log2FoldChange < 0 ~ "Down",
             TRUE ~ "NS"
           ),
           label = ifelse(sig, gene, NA_character_)) %>%
    ggplot(aes(log2FoldChange, -log10(pvalue + 1e-300), colour = direction, label = label)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(size = 2.5, max.overlaps = 25, na.rm = TRUE) +
    scale_colour_manual(values = c(NS = "grey70", Up = "firebrick", Down = "royalblue")) +
    geom_vline(xintercept = c(-1.25, 1.25), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = title_tag, subtitle = paste0(n_sig, " DEGs  |  positive = ", contrast[2]),
         x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() + theme(legend.position = "none")
  ggsave(file.path(pdir, "2_volcano.png"), p_vol, width = 7, height = 5, dpi = 200)

  # MA
  p_ma <- res_df %>% filter(!is.na(padj)) %>%
    mutate(sig = padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
    ggplot(aes(log10(baseMean + 1), log2FoldChange, colour = sig)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("grey70", "firebrick")) +
    labs(title = paste0("MA — ", title_tag),
         x = "log10(mean counts + 1)", y = "log2 Fold Change") +
    theme_classic() + theme(legend.position = "none")
  ggsave(file.path(pdir, "3_MA.png"), p_ma, width = 6, height = 4, dpi = 200)

  # Heatmap of top DEGs (VST z-scores)
  sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% mutate(score = abs(log2FoldChange) / padj)
  top_up    <- sig_genes %>% filter(log2FoldChange > 0) %>% slice_max(score, n = 30)
  top_down  <- sig_genes %>% filter(log2FoldChange < 0) %>% slice_max(score, n = 30)
  top_genes <- bind_rows(top_up, top_down) %>% pull(gene)
  if (length(top_genes) >= 5) {
    mat <- t(scale(t(assay(vsd)[top_genes, , drop = FALSE])))
    mat[mat >  3] <-  3; mat[mat < -3] <- -3
    ha <- HeatmapAnnotation(group = meta[colnames(mat), group_var],
                            col = list(group = group_cols),
                            annotation_name_side = "left")
    png(file.path(pdir, "4_heatmap_top_DEGs.png"), width = 1600, height = 2000, res = 200)
    draw(Heatmap(mat, name = "z-score",
                 col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")),
                 top_annotation = ha, show_row_names = TRUE, show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                 column_title = paste0("Top ", length(top_genes), " DEGs — ", title_tag),
                 clustering_method_columns = "ward.D2",
                 column_split = meta[colnames(mat), group_var]))
    dev.off()
  }

  # GSEA across collections
  sig_df <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% arrange(desc(log2FoldChange))
  if (nrow(sig_df) > 0) {
    gene_list <- setNames(sig_df$log2FoldChange, sig_df$gene)
    enrich_list <- lapply(gsea_sets, function(gs)
      tryCatch(clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
               error = function(e) NULL, warning = function(w) NULL))
    saveRDS(enrich_list, file.path(out_dir, "enrich.Rds"))
    for (nm in names(enrich_list)) {
      re <- enrich_list[[nm]]
      if (is.null(re) || nrow(re@result) == 0) next
      sig_n <- sum(re@result$p.adjust < 0.05, na.rm = TRUE)
      if (sig_n == 0) next
      ggsave(file.path(pdir, paste0("5_GSEA_", nm, ".pdf")),
             pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = re),
             width = 6, height = 10)
    }
  }

  # GSEA across collections — all DEGs regardless of p-value
  all_df <- res_df %>% filter(!is.na(log2FoldChange)) %>% arrange(desc(log2FoldChange))
  gene_list_all <- setNames(all_df$log2FoldChange, all_df$gene)
  enrich_list_all <- lapply(gsea_sets, function(gs)
    tryCatch(clusterProfiler::GSEA(gene_list_all, TERM2GENE = gs, verbose = FALSE),
             error = function(e) NULL, warning = function(w) NULL))
  saveRDS(enrich_list_all, file.path(out_dir, "enrich_allgenes.Rds"))
  for (nm in names(enrich_list_all)) {
    re <- enrich_list_all[[nm]]
    if (is.null(re) || nrow(re@result) == 0) next
    sig_n <- sum(re@result$p.adjust < 0.05, na.rm = TRUE)
    if (sig_n == 0) next
    ggsave(file.path(pdir, paste0("5b_GSEA_allDEGs_", nm, ".pdf")),
           pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = re),
           width = 6, height = 10)
  }

  # Summary
  summary_df <- data.frame(
    direction = c(paste0(contrast[2], " > ", contrast[3], " (lfc > 0.5)"),
                  paste0(contrast[3], " > ", contrast[2], " (lfc < -0.5)")),
    n_genes   = c(sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >  0.5),
                  sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5)))
  write.csv(summary_df, file.path(out_dir, "DEG_summary.csv"), row.names = FALSE)
  print(summary_df)
  invisible(list(dds = dds, res = res_df))
}

# ── Unpaired DT vs CB across ALL samples (design ~category_bin) ─────────────
run_unpaired_block(
  counts     = pb_counts_all,
  meta       = pb_meta_all,
  group_var  = "category_bin",
  contrast   = c("category_bin", "DT", "CB"),
  out_dir    = outdir,
  group_cols = c(CB = "steelblue", DT = "firebrick"),
  title_tag  = "DT vs CB — unpaired, tumour_normal_anno_srt (Tumour cells only)")

message("\nDone. Results in: ", outdir)
