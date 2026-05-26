#!/usr/bin/env Rscript
# 5.2.integrate_DESeq2_deg.R
# Cross-sample pseudobulk DESeq2: DT vs CB across all tumour samples.
# Pseudobulk unit = sample x subclone x category.
# Model: ~sample_subclone + category  (sample_subclone = paired blocking factor).

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
})
source("~/VisHD/functions.R")

# ── Config ────────────────────────────────────────────────────────────────────
samples  <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
              "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")
base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "5.2.DESeq2_results")
png_dir  <- file.path(outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))

min_cells <- 10   # minimum cells to keep a pseudobulk sample

# ── 1. Load all tumour objects ────────────────────────────────────────────────
srt_list <- list()
for (s in samples) {
  f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(f)) { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
  message("Loading ", s, " ...")
  obj <- qs_read(f)
  obj$sample_id       <- s
  obj$category_bin    <- ifelse(obj$category == "DT", "DT",
                          ifelse(grepl("^CB", obj$category), "CB", NA_character_))
  obj$sample_subclone <- paste0(s, "_sc", obj$subclone)
  srt_list[[s]] <- obj
}
cat(length(srt_list), "samples loaded\n")

# ── 2. Merge and save integrated object ───────────────────────────────────────
cat("Merging ...\n")
srt_merged <- merge(
  srt_list[[1]],
  y            = srt_list[-1],
  add.cell.ids = names(srt_list),
  merge.data   = FALSE
)
# Seurat v5: join layers so counts are accessible as a single matrix
srt_merged <- JoinLayers(srt_merged, assay = "Spatial")
cat("Merged:", ncol(srt_merged), "cells,", nrow(srt_merged), "genes\n")

qs_save(srt_merged, file.path(outdir, "merged_tumour_srt.qs2"))
cat("Integrated object saved.\n")

# ── 3. Pseudobulk aggregation ─────────────────────────────────────────────────
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

# ColData for DESeq2
pb_meta <- data.frame(
  group           = unique_groups,
  sample_subclone = sub("__.*", "", unique_groups),
  category_bin    = sub(".*__", "", unique_groups),
  n_cells         = as.integer(table(group_ids)[unique_groups]),
  row.names       = unique_groups,
  stringsAsFactors = FALSE
)

# Keep only groups with enough cells AND both DT+CB present for same sample_subclone
pb_meta <- pb_meta %>% filter(n_cells >= min_cells)
paired   <- names(which(table(pb_meta$sample_subclone) == 2))
pb_meta  <- pb_meta %>% filter(sample_subclone %in% paired)
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
cat(nrow(pb_meta), "pseudobulk samples retained (", length(paired), "pairs)\n")

if (nrow(pb_meta) < 4) stop("Too few pseudobulk samples for DESeq2 — check cell numbers.")

# Filter lowly expressed genes: ≥10 counts in ≥(n_pseudobulk / 5) samples
keep_genes <- rowSums(pb_counts >= 10) >= max(2, floor(nrow(pb_meta) / 5))
pb_counts  <- pb_counts[keep_genes, ]
cat(sum(keep_genes), "genes retained after count filtering\n")

# ── 4. DESeq2 ────────────────────────────────────────────────────────────────
pb_meta$category_bin    <- factor(pb_meta$category_bin, levels = c("CB", "DT"))
pb_meta$sample_subclone <- factor(pb_meta$sample_subclone)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_counts),
  colData   = pb_meta,
  design    = ~sample_subclone + category_bin
)
dds <- DESeq(dds, parallel = FALSE)
cat("DESeq2 done\n")

res <- results(dds, contrast = c("category_bin", "DT", "CB"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

saveRDS(dds,    file.path(outdir, "deseq2_dds.Rds"))
saveRDS(res_df, file.path(outdir, "deseq2_res_DT_vs_CB.Rds"))
write.csv(res_df, file.path(outdir, "deseq2_res_DT_vs_CB.csv"), row.names = FALSE)

n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5)
cat("Significant DEGs (padj<0.05, |lfc|>0.5):", n_sig, "\n")

# ── 5. Visualizations ─────────────────────────────────────────────────────────

# 5a. PCA of pseudobulk samples
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("category_bin", "sample_subclone"),
                    returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample_id <- sub("_sc.*", "", pca_data$sample_subclone)

p_pca <- ggplot(pca_data, aes(PC1, PC2,
                               colour = category_bin,
                               shape  = sample_id,
                               label  = sample_subclone)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c(CB = "steelblue", DT = "firebrick")) +
  scale_shape_manual(values = seq_len(length(unique(pca_data$sample_id)))) +
  labs(title = "PCA — pseudobulk samples (VST)",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       colour = "Category", shape = "Sample") +
  theme_classic()
ggsave(file.path(png_dir, "1_PCA_pseudobulk.png"),
       p_pca, width = 9, height = 6, dpi = 200)

# 5b. Volcano plot
p_vol <- res_df %>%
  filter(!is.na(pvalue)) %>%
  mutate(
    sig   = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5,
    label = ifelse(sig & rank(-abs(log2FoldChange)) <= 30, gene, NA_character_)
  ) %>%
  ggplot(aes(log2FoldChange, -log10(pvalue + 1e-300),
             colour = sig, label = label)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 25, na.rm = TRUE) +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(title = "DT vs CB — pseudobulk DESeq2 (all samples)",
       subtitle = paste0(n_sig, " DEGs  |  positive = DT"),
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "2_volcano_DT_vs_CB.png"),
       p_vol, width = 7, height = 5, dpi = 200)

# 5c. MA plot
p_ma <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
  ggplot(aes(log10(baseMean + 1), log2FoldChange, colour = sig)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  labs(title = "MA plot — DT vs CB",
       x = "log10(mean counts + 1)", y = "log2 Fold Change") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "3_MA_plot.png"),
       p_ma, width = 6, height = 4, dpi = 200)

# 5d. Heatmap of top DEGs (VST z-scores)
top_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  slice_min(padj, n = 60) %>%
  pull(gene)

if (length(top_genes) >= 5) {
  mat_hm <- assay(vsd)[top_genes, , drop = FALSE]
  mat_hm <- t(scale(t(mat_hm)))
  mat_hm[mat_hm >  3] <-  3
  mat_hm[mat_hm < -3] <- -3

  col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
  cat_col <- c(CB = "steelblue", DT = "firebrick")
  sc_lvls <- levels(pb_meta$sample_subclone)
  sc_col  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(sc_lvls)), sc_lvls)

  ha <- HeatmapAnnotation(
    Category       = pb_meta[colnames(mat_hm), "category_bin"],
    SampleSubclone = pb_meta[colnames(mat_hm), "sample_subclone"],
    col            = list(Category = cat_col, SampleSubclone = sc_col),
    annotation_name_side = "left"
  )

  png(file.path(png_dir, "4_heatmap_top_DEGs.png"),
      width = 1600, height = 2000, res = 200)
  draw(Heatmap(
    mat_hm,
    name              = "z-score",
    col               = col_fun,
    top_annotation    = ha,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7),
    column_title      = paste0("Top ", length(top_genes), " DEGs — DT vs CB"),
    clustering_method_columns = "ward.D2",
    column_split      = pb_meta[colnames(mat_hm), "category_bin"]
  ))
  dev.off()
  cat("Heatmap saved\n")
}

# ── 6. GSEA enrichment ────────────────────────────────────────────────────────
sig_df <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(log2FoldChange))

if (nrow(sig_df) > 0) {
  gene_list <- setNames(sig_df$log2FoldChange, sig_df$gene)

  enrich_list <- lapply(list(Hallmark = Hall, C6 = C6, C5 = C5), function(gs) {
    tryCatch(
      clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
  })
  saveRDS(enrich_list, file.path(outdir, "enrich_DT_vs_CB.Rds"))

  for (nm in names(enrich_list)) {
    res_e <- enrich_list[[nm]]
    if (is.null(res_e) || nrow(res_e@result) == 0) next
    sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
    if (sig_n == 0) next
    p_e <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e)
    ggsave(file.path(png_dir, paste0("5_GSEA_", nm, ".pdf")),
           p_e, width = 6, height = 10)
  }
  cat("Enrichment done\n")
} else {
  cat("No significant DEGs — skipping enrichment\n")
}

# ── 7. Summary ────────────────────────────────────────────────────────────────
summary_df <- data.frame(
  direction = c("DT > CB (lfc > 0.5)", "CB > DT (lfc < -0.5)"),
  n_genes   = c(
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >  0.5),
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5)
  )
)
write.csv(summary_df, file.path(outdir, "DEG_summary.csv"), row.names = FALSE)
print(summary_df)

# ── 8. Per-sample GSEA NES heatmap ───────────────────────────────────────────
# For each sample_subclone pair, rank genes by DT-CB difference in VST space,
# run GSEA, collect NES into a pathway × sample matrix, then plot as heatmap.

vst_mat <- assay(vsd)
pairs   <- unique(pb_meta$sample_subclone)

# Per-pair ranking: DT VST − CB VST (all genes)
per_pair_lfc <- vapply(pairs, function(p) {
  dt_col <- rownames(pb_meta)[pb_meta$sample_subclone == p & pb_meta$category_bin == "DT"]
  cb_col <- rownames(pb_meta)[pb_meta$sample_subclone == p & pb_meta$category_bin == "CB"]
  vst_mat[, dt_col[1]] - vst_mat[, cb_col[1]]
}, numeric(nrow(vst_mat)))
rownames(per_pair_lfc) <- rownames(vst_mat)

gsea_collections <- list(Hallmark = Hall, C6 = C6, C5 = C5)

for (coll_nm in names(gsea_collections)) {
  gs <- gsea_collections[[coll_nm]]

  # Run GSEA per pair (pvalueCutoff=1 to capture all NES values)
  gsea_per_pair <- setNames(lapply(pairs, function(p) {
    lfc_vec <- sort(per_pair_lfc[, p], decreasing = TRUE)
    tryCatch(
      clusterProfiler::GSEA(lfc_vec, TERM2GENE = gs, verbose = FALSE,
                             minGSSize = 5, maxGSSize = 500, pvalueCutoff = 1),
      error = function(e) NULL
    )
  }), pairs)

  # Collect pathway universe and build NES + padj matrices
  all_pathways <- unique(unlist(lapply(gsea_per_pair, function(x) {
    if (is.null(x) || nrow(x@result) == 0) character(0) else x@result$ID
  })))
  if (length(all_pathways) == 0) { cat(coll_nm, ": no pathways\n"); next }

  nes_mat  <- matrix(NA_real_, nrow = length(all_pathways), ncol = length(pairs),
                     dimnames = list(all_pathways, pairs))
  padj_mat <- nes_mat
  for (p in pairs) {
    r <- gsea_per_pair[[p]]
    if (is.null(r) || nrow(r@result) == 0) next
    idx <- match(r@result$ID, all_pathways)
    nes_mat [idx, p] <- r@result$NES
    padj_mat[idx, p] <- r@result$p.adjust
  }

  # Keep pathways significant (padj < 0.05) in at least one sample
  sig_rows <- rowSums(!is.na(padj_mat) & padj_mat < 0.05) >= 1
  if (sum(sig_rows) == 0) { cat(coll_nm, ": no significant pathways\n"); next }

  nes_plot <- nes_mat[sig_rows, , drop = FALSE]

  # For large collections cap at top 50 rows by breadth then mean |NES|
  if (nrow(nes_plot) > 50) {
    breadth  <- rowSums(!is.na(padj_mat[sig_rows, , drop = FALSE]) &
                          padj_mat[sig_rows, , drop = FALSE] < 0.05)
    mean_abs <- rowMeans(abs(nes_plot), na.rm = TRUE)
    ord      <- order(-breadth, -mean_abs)
    nes_plot <- nes_plot[ord[seq_len(50)], , drop = FALSE]
  }

  nes_plot[is.na(nes_plot)] <- 0  # NA → 0 for display (pathway not tested in that sample)

  # Clean row names
  rn <- rownames(nes_plot)
  rn <- sub("^HALLMARK_", "", rn)
  rn <- sub("^KEGG_|^REACTOME_|^GOBP_|^GOCC_|^GOMF_", "", rn)
  rn <- gsub("_", " ", rn)
  rownames(nes_plot) <- rn

  # Column annotation: sample_id + subclone
  col_sample   <- sub("_sc.*", "", colnames(nes_plot))
  col_subclone <- sub(".*_sc", "sc", colnames(nes_plot))
  samp_lvls    <- unique(col_sample)
  samp_col     <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(samp_lvls)), samp_lvls)

  ha_col <- HeatmapAnnotation(
    Sample   = col_sample,
    Subclone = col_subclone,
    col      = list(Sample = samp_col),
    annotation_name_side = "left"
  )

  nes_lim <- min(max(abs(nes_plot), na.rm = TRUE), 3)
  col_nes <- colorRamp2(c(-nes_lim, 0, nes_lim), c("#2166AC", "white", "#B2182B"))

  png_h <- max(800,  35 * nrow(nes_plot) + 300)
  png_w <- max(1000, 55 * ncol(nes_plot) + 500)

  png(file.path(png_dir, paste0("6_GSEA_NES_heatmap_", coll_nm, ".png")),
      width = png_w, height = png_h, res = 150)
  draw(Heatmap(
    nes_plot,
    name              = "NES",
    col               = col_nes,
    top_annotation    = ha_col,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 8),
    column_names_gp   = gpar(fontsize = 7),
    column_title      = paste0("GSEA NES — ", coll_nm, "  (DT vs CB per subclone)"),
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    clustering_method_rows    = "ward.D2",
    clustering_method_columns = "ward.D2",
    rect_gp           = gpar(col = "grey85", lwd = 0.5)
  ))
  dev.off()
  cat("GSEA NES heatmap:", coll_nm, "—", nrow(nes_plot), "pathways ×",
      ncol(nes_plot), "samples\n")
}

message("\nDone. Results in: ", outdir)
