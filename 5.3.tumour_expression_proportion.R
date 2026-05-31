#!/usr/bin/env Rscript
# 5.3.tumour_expression_proportion.R
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
outdir   <- file.path(base_dir, "5.3.expression_proportion")
png_dir  <- file.path(outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

# ── Gene set: AR + FOLH1 + all DEGs (padj<0.05, |lfc|>1.25) from 5.2 ─────────
deg_path <- file.path(base_dir, "5.2.DESeq2_results", "deseq2_res_DT_vs_CB.Rds")
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
genes <- unique(c("AR", "FOLH1", top_degs))
cat("Testing", length(genes), "genes total\n")

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

# ── Bar plot - per-sample positive proportion ─────────────────────────────────
n_genes <- length(unique(prop_df$gene))
n_col   <- min(4, n_genes)
n_row   <- ceiling(n_genes / n_col)

p_bar <- ggplot(prop_df, aes(sample, prop, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = n_col) +
  scale_fill_manual(values = c(CB = "steelblue", DT = "firebrick")) +
  labs(title = "Positive expression proportion - DT vs CB (per sample)",
       x = NULL, y = "Positive proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text  = element_text(size = 9, face = "bold"))

ggsave(file.path(png_dir, "1_barplot_positive_proportion.png"), p_bar,
       width = n_col * 3.5, height = n_row * 2.5, dpi = 200, limitsize = FALSE)

# ── Box plot - DT vs CB across samples per gene with paired Wilcoxon p ───────
ann_y <- prop_df %>% group_by(gene) %>%
  summarise(y = max(prop, na.rm = TRUE) * 1.08, .groups = "drop")
ann_df <- left_join(wilcox_df, ann_y, by = "gene")

p_box <- ggplot(prop_df, aes(category, prop, fill = category)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  geom_text(data = ann_df, aes(x = 1.5, y = y, label = p_lab),
            inherit.aes = FALSE, size = 2.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = n_col) +
  scale_fill_manual(values = c(CB = "steelblue", DT = "firebrick")) +
  labs(title = "DT vs CB positive expression proportion (paired Wilcoxon)",
       x = NULL, y = "Positive proportion") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(png_dir, "2_boxplot_DT_vs_CB.png"), p_box,
       width = n_col * 2.8, height = n_row * 2.8, dpi = 200, limitsize = FALSE)

# ── Genes with significant Wilcoxon (p < 0.05) ───────────────────────────────
sig_wilcox <- wilcox_df %>% filter(!is.na(p), p < 0.05) %>% arrange(p)
writeLines(as.character(sig_wilcox$gene),
           file.path(outdir, "significant_wilcoxon_genes.txt"))
write.csv(sig_wilcox, file.path(outdir, "significant_wilcoxon_genes.csv"),
          row.names = FALSE)
cat("\nSignificant Wilcoxon genes (p<0.05):", nrow(sig_wilcox), "\n")
print(sig_wilcox)

message("\nDone. Results in: ", outdir)
