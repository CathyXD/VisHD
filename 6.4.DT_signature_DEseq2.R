#!/usr/bin/env Rscript
# 6.4.DT_signature_DEseq2.R   (run-once, after 6.4.DT_signature_analysis.r array 1-8)
# Cross-sample pseudobulk DESeq2 on the DT groupdeg Module_group labels.
# Pseudobulk unit = sample x Module_group. Model: ~sample_id + Module_group.
# Every Module_group level (except "Neg" itself) is contrasted against "Neg"
# cells, pulled from a single DESeq2 fit (mirrors 5.3.integrate_DESeq2_deg.R's
# pseudobulk-then-contrast pattern, adapted to Module_group instead of category).
#
# Reads:
#   <sample>/tumour/tumour_srt.qs2                    (raw Spatial counts)
#   6.4.DT_signature_analysis/<sample>/meta.Rds       (per-cell Module_group)
# Writes (to 6.4.DT_signature_analysis/DESeq2/):
#   merged_tumour_srt.qs2                 cached merged object (all samples)
#   pseudobulk_meta.csv                   pseudobulk sample x Module_group table
#   deseq2_dds.Rds                        fitted DESeqDataSet
#   deseq2_res_<group>_vs_Neg.{Rds,csv}   per-contrast results, one per Module_group level
#   enrich_<group>_vs_Neg.Rds             per-contrast GSEA (Hallmark/C6/C5/Archetype/MetaProgMalignant)
#   DEG_summary.csv                       significant-DEG count per contrast
#   png/0_PCA_pseudobulk.png
#   png/volcano_<group>_vs_Neg.png        one per contrast
#   png/heatmap_<group>_vs_Neg.png        top-DEG VST z-score heatmap, one per contrast
#   png/GSEA_<group>_vs_Neg_<geneset>.pdf one per contrast x gene-set collection
#
#   Rscript 6.4.DT_signature_DEseq2.R

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tibble)
  library(DESeq2, lib.loc = "~/R_Library/4.5")
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(readxl)
  library(purrr)
})

# ── Config ────────────────────────────────────────────────────────────────────
samples  <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
              "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")
base_dir <- "~/VisHD"
sig_dir  <- file.path(base_dir, "6.4.DT_signature_analysis")
outdir   <- file.path(sig_dir, "DESeq2")
png_dir  <- file.path(outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(base_dir, "functions.R"))   # pathwayenrich_plot()

min_cells       <- 10   # minimum cells to keep a pseudobulk sample
min_samples_grp <- 2    # a Module_group level needs pseudobulks in >=2 samples to be tested

# Gene sets for per-contrast GSEA (same collections as 5.3.integrate_DESeq2_deg.R)
Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))
archetype_module <- readRDS(file.path(base_dir, "public_signature/clean_module.Rds"))
meta_programs <- read_excel(
  file.path(base_dir, "public_signature/meta_programs_2025-01-29.xlsx"),
  sheet = "Malignant", col_names = TRUE
) |>
  map(~ as.character(na.omit(.x)))

list_to_term2gene <- function(lst) {
  data.frame(
    term = rep(names(lst), lengths(lst)),
    gene = unlist(lst, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}
Archetype <- list_to_term2gene(archetype_module)
MetaProgM <- list_to_term2gene(meta_programs)
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5,
                  Archetype = Archetype, MetaProgMalignant = MetaProgM)

# ── 1. Load all tumour objects + join per-cell Module_group ──────────────────
srt_list <- list()
for (s in samples) {
  srt_f  <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  meta_f <- file.path(sig_dir, s, "meta.Rds")
  if (!file.exists(srt_f))  { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
  if (!file.exists(meta_f)) { message("Skipping ", s, " — missing 6.4.DT_signature_analysis meta.Rds"); next }
  message("Loading ", s, " ...")
  obj  <- qs_read(srt_f)
  meta <- readRDS(meta_f)
  obj$Module_group <- as.character(meta$Module_group[match(colnames(obj), rownames(meta))])
  obj$sample_id    <- s
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

# ── 3. Pseudobulk aggregation (sample x Module_group) ─────────────────────────
meta        <- srt_merged@meta.data %>% filter(!is.na(Module_group))
valid_cells <- rownames(meta)

counts_raw <- GetAssayData(srt_merged, assay = "Spatial", layer = "counts")[, valid_cells]

group_ids     <- setNames(paste0(meta$sample_id, "__", meta$Module_group), valid_cells)
unique_groups <- unique(group_ids)
cat(length(unique_groups), "pseudobulk groups before filtering\n")

# Sum counts per group (sparse-friendly)
pb_counts <- vapply(unique_groups, function(g) {
  cells_g <- names(group_ids)[group_ids == g]
  Matrix::rowSums(counts_raw[, cells_g, drop = FALSE])
}, numeric(nrow(counts_raw)))

pb_meta <- data.frame(
  group        = unique_groups,
  sample_id    = sub("__.*", "", unique_groups),
  Module_group = sub(".*__", "", unique_groups),
  n_cells      = as.integer(table(group_ids)[unique_groups]),
  row.names    = unique_groups,
  stringsAsFactors = FALSE
)
pb_meta   <- pb_meta %>% filter(n_cells >= min_cells)
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
cat(nrow(pb_meta), "pseudobulks retained (>=", min_cells, "cells)\n")

# Keep only Module_group levels seen in >= min_samples_grp distinct samples —
# a level present in only one sample would alias with that sample's blocking
# term and make the ~sample_id + Module_group design matrix rank-deficient.
grp_n_samples <- pb_meta %>%
  dplyr::distinct(sample_id, Module_group) %>%
  dplyr::count(Module_group, name = "n_samples")
keep_groups <- grp_n_samples %>% dplyr::filter(n_samples >= min_samples_grp) %>% dplyr::pull(Module_group)
if (!"Neg" %in% keep_groups) stop("'Neg' Module_group not present in enough samples — cannot run contrasts.")

pb_meta   <- pb_meta %>% filter(Module_group %in% keep_groups)
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
cat(length(keep_groups), "Module_group levels retained:", paste(keep_groups, collapse = ", "), "\n")

if (nrow(pb_meta) < 4) stop("Too few pseudobulk samples for DESeq2 — check cell numbers.")

# Filter lowly expressed genes: ≥10 counts in ≥(n_pseudobulk / 5) samples
keep_genes <- rowSums(pb_counts >= 10) >= max(2, floor(nrow(pb_meta) / 5))
pb_counts  <- pb_counts[keep_genes, ]
cat(sum(keep_genes), "genes retained after count filtering\n")

# ── 4. DESeq2 (single fit; contrasts pulled per Module_group level below) ─────
pb_meta$Module_group <- factor(pb_meta$Module_group, levels = c("Neg", setdiff(keep_groups, "Neg")))
pb_meta$sample_id    <- factor(pb_meta$sample_id)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_counts),
  colData   = pb_meta,
  design    = ~sample_id + Module_group
)
dds <- DESeq(dds, parallel = FALSE)
cat("DESeq2 done\n")

saveRDS(dds, file.path(outdir, "deseq2_dds.Rds"))
write.csv(pb_meta, file.path(outdir, "pseudobulk_meta.csv"), row.names = FALSE)

# PCA of pseudobulk samples
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("Module_group", "sample_id"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = Module_group, shape = sample_id)) +
  geom_point(size = 3) +
  labs(title = "PCA — pseudobulk samples (VST)",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       colour = "Module_group", shape = "Sample") +
  theme_classic()
ggsave(file.path(png_dir, "0_PCA_pseudobulk.png"), p_pca, width = 9, height = 6, dpi = 200)

# ── 5. Module_group vs Neg contrasts ──────────────────────────────────────────
contrast_levels <- setdiff(levels(pb_meta$Module_group), "Neg")
summary_list <- list()

for (lvl in contrast_levels) {
  lvl_safe <- gsub("[^A-Za-z0-9]+", "_", lvl)
  cat("\nContrast:", lvl, "vs Neg\n")

  res    <- results(dds, contrast = c("Module_group", lvl, "Neg"), alpha = 0.05)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)

  saveRDS(res_df, file.path(outdir, sprintf("deseq2_res_%s_vs_Neg.Rds", lvl_safe)))
  write.csv(res_df, file.path(outdir, sprintf("deseq2_res_%s_vs_Neg.csv", lvl_safe)),
            row.names = FALSE)

  n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5)
  cat("  significant DEGs (padj<0.05, |lfc|>0.5):", n_sig, "\n")
  summary_list[[lvl]] <- data.frame(Module_group = lvl, n_sig_DEGs = n_sig)

  # Volcano plot
  p_vol <- res_df %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      sig       = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.25,
      direction = case_when(
        sig & log2FoldChange > 0 ~ "Up",
        sig & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(sig, gene, NA_character_)
    ) %>%
    ggplot(aes(log2FoldChange, -log10(pvalue + 1e-300), colour = direction, label = label)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(size = 2.5, max.overlaps = 25, na.rm = TRUE) +
    scale_colour_manual(values = c(NS = "grey70", Up = "firebrick", Down = "royalblue")) +
    geom_vline(xintercept = c(-1.25, 1.25), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = sprintf("%s vs Neg — pseudobulk DESeq2", lvl),
         subtitle = paste0(n_sig, " DEGs  |  positive = ", lvl),
         x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() + theme(legend.position = "none")
  ggsave(file.path(png_dir, sprintf("volcano_%s_vs_Neg.png", lvl_safe)),
         p_vol, width = 7, height = 5, dpi = 200)

  # Heatmap of top DEGs (VST z-scores), restricted to this contrast's Neg + lvl samples
  sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% mutate(score = abs(log2FoldChange) / padj)
  top_up    <- sig_genes %>% filter(log2FoldChange > 0) %>% slice_max(score, n = 30)
  top_down  <- sig_genes %>% filter(log2FoldChange < 0) %>% slice_max(score, n = 30)
  top_genes <- bind_rows(top_up, top_down) %>% pull(gene)

  if (length(top_genes) >= 5) {
    hm_samples <- rownames(pb_meta)[pb_meta$Module_group %in% c("Neg", lvl)]
    mat_hm <- assay(vsd)[top_genes, hm_samples, drop = FALSE]
    mat_hm <- t(scale(t(mat_hm)))
    mat_hm[mat_hm >  3] <-  3
    mat_hm[mat_hm < -3] <- -3

    col_fun   <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
    grp_col   <- setNames(c("grey70", "firebrick"), c("Neg", lvl))
    samp_lvls <- levels(pb_meta$sample_id)
    samp_col  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(samp_lvls)), samp_lvls)

    ha <- HeatmapAnnotation(
      Module_group = pb_meta[hm_samples, "Module_group"],
      Sample       = pb_meta[hm_samples, "sample_id"],
      col          = list(Module_group = grp_col, Sample = samp_col),
      annotation_name_side = "left"
    )

    png(file.path(png_dir, sprintf("heatmap_%s_vs_Neg.png", lvl_safe)),
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
      column_title      = sprintf("Top %d DEGs — %s vs Neg", length(top_genes), lvl),
      clustering_method_columns = "ward.D2",
      column_split      = droplevels(pb_meta[hm_samples, "Module_group"])
    ))
    dev.off()
    cat("  heatmap saved\n")
  }

  # GSEA enrichment on this contrast's significant genes
  sig_df <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% arrange(desc(log2FoldChange))
  if (nrow(sig_df) > 0) {
    gene_list   <- setNames(sig_df$log2FoldChange, sig_df$gene)
    enrich_list <- lapply(gene_sets, function(gs) {
      tryCatch(
        clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
        error   = function(e) NULL,
        warning = function(w) NULL
      )
    })
    saveRDS(enrich_list, file.path(outdir, sprintf("enrich_%s_vs_Neg.Rds", lvl_safe)))

    for (nm in names(enrich_list)) {
      res_e <- enrich_list[[nm]]
      if (is.null(res_e) || nrow(res_e@result) == 0) next
      sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
      if (sig_n == 0) next
      p_e <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e)
      ggsave(file.path(png_dir, sprintf("GSEA_%s_vs_Neg_%s.pdf", lvl_safe, nm)),
             p_e, width = 6, height = 10)
    }
    cat("  GSEA done\n")
  } else {
    cat("  no significant DEGs — skipping GSEA\n")
  }
}

summary_df <- bind_rows(summary_list) %>% arrange(desc(n_sig_DEGs))
write.csv(summary_df, file.path(outdir, "DEG_summary.csv"), row.names = FALSE)
print(summary_df)

cat("\nDone. Outputs in", outdir, "\n")
