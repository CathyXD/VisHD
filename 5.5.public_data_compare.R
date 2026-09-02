#!/usr/bin/env Rscript
# 5.5.public_data_compare.R
# Public RP cohort comparison: combine per-patient pre/post HTSeq counts from
# RP_data/ and run paired DESeq2 (Post vs Pre).
# Design: ~sample + group  (sample = patient, paired blocking factor).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(DESeq2, lib.loc = "~/R_Library/4.5")
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(readxl)
  library(purrr)
  library(pals)
  library(VennDiagram, lib.loc = "~/R_Library/4.5")
})
source("~/VisHD/functions.R")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# ── Config ────────────────────────────────────────────────────────────────────
base_dir <- "~/VisHD"
rp_dir   <- file.path(base_dir, "RP_data")
outdir   <- file.path(base_dir, "5.5.public_data_compare")
png_dir  <- file.path(outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))

# Public signature gene sets (same collections as 5.3.integrate_DESeq2_deg.R)
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

gsea_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5,
                   Archetype = Archetype, MetaProgMalignant = MetaProgM)

# ── 1. Load & combine HTSeq counts ────────────────────────────────────────────
count_files <- list.files(rp_dir, pattern = "\\.sorted\\.txt$", full.names = TRUE)
cat(length(count_files), "count files found\n")

read_htseq <- function(f) {
  d <- read.table(f, header = FALSE, sep = "\t",
                   col.names = c("gene_id", "count"), stringsAsFactors = FALSE)
  d[!grepl("^__", d$gene_id), ]
}

sample_group <- data.frame(
  file   = count_files,
  sample = sub("-(Pre|Post)\\.sorted\\.txt$", "", basename(count_files)),
  group  = sub(".*-(Pre|Post)\\.sorted\\.txt$", "\\1", basename(count_files)),
  stringsAsFactors = FALSE
)
sample_group$colname <- paste0(sample_group$sample, "_", sample_group$group)

count_list <- lapply(count_files, read_htseq)
gene_ids   <- count_list[[1]]$gene_id
stopifnot(all(vapply(count_list, function(d) identical(d$gene_id, gene_ids), logical(1))))

counts_mat <- do.call(cbind, lapply(count_list, `[[`, "count"))
rownames(counts_mat) <- gene_ids
colnames(counts_mat)  <- sample_group$colname
storage.mode(counts_mat) <- "integer"
cat("Combined count matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")

# ── 2. colData ────────────────────────────────────────────────────────────────
col_data <- data.frame(
  sample = factor(sample_group$sample),
  group  = factor(sample_group$group, levels = c("Pre", "Post")),
  row.names = sample_group$colname
)
cat(nrow(col_data), "samples:",
    paste(rownames(col_data), collapse = ", "), "\n")

# ── 3. Low-count gene filter ──────────────────────────────────────────────────
keep_genes <- rowSums(counts_mat >= 10) >= max(2, floor(nrow(col_data) / 5))
counts_mat <- counts_mat[keep_genes, ]
cat(sum(keep_genes), "genes retained after count filtering\n")

# ── 4. DESeq2 ─────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~sample + group
)
dds <- DESeq(dds, parallel = FALSE)
cat("DESeq2 done\n")

res <- results(dds, contrast = c("group", "Post", "Pre"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

saveRDS(dds,    file.path(outdir, "deseq2_dds.Rds"))
saveRDS(res_df, file.path(outdir, "deseq2_res_Post_vs_Pre.Rds"))
write.csv(res_df, file.path(outdir, "deseq2_res_Post_vs_Pre.csv"), row.names = FALSE)

n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.25)
cat("Significant DEGs (padj<0.05, |lfc|>1.25):", n_sig, "\n")

# ── 5. Visualizations ─────────────────────────────────────────────────────────

# 5a. PCA
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("group", "sample"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = group, shape = sample, label = sample)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c(Pre = "steelblue", Post = "firebrick")) +
  scale_shape_manual(values = seq_len(length(unique(pca_data$sample)))) +
  labs(title = "PCA — RP cohort (VST)",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       colour = "Group", shape = "Sample") +
  theme_classic()
ggsave(file.path(png_dir, "1_PCA.png"), p_pca, width = 9, height = 6, dpi = 200)

# 5b. Volcano
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
  labs(title = "Post vs Pre — RP cohort DESeq2",
       subtitle = paste0(n_sig, " DEGs  |  positive = Post"),
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "2_volcano_Post_vs_Pre.png"), p_vol, width = 7, height = 5, dpi = 200)

# 5c. MA plot
p_ma <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1.25) %>%
  ggplot(aes(log10(baseMean + 1), log2FoldChange, colour = sig)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = c(-1.25, 1.25), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  labs(title = "MA plot — Post vs Pre",
       x = "log10(mean counts + 1)", y = "log2 Fold Change") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "3_MA_plot.png"), p_ma, width = 6, height = 4, dpi = 200)

# 5d. Heatmap of top DEGs (VST z-scores) — top 30 up / top 30 down by abs(lfc)/padj
sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% mutate(score = abs(log2FoldChange) / padj)
top_up    <- sig_genes %>% filter(log2FoldChange > 0) %>% slice_max(score, n = 30)
top_down  <- sig_genes %>% filter(log2FoldChange < 0) %>% slice_max(score, n = 30)
top_genes <- bind_rows(top_up, top_down) %>% pull(gene)

if (length(top_genes) >= 5) {
  mat_hm <- assay(vsd)[top_genes, , drop = FALSE]
  mat_hm <- t(scale(t(mat_hm)))
  mat_hm[mat_hm >  3] <-  3
  mat_hm[mat_hm < -3] <- -3

  col_fun   <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
  grp_col   <- c(Pre = "steelblue", Post = "firebrick")
  samp_lvls <- levels(col_data$sample)
  samp_col  <- setNames(as.vector(pals::polychrome(length(samp_lvls))), samp_lvls)

  ha <- HeatmapAnnotation(
    Group  = col_data[colnames(mat_hm), "group"],
    Sample = col_data[colnames(mat_hm), "sample"],
    col    = list(Group = grp_col, Sample = samp_col),
    annotation_name_side = "left"
  )

  png(file.path(png_dir, "4_heatmap_top_DEGs.png"), width = 1600, height = 2000, res = 200)
  draw(Heatmap(
    mat_hm,
    name              = "z-score",
    col               = col_fun,
    top_annotation    = ha,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7),
    column_title      = paste0("Top ", length(top_genes), " DEGs — Post vs Pre"),
    clustering_method_columns = "ward.D2",
    column_split      = col_data[colnames(mat_hm), "group"]
  ))
  dev.off()
  cat("Heatmap saved\n")
}

# ── 6. GSEA enrichment — significant genes ────────────────────────────────────
sig_df <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(log2FoldChange))

if (nrow(sig_df) > 0) {
  gene_list <- setNames(sig_df$log2FoldChange, sig_df$gene)

  enrich_list <- lapply(gsea_sets, function(gs) {
    tryCatch(
      clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
  })
  saveRDS(enrich_list, file.path(outdir, "enrich_Post_vs_Pre.Rds"))

  for (nm in names(enrich_list)) {
    res_e <- enrich_list[[nm]]
    if (is.null(res_e) || nrow(res_e@result) == 0) next
    sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
    if (sig_n == 0) next
    p_e <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e)
    ggsave(file.path(png_dir, paste0("5_GSEA_", nm, ".pdf")), p_e, width = 6, height = 10)
  }
  cat("Enrichment done\n")
} else {
  cat("No significant DEGs — skipping enrichment\n")
}

# ── 6b. GSEA enrichment — all genes regardless of p-value ────────────────────
all_df <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

gene_list_all <- setNames(all_df$log2FoldChange, all_df$gene)

enrich_list_all <- lapply(gsea_sets, function(gs) {
  tryCatch(
    clusterProfiler::GSEA(gene_list_all, TERM2GENE = gs, verbose = FALSE),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
})
saveRDS(enrich_list_all, file.path(outdir, "enrich_Post_vs_Pre_allgenes.Rds"))

for (nm in names(enrich_list_all)) {
  res_e <- enrich_list_all[[nm]]
  if (is.null(res_e) || nrow(res_e@result) == 0) next
  sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
  if (sig_n == 0) next
  p_e <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e)
  ggsave(file.path(png_dir, paste0("5b_GSEA_allDEGs_", nm, ".pdf")), p_e, width = 6, height = 10)
}
cat("Enrichment (all genes, no p-value filter) done\n")

# ── 7. Summary ────────────────────────────────────────────────────────────────
summary_df <- data.frame(
  direction = c("Post > Pre (lfc > 1.25)", "Pre > Post (lfc < -1.25)"),
  n_genes   = c(
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >  1.25),
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1.25)
  )
)
write.csv(summary_df, file.path(outdir, "DEG_summary.csv"), row.names = FALSE)
print(summary_df)

# ── 8. Venn overlap vs 5.3.2 unpaired tn (DT vs CB) ───────────────────────────
# Post -> DT, Pre -> CB: compare RP cohort (Post vs Pre) significant DEGs
# against the internal tumour cohort (5.3.2, unpaired DT vs CB on
# tumour_normal_anno_srt-derived pseudobulks) significant DEGs.
dt_res_file <- file.path(base_dir, "5.3.2.tn_DESeq2", "deseq2_res.Rds")

if (file.exists(dt_res_file)) {
  res_dt_df <- readRDS(dt_res_file)

  post_up <- res_df    %>% filter(!is.na(padj), padj < 0.05, log2FoldChange >  1.25) %>% pull(gene)
  pre_up  <- res_df    %>% filter(!is.na(padj), padj < 0.05, log2FoldChange < -1.25) %>% pull(gene)
  dt_up   <- res_dt_df %>% filter(!is.na(padj), padj < 0.05, log2FoldChange >  1.25) %>% pull(gene)
  cb_up   <- res_dt_df %>% filter(!is.na(padj), padj < 0.05, log2FoldChange < -1.25) %>% pull(gene)

  venn.diagram(
    x         = list(`Post (RP)` = post_up, `DT (LUT)` = dt_up),
    filename  = file.path(png_dir, "6_venn_Post_vs_DT.png"),
    imagetype = "png", height = 800, width = 800, resolution = 200,
    fill = c("firebrick", "tomato"), alpha = 0.5,
    cat.cex = 1, cex = 1.2, main = "Up in Post vs Up in DT"
  )

  venn.diagram(
    x         = list(`Pre (RP)` = pre_up, `CB (LUT)` = cb_up),
    filename  = file.path(png_dir, "7_venn_Pre_vs_CB.png"),
    imagetype = "png", height = 800, width = 800, resolution = 200,
    fill = c("steelblue", "lightblue"), alpha = 0.5,
    cat.cex = 1, cex = 1.2, main = "Up in Pre vs Up in CB"
  )

  overlap_summary <- data.frame(
    comparison = c("Post vs DT (up)", "Pre vs CB (up)"),
    n_set1     = c(length(post_up), length(pre_up)),
    n_set2     = c(length(dt_up), length(cb_up)),
    n_overlap  = c(length(intersect(post_up, dt_up)), length(intersect(pre_up, cb_up)))
  )
  write.csv(overlap_summary, file.path(outdir, "venn_overlap_summary.csv"), row.names = FALSE)
  print(overlap_summary)
  cat("Venn diagrams saved\n")

  # Per-gene overlap table: which gene set(s) each gene came from, and whether it overlaps
  build_gene_overlap_table <- function(set1, set2, name1, name2) {
    all_genes <- union(set1, set2)
    df <- data.frame(
      gene     = all_genes,
      gene_set = ifelse(all_genes %in% set1 & all_genes %in% set2, paste0(name1, " & ", name2),
                  ifelse(all_genes %in% set1, name1, name2)),
      overlap  = ifelse(all_genes %in% set1 & all_genes %in% set2, "overlap", "not_overlap"),
      stringsAsFactors = FALSE
    )
    arrange(df, desc(overlap))
  }

  gene_overlap_post_dt <- build_gene_overlap_table(post_up, dt_up, "Post (RP)", "DT (LUT)")
  gene_overlap_pre_cb  <- build_gene_overlap_table(pre_up, cb_up, "Pre (RP)", "CB (LUT)")

  write.csv(gene_overlap_post_dt, file.path(outdir, "gene_overlap_Post_vs_DT.csv"), row.names = FALSE)
  write.csv(gene_overlap_pre_cb,  file.path(outdir, "gene_overlap_Pre_vs_CB.csv"), row.names = FALSE)
  cat("Per-gene overlap tables saved\n")
} else {
  cat("5.3.2 results not found (", dt_res_file, ") — skipping Venn comparison\n", sep = "")
}

message("\nDone. Results in: ", outdir)
