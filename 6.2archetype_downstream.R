#!/usr/bin/env Rscript
# 6.2archetype_downstream.R
# Cross-sample downstream analysis of archetypal analysis results.
# Reads per-sample archetype_result/ CSVs + native tumour_srt.qs2, then:
#   - Aggregates archetype expression and pathway enrichment across all 8 samples
#   - Computes archetype–archetype correlations (expression + pathway)
#   - Derives recurrent gene expression modules via hierarchical clustering
#   - Produces ComplexHeatmap visualisations
#   - Imports per-cell archetype weights + argmax group into each tumour Seurat,
#     scores each archetype's top DE genes (AddModuleScore), and renders
#     cross-sample spatial maps (ImageDimPlot/ImageFeaturePlot); exports metadata
# Output: VisHD/6.2archetype_downstream_tumour/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(jsonlite)
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
  library(qs2)
})

# ── Config ────────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir  <- "~/VisHD"
outdir    <- file.path(base_dir, "6.2archetype_downstream_tumour")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

top_n_genes  <- 50    # genes per archetype for module discovery
k_modules    <- 4     # number of recurrent archetype modules to cut
n_display    <- 15    # top genes per module shown in heatmap
top_n_de     <- 30    # top per-archetype DE genes used for AddModuleScore

# Per-cell archetype spatial figures go in a sub-dir
spdir <- file.path(outdir, "spatial")
dir.create(spdir, showWarnings = FALSE, recursive = TRUE)

# Archetype colour palette (supports up to 8 per-sample archetypes)
max_arch  <- 8
arch_cols <- setNames(brewer.pal(max_arch, "Set1"), paste0("Archetype_", 0:(max_arch - 1)))

# ── 1. Load per-sample results ─────────────────────────────────────────────────
arch_expr_list    <- list()
pw_est_list       <- list()
pw_pval_list      <- list()
arch_idx_top_vec  <- c()

# Per-cell archetype accumulators (filled in the per-sample loop below)
group_plots  <- list()   # ImageDimPlot of archetype_group, one per sample
weight_plots <- list()   # weight_plots[[archetype_idx]][[sample]] — cell-weight maps
score_plots  <- list()   # score_plots[[archetype_idx]][[sample]]  — DEG module scores
meta_list    <- list()   # per-sample archetype meta.data for export

for (i in seq_along(samples)) {
  s          <- samples[i]
  result_dir <- file.path(base_dir, s, "tumour", "archetype_result")

  expr_f  <- file.path(result_dir, "archetype_expression.csv")
  pw_e_f  <- file.path(result_dir, "pathway_enrichment_est.csv")
  pw_p_f  <- file.path(result_dir, "pathway_enrichment_pval.csv")
  top_f   <- file.path(result_dir, "arch_idx_top.json")

  if (!file.exists(expr_f)) {
    message("Skipping ", s, " — missing archetype_expression.csv")
    next
  }

  expr   <- read.csv(expr_f,  row.names = 1, check.names = FALSE)
  pw_est <- read.csv(pw_e_f,  row.names = 1, check.names = FALSE)
  pw_p   <- read.csv(pw_p_f,  row.names = 1, check.names = FALSE)

  # Label: {SampleID}_A{arch_idx}
  arch_label <- paste0(s, "_A", rownames(expr))
  rownames(expr)   <- arch_label
  rownames(pw_est) <- arch_label
  rownames(pw_p)   <- arch_label

  arch_expr_list[[s]] <- expr
  pw_est_list[[s]]    <- pw_est
  pw_pval_list[[s]]   <- pw_p

  if (file.exists(top_f)) {
    top_info <- fromJSON(top_f)
    arch_idx_top_vec[s] <- paste0(s, "_A", top_info$arch_idx_top)
  }

  # ── Per-cell archetype import + spatial visualisation ──────────────────────
  # K (chosen archetype count) = number of per-archetype DE tables; the matching
  # per-cell weights live in AA_cell_weights_n{K}.csv.
  tumour_srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  K           <- length(list.files(result_dir, pattern = "^DE_Archetype_[0-9]+_vs_rest\\.csv$"))
  weights_f   <- file.path(result_dir, sprintf("AA_cell_weights_n%d.csv", K))

  if (K == 0 || !file.exists(tumour_srt_f) || !file.exists(weights_f)) {
    message("  Skipping per-cell archetype import for ", s,
            " — missing tumour_srt.qs2 or AA_cell_weights_n", K, ".csv")
    next
  }

  message("Reading tumour_srt.qs2 for ", s, " (K=", K, " archetypes) ...")
  srt <- qs_read(tumour_srt_f)

  # AA weights are written in adata/tumour_srt cell order (their integer obs_names
  # are positional labels, not barcodes) → align by position to the Seurat cells.
  w <- read.csv(weights_f, row.names = 1, check.names = FALSE)
  if (nrow(w) != ncol(srt)) {
    message("  ", s, ": weight rows (", nrow(w), ") != cells (", ncol(srt),
            ") — skipping archetype import")
    rm(srt); gc(); next
  }
  w_norm <- as.data.frame(w / rowSums(w))          # column-normalised → per-cell membership
  w_cols <- paste0("archW_", 0:(K - 1))
  colnames(w_norm) <- w_cols
  rownames(w_norm) <- colnames(srt)
  srt <- AddMetaData(srt, w_norm)
  srt$archetype_group <- factor(paste0("Archetype_", max.col(w_norm) - 1),
                                levels = paste0("Archetype_", 0:(K - 1)))

  # AddModuleScore from each archetype's top DE genes (per sample)
  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  score_cols <- character(0)
  for (i in 0:(K - 1)) {
    de    <- read.csv(file.path(result_dir, sprintf("DE_Archetype_%d_vs_rest.csv", i)),
                      check.names = FALSE)
    de    <- de[!grepl("^MT-", de$names) & de$pvals_adj < 0.05 & de$logfoldchanges > 0, ]
    de    <- de[order(-de$scores), ]
    genes <- intersect(head(de$names, top_n_de), rownames(srt))
    if (length(genes) < 2) next
    srt <- AddModuleScore(srt, features = list(genes), name = sprintf("archDEG_%d", i))
    names(srt@meta.data)[names(srt@meta.data) == sprintf("archDEG_%d1", i)] <-
      sprintf("archDEG_%d", i)
    score_cols <- c(score_cols, sprintf("archDEG_%d", i))
  }

  # Per-sample spatial plots (FOV-based Visium HD → Image* family)
  group_plots[[s]] <- ImageDimPlot(srt, group.by = "archetype_group",
                                   cols = arch_cols, size = 0.4) +
    ggtitle(s) + theme(plot.title = element_text(size = 9), legend.position = "right")

  for (i in 0:(K - 1)) {
    weight_plots[[as.character(i)]][[s]] <-
      ImageFeaturePlot(srt, features = sprintf("archW_%d", i), size = 0.4) +
      ggtitle(paste0(s, "  A", i, " weight")) + theme(plot.title = element_text(size = 8))
    sc <- sprintf("archDEG_%d", i)
    if (sc %in% colnames(srt@meta.data))
      score_plots[[as.character(i)]][[s]] <-
        ImageFeaturePlot(srt, features = sc, size = 0.4) +
        ggtitle(paste0(s, "  A", i, " DEG score")) + theme(plot.title = element_text(size = 8))
  }

  # Collect archetype metadata for export, then free the object
  md         <- srt@meta.data[, c("archetype_group", w_cols, score_cols), drop = FALSE]
  md$barcode <- rownames(md)
  md$slide   <- s
  meta_list[[s]] <- md
  rm(srt); gc()
}

# ── 2. Aggregate across samples ───────────────────────────────────────────────
# Intersect genes across samples for a consistent feature space
common_genes <- Reduce(intersect, lapply(arch_expr_list, colnames))
common_genes <- common_genes[!grepl("^MT-", common_genes)]
message(length(common_genes), " common genes across all samples (MT genes excluded)")

arch_expr_all <- do.call(rbind, lapply(unname(arch_expr_list), function(df) df[, common_genes]))
pw_est_all    <- do.call(rbind, unname(pw_est_list))
pw_pval_all   <- do.call(rbind, unname(pw_pval_list))

# Sample annotation vector aligned to combined rows
sample_annot <- rep(samples, times = sapply(arch_expr_list, nrow))
names(sample_annot) <- rownames(arch_expr_all)

sample_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(samples)),
  samples
)

# ── 3. Archetype–archetype correlations ───────────────────────────────────────
cor_expr    <- cor(t(as.matrix(arch_expr_all)), method = "spearman")
cor_pathway <- cor(t(as.matrix(pw_est_all)),    method = "spearman", use = "pairwise.complete.obs")

# ── 4. Heatmap: correlation by gene expression ────────────────────────────────
col_cor <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

row_ha <- rowAnnotation(
  Sample = sample_annot,
  col    = list(Sample = sample_colors),
  show_annotation_name = FALSE
)
col_ha <- HeatmapAnnotation(
  Sample = sample_annot,
  col    = list(Sample = sample_colors),
  show_annotation_name = FALSE
)

png(file.path(outdir, "1_archetype_cor_expression.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(
  cor_expr,
  name                    = "Spearman r",
  col                     = col_cor,
  left_annotation         = row_ha,
  top_annotation          = col_ha,
  show_row_names          = TRUE,
  show_column_names       = TRUE,
  row_names_gp            = gpar(fontsize = 7),
  column_names_gp         = gpar(fontsize = 7),
  column_title            = "Archetype–archetype correlation (gene expression)",
  clustering_method_rows  = "ward.D2",
  clustering_method_columns = "ward.D2"
))
dev.off()

# ── 5. Heatmap: correlation by pathway enrichment ─────────────────────────────
png(file.path(outdir, "2_archetype_cor_pathway.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(
  cor_pathway,
  name                    = "Spearman r",
  col                     = col_cor,
  left_annotation         = row_ha,
  top_annotation          = col_ha,
  show_row_names          = TRUE,
  show_column_names       = TRUE,
  row_names_gp            = gpar(fontsize = 7),
  column_names_gp         = gpar(fontsize = 7),
  column_title            = "Archetype–archetype correlation (pathway enrichment)",
  clustering_method_rows  = "ward.D2",
  clustering_method_columns = "ward.D2"
))
dev.off()

# ── 6. Side-by-side: expression vs pathway clustering comparison ───────────────
dend_expr    <- as.dendrogram(hclust(as.dist(1 - cor_expr),    method = "ward.D2"))
dend_pathway <- as.dendrogram(hclust(as.dist(1 - cor_pathway), method = "ward.D2"))

ht_expr <- Heatmap(
  cor_expr,
  name                   = "r (expr)",
  col                    = col_cor,
  left_annotation        = row_ha,
  top_annotation         = col_ha,
  cluster_rows           = dend_expr,
  cluster_columns        = dend_expr,
  show_row_names         = TRUE,
  show_column_names      = FALSE,
  row_names_gp           = gpar(fontsize = 7),
  column_title           = "Gene expression",
  row_dend_width         = unit(15, "mm"),
  column_dend_height     = unit(15, "mm")
)

ht_pathway <- Heatmap(
  cor_pathway,
  name                   = "r (pathway)",
  col                    = col_cor,
  top_annotation         = HeatmapAnnotation(Sample = sample_annot, col = list(Sample = sample_colors), show_legend = FALSE, show_annotation_name = FALSE),
  cluster_rows           = dend_pathway,
  cluster_columns        = dend_pathway,
  show_row_names         = TRUE,
  show_column_names      = FALSE,
  row_names_gp           = gpar(fontsize = 7),
  column_title           = "Pathway enrichment",
  row_dend_width         = unit(15, "mm"),
  column_dend_height     = unit(15, "mm")
)

png(file.path(outdir, "3_archetype_cor_comparison.png"), width = 2800, height = 1300, res = 150)
draw(ht_expr + ht_pathway, ht_gap = unit(8, "mm"), merge_legends = FALSE)
dev.off()

# ── 7. Heatmap: pathway enrichment (all archetypes, significant pathways) ──────

sig_paths <- colnames(pw_est_all)[
  apply(pw_pval_all, 2, function(x) any(x < 0.05, na.rm = TRUE))
]
if (length(sig_paths) == 0) sig_paths <- colnames(pw_est_all)

pw_mat    <- as.matrix(pw_est_all[, sig_paths])
pw_p_mat  <- as.matrix(pw_pval_all[, sig_paths])
col_pw    <- colorRamp2(c(-max(abs(pw_mat)), 0, max(abs(pw_mat))), c("#2166AC", "white", "#B2182B"))

png(file.path(outdir, "3_pathway_enrichment_heatmap.png"), width = 1000, height = 1600, res = 150)
draw(Heatmap(
  pw_mat,
  name                    = "ULM score",
  col                     = col_pw,
  left_annotation         = row_ha,
  show_row_names          = TRUE,
  show_column_names       = TRUE,
  row_names_gp            = gpar(fontsize = 7),
  column_names_gp         = gpar(fontsize = 10),
  column_title            = "Pathway enrichment per archetype (* p < 0.05)",
  clustering_method_rows  = "ward.D2",
  clustering_method_columns = "ward.D2",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(pw_p_mat[i, j]) && pw_p_mat[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
))
dev.off()

# ── 8. Expression density — inspect to choose cutoff for module gene selection ─
png(file.path(outdir, "4a_archetype_expr_density.png"), width = 900, height = 600, res = 150)
plot(density(as.vector(as.matrix(arch_expr_all))),
     main = "Archetype expression value density",
     xlab = "Expression (z-score)", ylab = "Density")
abline(v = 0, col = "red", lty = 2)
dev.off()

# ── 8. Recurrent gene expression modules ─────────────────────────────────────
# Cluster archetypes by expression; find consensus top genes per cluster
k_modules = 3
hc         <- hclust(as.dist(1 - cor_expr), method = "ward.D2")
arch_group <- cutree(hc, k = k_modules)

# Genes > 0 per archetype
pos_genes_list <- apply(as.matrix(arch_expr_all), 1, function(x)
  names(x[x > 0]), simplify = FALSE
)

recurrent_modules <- lapply(seq_len(k_modules), function(g) {
  members <- names(arch_group)[arch_group == g]
  gene_freq <- sort(table(unlist(pos_genes_list[members])), decreasing = TRUE)
  threshold <- max(1, ceiling(length(members) * 0.5))
  gene_freq[gene_freq >= threshold]
})
names(recurrent_modules) <- paste0("Module_", seq_len(k_modules))

# Remove genes shared across modules or with negative average expression in module members
shared_genes <- names(which(table(unlist(lapply(recurrent_modules, names))) > 1))
recurrent_modules <- lapply(names(recurrent_modules), function(nm) {
  m       <- recurrent_modules[[nm]]
  m       <- m[!names(m) %in% shared_genes]
  members <- names(arch_group)[arch_group == which(names(recurrent_modules) == nm)]
  genes   <- names(m)[names(m) %in% colnames(arch_expr_all)]
  avg     <- colMeans(as.matrix(arch_expr_all[members, genes, drop = FALSE]))
  m[names(m) %in% names(avg[avg > 0])]
})
names(recurrent_modules) <- paste0("Module_", seq_len(k_modules))

# Print top 10 genes per module
for (nm in names(recurrent_modules)) {
  grp     <- which(names(recurrent_modules) == nm)
  members <- names(arch_group)[arch_group == grp]
  genes   <- names(recurrent_modules[[nm]])
  avg_expr <- sort(colMeans(as.matrix(arch_expr_all[members, genes[genes %in% colnames(arch_expr_all)], drop = FALSE])), decreasing = TRUE)
  cat("\n===", nm, "(n =", length(members), "archetypes) ===\n")
  print(head(avg_expr, 10))
}

# ── 9. Heatmap: recurrent module genes across all archetypes ──────────────────
module_gene_list <- lapply(names(recurrent_modules), function(nm) {
  grp     <- which(names(recurrent_modules) == nm)
  members <- names(arch_group)[arch_group == grp]
  genes   <- names(recurrent_modules[[nm]])
  genes   <- genes[genes %in% colnames(arch_expr_all)]
  avg_expr <- colMeans(as.matrix(arch_expr_all[members, genes, drop = FALSE]))
  names(sort(avg_expr, decreasing = TRUE))[seq_len(min(n_display, length(avg_expr)))]
})
names(module_gene_list) <- names(recurrent_modules)
all_mod_genes <- unique(unlist(module_gene_list))

gene_module_map <- unlist(unname(
  mapply(function(genes, nm) setNames(rep(nm, length(genes)), genes),
         module_gene_list, names(module_gene_list), SIMPLIFY = FALSE)
))

mod_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(k_modules), names(recurrent_modules))

col_z <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

gene_ha <- rowAnnotation(
  Module = gene_module_map[all_mod_genes],
  col    = list(Module = mod_colors),
  show_annotation_name = FALSE
)
arch_ha <- HeatmapAnnotation(
  Sample  = sample_annot,
  Module  = paste0("Module_", arch_group),
  col     = list(Sample = sample_colors, Module = mod_colors),
  show_annotation_name = TRUE
)

expr_sub <- t(as.matrix(arch_expr_all[, all_mod_genes]))

png(file.path(outdir, "4_recurrent_modules_heatmap.png"), width = 1600, height = 900, res = 150)
draw(Heatmap(
  expr_sub,
  name                    = "z-score",
  col                     = col_z,
  left_annotation         = gene_ha,
  top_annotation          = arch_ha,
  show_row_names          = TRUE,
  show_column_names       = TRUE,
  row_names_gp            = gpar(fontsize = 7),
  column_names_gp         = gpar(fontsize = 7),
  column_title            = "Recurrent gene expression modules",
  cluster_rows            = FALSE,
  clustering_method_columns = "ward.D2",
  row_split               = gene_module_map[all_mod_genes],
  column_split            = arch_group
))
dev.off()

# ── 10. Combined: top-variable genes + pathways side-by-side ─────────────────
gene_var      <- apply(as.matrix(arch_expr_all), 2, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(100)]
top_var_genes <- top_var_genes[top_var_genes %in% colnames(arch_expr_all)]

pw_scaled <- scale(pw_mat)

png(file.path(outdir, "5_combined_expression_pathway.png"), width = 2200, height = 1000, res = 150)
draw(
  Heatmap(
    as.matrix(arch_expr_all[, top_var_genes]),
    name                    = "Gene z-score",
    col                     = col_z,
    left_annotation         = rowAnnotation(
      Sample = sample_annot, Module = paste0("Module_", arch_group),
      col = list(Sample = sample_colors, Module = mod_colors),
      show_annotation_name = FALSE
    ),
    show_row_names          = TRUE,
    show_column_names       = FALSE,
    row_names_gp            = gpar(fontsize = 7),
    column_title            = "Top 100 variable genes",
    clustering_method_rows  = "ward.D2",
    row_split               = arch_group
  ) +
  Heatmap(
    pw_scaled,
    name                    = "Pathway\n(scaled)",
    col                     = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    show_row_names          = FALSE,
    show_column_names       = TRUE,
    column_names_gp         = gpar(fontsize = 9),
    column_title            = "PROGENy pathways",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(pw_p_mat[i, sig_paths[j]]) && pw_p_mat[i, sig_paths[j]] < 0.05)
        grid.text("*", x, y, gp = gpar(fontsize = 9, fontface = "bold"))
    }
  ),
  merge_legends = TRUE
)
dev.off()

# ── 11. Save summary tables ───────────────────────────────────────────────────
write.csv(arch_expr_all, file.path(outdir, "archetype_expression_all_samples.csv"))
write.csv(pw_est_all,    file.path(outdir, "pathway_enrichment_all_samples.csv"))
write.csv(cor_expr,      file.path(outdir, "archetype_cor_expression.csv"))
write.csv(cor_pathway,   file.path(outdir, "archetype_cor_pathway.csv"))

module_df <- data.frame(
  archetype = names(arch_group),
  module    = paste0("Module_", arch_group),
  sample    = sample_annot[names(arch_group)],
  is_top    = names(arch_group) %in% arch_idx_top_vec
)
write.csv(module_df, file.path(outdir, "archetype_module_assignments.csv"), row.names = FALSE)

saveRDS(recurrent_modules, file.path(outdir, "recurrent_modules.Rds"))

# ── 12. Per-cell archetype spatial maps + metadata export ─────────────────────
# NOTE: archetype indices are NOT aligned across samples (archetypes are fit
# per-sample), so each per-archetype grid shows that sample's own A{i}.
if (length(group_plots) > 0) {
  ggsave(file.path(spdir, "6_archetype_group_spatial_allsamples.png"),
         wrap_plots(group_plots, ncol = 4),
         width = 22, height = 11, dpi = 150, limitsize = FALSE)
}

arch_idx_all <- sort(unique(as.integer(c(names(weight_plots), names(score_plots)))))
for (i in arch_idx_all) {
  wp <- weight_plots[[as.character(i)]]
  if (length(wp) > 0)
    ggsave(file.path(spdir, sprintf("7_archetype%d_weight_spatial_allsamples.png", i)),
           wrap_plots(wp, ncol = 4), width = 22, height = 11, dpi = 150, limitsize = FALSE)
  sp <- score_plots[[as.character(i)]]
  if (length(sp) > 0)
    ggsave(file.path(spdir, sprintf("8_archetype%d_DEGscore_spatial_allsamples.png", i)),
           wrap_plots(sp, ncol = 4), width = 22, height = 11, dpi = 150, limitsize = FALSE)
}

# Archetype metadata imported into each per-sample Seurat object
if (length(meta_list) > 0) {
  write.csv(dplyr::bind_rows(meta_list),
            file.path(outdir, "archetype_cell_metadata_all_samples.csv"), row.names = FALSE)
  saveRDS(meta_list, file.path(outdir, "archetype_cell_metadata.Rds"))
  for (s in names(meta_list))
    write.csv(meta_list[[s]],
              file.path(outdir, sprintf("archetype_cell_metadata_%s.csv", s)), row.names = FALSE)
}

message("\nDone. Results saved to ", outdir)
