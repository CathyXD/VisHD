#!/usr/bin/env Rscript
# 6.2archetype_downstream.R
# Cross-sample downstream analysis of archetypal analysis results.
# Reads per-sample archetype_result/ CSVs and h5ad files, then:
#   - Converts AnnData → Seurat v5 (per sample)
#   - Aggregates archetype expression and pathway enrichment across all 8 samples
#   - Computes archetype–archetype correlations (expression + pathway)
#   - Derives recurrent gene expression modules via hierarchical clustering
#   - Produces ComplexHeatmap visualisations
# Output: VisHD/archetype_downstream/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(jsonlite)
  library(Seurat)
  library(SeuratObject)
  library(anndataR, lib.loc = "~/R_Library/4.5")
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

# ── Helper: AnnData → Seurat v5 ───────────────────────────────────────────────
h5ad_to_seurat <- function(h5ad_path) {
  adata <- anndataR::read_h5ad(h5ad_path)

  # Raw counts (cells × genes → genes × cells)
  counts <- t(as.matrix(adata$layers[["counts"]]))
  rownames(counts) <- adata$var_names
  colnames(counts) <- adata$obs_names

  # SpaNorm log-normalised data
  lognorm <- t(as.matrix(adata$X))
  rownames(lognorm) <- adata$var_names
  colnames(lognorm) <- adata$obs_names

  srt <- CreateSeuratObject(counts = counts, project = basename(dirname(h5ad_path)))

  # Add normalised assay layer
  srt[["Spatial"]] <- CreateAssay5Object(counts = counts, data = lognorm)
  DefaultAssay(srt) <- "Spatial"

  # Cell metadata (obs)
  meta_cols <- setdiff(colnames(adata$obs), colnames(srt@meta.data))
  srt <- AddMetaData(srt, metadata = adata$obs[, meta_cols, drop = FALSE])

  # Dimensionality reductions
  if (!is.null(adata$obsm[["X_pca"]])) {
    pca_mat <- as.matrix(adata$obsm[["X_pca"]])
    rownames(pca_mat) <- adata$obs_names
    srt[["pca"]] <- CreateDimReducObject(embeddings = pca_mat, key = "PC_", assay = "Spatial")
  }
  if (!is.null(adata$obsm[["X_umap"]])) {
    umap_mat <- as.matrix(adata$obsm[["X_umap"]])
    rownames(umap_mat) <- adata$obs_names
    colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
    srt[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = "Spatial")
  }

  srt
}

# ── 1. Load per-sample results ─────────────────────────────────────────────────
arch_expr_list    <- list()
pw_est_list       <- list()
pw_pval_list      <- list()
arch_idx_top_vec  <- c()
seurat_list       <- list()

for (i in seq_along(samples)) {
  s          <- samples[i]
  result_dir <- file.path(base_dir, s, "tumour", "archetype_result")

  expr_f  <- file.path(result_dir, "archetype_expression.csv")
  pw_e_f  <- file.path(result_dir, "pathway_enrichment_est.csv")
  pw_p_f  <- file.path(result_dir, "pathway_enrichment_pval.csv")
  h5ad_f  <- file.path(result_dir, "archetype_adata.h5ad")
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

  # AnnData → Seurat (skip if cached .rds already exists)
  if (file.exists(h5ad_f)) {
    srt_cache <- file.path(result_dir, "seurat_from_h5ad.rds")
    if (file.exists(srt_cache)) {
      message("Loading cached Seurat for ", s, " ...")
      seurat_list[[s]] <- readRDS(srt_cache)
    } else {
      message("Converting ", s, " h5ad → Seurat ...")
      seurat_list[[s]] <- tryCatch(h5ad_to_seurat(h5ad_f), error = function(e) {
        message("  Skipped (", conditionMessage(e), ")")
        NULL
      })
      if (!is.null(seurat_list[[s]])) saveRDS(seurat_list[[s]], srt_cache)
    }
  }
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
k_modules = 2
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
if (length(seurat_list) > 0) saveRDS(seurat_list, file.path(outdir, "seurat_list.Rds"))

message("\nDone. Results saved to ", outdir)
