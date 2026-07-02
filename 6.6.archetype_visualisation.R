#!/usr/bin/env Rscript
# 6.6.archetype_visualisation.R
# Visualise the integrated-tumour archetypal analysis results produced by
# 6.5.integrate_tumour_archetype.ipynb (partipy, n = 4 archetypes) against the
# actual Seurat object, in the visual style of 6.2.archetype_downstream.R.
#
# Renders three result families:
#   1. Archetype weights   (per-cell membership on UMAP + spatial + composition)
#   2. DE genes            (top markers per archetype, z-score heatmap)
#   3. Enrichment          (PROGENy pathway, category, pearson cluster)
#
# Input : ~/VisHD/4.5.integrate_tumour_anno/integrated_pearson_srt.qs2
#         ~/VisHD/6.3.integrate_tumour_archetype/*.csv
# Output: ~/VisHD/6.3.integrate_tumour_archetype/viz/

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(qs2)
})

# ── Config ────────────────────────────────────────────────────────────────────
base_dir <- "~/VisHD"
res_dir  <- file.path(base_dir, "6.3.integrate_tumour_archetype")
srt_path <- file.path(base_dir, "4.5.integrate_tumour_anno", "integrated_pearson_srt.qs2")
outdir   <- file.path(res_dir, "viz")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

archetypes <- paste0("Archetype_", 0:3)
top_n_de   <- 15

arch_cols  <- setNames(brewer.pal(length(archetypes), "Set1"), archetypes)

# ── 1. Load Seurat + archetype weights, join by barcode ───────────────────────
message("Loading Seurat object ...")
srt <- qs_read(srt_path)

w <- read.csv(file.path(res_dir, "AA_cell_weights_n4.csv"),
              row.names = 1, check.names = FALSE)

cells <- intersect(colnames(srt), rownames(w))
message(length(cells), " cells overlap between Seurat object and AA weights")
if (length(cells) < 0.99 * ncol(srt))
  warning("Fewer than 99% of Seurat cells matched the AA weight table — check barcodes.")

srt <- subset(srt, cells = cells)
w   <- w[cells, archetypes, drop = FALSE]

# AA_cell_weights_n4.csv is column-normalised; row-normalise to per-cell membership
w_norm <- as.data.frame(w / rowSums(w))

srt <- AddMetaData(srt, w_norm)
srt$archetype_group <- factor(colnames(w_norm)[max.col(w_norm)], levels = archetypes)
srt$max_weight      <- apply(w_norm, 1, max)

# Reduction to plot on (AA was run on pearsonbatchpca)
red <- if ("pearsonbatchumap" %in% Reductions(srt)) "pearsonbatchumap" else "pearsonumap"
message("Using reduction: ", red)

# ════════════════════════════════════════════════════════════════════════════
# Section 1 — Weights
# ════════════════════════════════════════════════════════════════════════════
message("Section 1: weights ...")

# 1a. Per-archetype membership on UMAP
p_feat <- FeaturePlot(srt, features = archetypes, reduction = red,
                      order = TRUE, ncol = 3) &
  scale_colour_gradientn(colours = c("lightgrey", "#FDBB84", "#B2182B")) &
  theme(plot.title = element_text(size = 10))
ggsave(file.path(outdir, "1_weights_umap_featureplot.png"),
       p_feat, width = 12, height = 8, dpi = 150)

# 1b. Argmax archetype group on UMAP
p_grp <- DimPlot(srt, reduction = red, group.by = "archetype_group",
                 cols = arch_cols, raster = FALSE) +
  ggtitle("Archetype (argmax membership)")
ggsave(file.path(outdir, "2_archetype_group_umap.png"),
       p_grp, width = 7, height = 6, dpi = 150)

# 1c. Spatial — faceted ggplot by slide
df <- srt@meta.data[, c("x_centroid", "y_centroid", "slide", "archetype_group", "max_weight")]

p_sp <- ggplot(df, aes(x_centroid, y_centroid, colour = archetype_group)) +
  geom_point(size = 0.1, shape = 16) +
  facet_wrap(~ slide, scales = "free") +
  scale_colour_manual(values = arch_cols) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(legend.position = "right", strip.text = element_text(size = 9)) +
  labs(title = "Archetype (argmax) — spatial", colour = "Archetype")
ggsave(file.path(outdir, "3_archetype_group_spatial.png"),
       p_sp, width = 14, height = 8, dpi = 150)

p_spw <- ggplot(df, aes(x_centroid, y_centroid, colour = max_weight)) +
  geom_point(size = 0.1, shape = 16) +
  facet_wrap(~ slide, scales = "free") +
  scale_colour_viridis_c(option = "magma") +
  theme_void() +
  theme(legend.position = "right", strip.text = element_text(size = 9)) +
  labs(title = "Max archetype membership — spatial", colour = "Max weight")
ggsave(file.path(outdir, "3b_max_weight_spatial.png"),
       p_spw, width = 14, height = 8, dpi = 150)

# 1d. Composition per slide
p_comp <- ggplot(srt@meta.data, aes(x = slide, fill = archetype_group)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = arch_cols) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells", fill = "Archetype",
       title = "Archetype composition per slide")
ggsave(file.path(outdir, "4_archetype_group_by_slide.png"),
       p_comp, width = 8, height = 5, dpi = 150)

# ════════════════════════════════════════════════════════════════════════════
# Section 2 — DE genes
# ════════════════════════════════════════════════════════════════════════════
message("Section 2: DE genes ...")

# 2a. Top genes per archetype from the per-archetype DE tables
de_top <- lapply(seq_along(archetypes) - 1, function(i) {
  f  <- file.path(res_dir, sprintf("DE_Archetype_%d_vs_rest.csv", i))
  de <- read.csv(f, check.names = FALSE)
  de <- de[!grepl("^MT-", de$names), ]
  de <- de[de$pvals_adj < 0.05 & de$logfoldchanges > 0, ]
  de <- de[order(-de$scores), ]
  head(de$names, top_n_de)
})
names(de_top) <- archetypes
de_genes <- unique(unlist(de_top))

# Persist the per-archetype top-gene table
de_tbl <- do.call(rbind, lapply(archetypes, function(a)
  data.frame(archetype = a, rank = seq_along(de_top[[a]]), gene = de_top[[a]])))
write.csv(de_tbl, file.path(outdir, "top_DE_genes.csv"), row.names = FALSE)

# 2b. Archetype-weighted mean z-score expression (rows 0..4 → Archetype_0..4)
arch_expr <- read.csv(file.path(res_dir, "archetype_expression.csv"),
                      row.names = 1, check.names = FALSE)
rownames(arch_expr) <- paste0("Archetype_", rownames(arch_expr))
de_genes <- de_genes[de_genes %in% colnames(arch_expr)]

# Map each gene to the archetype it is a top marker of (first occurrence)
gene_module <- setNames(rep(NA_character_, length(de_genes)), de_genes)
for (a in archetypes) {
  g <- intersect(de_top[[a]], de_genes)
  gene_module[g[is.na(gene_module[g])]] <- a
}
gene_module <- factor(gene_module, levels = archetypes)

expr_mat <- t(as.matrix(arch_expr[archetypes, de_genes]))   # genes × archetypes
col_z    <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

gene_ha <- rowAnnotation(
  Archetype = gene_module,
  col = list(Archetype = arch_cols),
  show_annotation_name = FALSE
)

png(file.path(outdir, "5_DE_top_genes_heatmap.png"),
    width = 1100, height = 1600, res = 150)
draw(Heatmap(
  expr_mat,
  name              = "z-score",
  col               = col_z,
  left_annotation   = gene_ha,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_split         = gene_module,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 6),
  column_names_gp   = gpar(fontsize = 9),
  column_title      = sprintf("Top %d DE genes per archetype (weighted mean z-score)", top_n_de)
))
dev.off()

# ════════════════════════════════════════════════════════════════════════════
# Section 3 — Enrichment
# ════════════════════════════════════════════════════════════════════════════
message("Section 3: enrichment ...")

# 3a. PROGENy pathway enrichment (ULM score + * where p < 0.05)
pw_est <- read.csv(file.path(res_dir, "pathway_enrichment_est.csv"),
                   row.names = 1, check.names = FALSE)
pw_p   <- read.csv(file.path(res_dir, "pathway_enrichment_pval.csv"),
                   row.names = 1, check.names = FALSE)
rownames(pw_est) <- archetypes
rownames(pw_p)   <- archetypes

pw_mat   <- t(as.matrix(pw_est))    # pathways × archetypes
pw_p_mat <- t(as.matrix(pw_p))
col_pw   <- colorRamp2(c(-max(abs(pw_mat)), 0, max(abs(pw_mat))),
                       c("#2166AC", "white", "#B2182B"))

png(file.path(outdir, "6_pathway_enrichment_heatmap.png"),
    width = 900, height = 800, res = 150)
draw(Heatmap(
  pw_mat,
  name              = "ULM score",
  col               = col_pw,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 9),
  column_names_gp   = gpar(fontsize = 9),
  column_title      = "PROGENy pathway enrichment per archetype (* p < 0.05)",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(pw_p_mat[i, j]) && pw_p_mat[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 11, fontface = "bold"))
  }
))
dev.off()

# Helper: proportion heatmap (rows = archetypes; each row sums to 1)
prop_heatmap <- function(csv, title, out_png, w_px) {
  m <- read.csv(file.path(res_dir, csv), row.names = 1, check.names = FALSE)
  rownames(m) <- archetypes
  m <- as.matrix(m)
  col_p <- colorRamp2(c(0, max(m)), c("white", "#08519C"))
  png(file.path(outdir, out_png), width = w_px, height = 600, res = 150)
  draw(Heatmap(
    m,
    name              = "proportion",
    col               = col_p,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 9),
    column_names_gp   = gpar(fontsize = 9),
    column_title      = title,
    cell_fun = function(j, i, x, y, width, height, fill)
      grid.text(sprintf("%.2f", m[i, j]), x, y, gp = gpar(fontsize = 7))
  ))
  dev.off()
}

# 3b. Category enrichment (per-archetype proportions, row-normalised)
prop_heatmap("category_enrichment.csv",
             "Category enrichment per archetype (row-normalised proportions)",
             "7_category_enrichment_heatmap.png", w_px = 700)

# 3c. Pearson-cluster enrichment (per-archetype proportions, row-normalised)
prop_heatmap("cluster_enrichment.csv",
             "Pearson-cluster enrichment per archetype (row-normalised proportions)",
             "8_cluster_enrichment_heatmap.png", w_px = 1400)

# ── Save joined per-cell weights for downstream use ───────────────────────────
out_weights <- cbind(w_norm,
                     archetype_group = as.character(srt$archetype_group),
                     slide = srt$slide)
write.csv(out_weights,
          file.path(outdir, "cell_archetype_weights_normalised.csv"))

message("\nDone. Figures and tables written to ", outdir)
