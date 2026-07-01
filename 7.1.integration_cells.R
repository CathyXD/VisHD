library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(scPearsonPCA, lib.loc = "~/R_Library/4.5")
library(qs,          lib.loc = "~/R_Library/4.5")
library(qs2)

# ── Paths ──────────────────────────────────────────────────────────────────
paths   <- system("realpath ~/VisHD/LUT-245-*/tumour_subclone_srt.qs2", intern = TRUE)
slides  <- basename(dirname(paths))
out_dir <- path.expand("~/VisHD/integration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pearson_path <- file.path(out_dir, "integrated_pearson_srt.qs2")

cat("Found", length(paths), "slides:", paste(slides, collapse = ", "), "\n")

# ── Archetype modules ──────────────────────────────────────────────────────
archetype_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(archetype_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
module_names <- names(archetype_module)
module_cols  <- module_names

# ── Normal / tumour marker gene sets ──────────────────────────────────────
source("~/VisHD/functions.R")      # srt2anndata
source("~/VisHD/normal_markers.R")
tumour_markers <- c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")
ucell_names <- c("tumour_score", "normal_score")
ucell_cols  <- ucell_names

if (file.exists(pearson_path)) {
  cat("Loading precomputed integrated srt from:", pearson_path, "\n")
  srt <- qs_read(pearson_path)
} else {
  # ── Load and merge counts only ───────────────────────────────────────────
  cat("Loading and merging", length(paths), "slides...\n")
  srt_list <- lapply(seq_along(paths), function(i) {
    cat("  Loading", slides[i], "\n")
    srt_full <- qs_read(paths[i])
    srt_full <- subset(srt_full, subset = tumour_anno != "Removed")
    srt_full <- JoinLayers(srt_full, assay = "Spatial")
    counts   <- GetAssayData(srt_full, assay = "Spatial", layer = "counts")
    meta     <- srt_full@meta.data
    meta$slide <- slides[i]
    # GetTissueCoordinates rows are ordered as colnames() with integer rownames
    # (1:n); the real barcode is in coords$cell. Map by cell id — indexing by
    # rownames(meta) (barcodes) silently mismatches/yields NA.
    coords   <- GetTissueCoordinates(srt_full, which = "centroids")
    idx <- match(rownames(meta), coords$cell)
    stopifnot(!anyNA(idx))
    meta$x_centroid <- coords$x[idx]
    meta$y_centroid <- coords$y[idx]
    srt_mini <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
    rm(srt_full); gc()
    srt_mini
  })
  srt <- merge(srt_list[[1]], y = srt_list[-1], add.cell.ids = slides)
  rm(srt_list); gc()
  cat("Merged:", ncol(srt), "cells\n")

  DefaultAssay(srt) <- "Spatial"
  srt <- JoinLayers(srt, assay = "Spatial")
  srt <- NormalizeData(srt)

  srt <- AddModuleScore(srt, features = archetype_module, name = "module_score")
  old_cols <- paste0("module_score", seq_along(module_names))
  colnames(srt@meta.data)[match(old_cols, colnames(srt@meta.data))] <- module_names

  # ── Pearson residual PCA ─────────────────────────────────────────────────
  counts_mat <- GetAssayData(srt, assay = "Spatial", layer = "counts")
  tc         <- Matrix::colSums(counts_mat)
  srt        <- FindVariableFeatures(srt, nfeatures = 5000)
  hvgs       <- VariableFeatures(srt)

  # Without batch correction
  pcaobj <- sparse_quasipoisson_pca_seurat(
    counts_mat[hvgs, ],
    totalcounts = tc,
    grate       = gene_frequency(counts_mat)[hvgs],
    scale.max   = 10, do.scale = TRUE, do.center = TRUE
  )
  umapobj <- make_umap(pcaobj)
  srt[["pearsonpca"]]   <- pcaobj$reduction.data
  srt[["pearsonumap"]]  <- umapobj$ump
  srt[["pearsongraph"]] <- Seurat::as.Graph(umapobj$grph)
  srt <- FindClusters(srt, graph = "pearsongraph")
  srt@meta.data$pearson_clusters <- srt@meta.data$seurat_clusters

  # With batch correction (batch = "slide")
  srt@meta.data$cell_ID <- rownames(srt@meta.data)
  obs <- data.table(srt@meta.data)[, .(cell_ID, slide)]

  genefreq_batch <- gene_frequency(
    counts_mat,
    obs            = obs,
    cellid_colname = "cell_ID",
    batch_variable = "slide"
  )
  pcaobj_batch <- sparse_quasipoisson_pca_seurat_batch(
    counts_mat[hvgs, ],
    totalcounts    = tc,
    grate          = genefreq_batch[hvgs, ],
    obs            = obs,
    batch_variable = "slide",
    cellid_colname = "cell_ID",
    scale.max      = 10, do.scale = TRUE, do.center = TRUE
  )
  umapobj_batch <- make_umap(pcaobj_batch)
  srt[["pearsonbatchpca"]]   <- pcaobj_batch$reduction.data
  srt[["pearsonbatchumap"]]  <- umapobj_batch$ump
  srt[["pearsonbatchgraph"]] <- Seurat::as.Graph(umapobj_batch$grph)
  srt <- FindClusters(srt, graph = "pearsonbatchgraph")
  srt@meta.data$pearson_clusters_batch <- srt@meta.data$seurat_clusters

  qs_save(srt, pearson_path)
}

# ── Tumour / normal module scores ──────────────────────────────────────────
if (!"tumour_score" %in% colnames(srt@meta.data)) {
  normal_modules <- intersect(unlist(all_marker), rownames(srt))
  srt <- AddModuleScore(srt,
    features = list(tumour_markers, normal_modules),
    name     = "marker_score")
  colnames(srt@meta.data)[match(c("marker_score1", "marker_score2"),
                                colnames(srt@meta.data))] <- ucell_names
  qs_save(srt, pearson_path)
}

# ── QC plots ───────────────────────────────────────────────────────────────
qc_vars <- c("nCount_Spatial", "nFeature_Spatial")

qc_fp <- wrap_plots(
  lapply(qc_vars, function(v)
    FeaturePlot(srt, features = v, reduction = "pearsonumap", order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue") +
      ggtitle(v) +
      theme(plot.title = element_text(size = 10), legend.key.width = unit(0.3, "cm"))
  ), nrow = 1) +
  plot_annotation(title = "Integrated — QC metrics")
ggsave(file.path(out_dir, "0a_FeaturePlot_QC.png"), qc_fp,
       width = 8, height = 3.5, dpi = 400)

n_clusters <- nlevels(srt@meta.data$pearson_clusters)
qc_vln <- wrap_plots(
  lapply(qc_vars, function(v)
    VlnPlot(srt, features = v, group.by = "pearson_clusters", pt.size = 0) +
      ggtitle(v) +
      theme(plot.title = element_text(size = 10), axis.title.x = element_blank(),
            legend.position = "none")
  ), ncol = 1) +
  plot_annotation(title = "Integrated — QC by cluster")
ggsave(file.path(out_dir, "0b_VlnPlot_QC.png"), qc_vln,
       width = max(6, n_clusters * 0.4 + 2), height = 7, dpi = 400)

# ── Gene expression plots ──────────────────────────────────────────────────
goi         <- c("AR", "FOLH1", "ASCL1")
goi_present <- intersect(goi, rownames(srt))
if (length(goi_present) < length(goi))
  message("Not found: ", paste(setdiff(goi, goi_present), collapse = ", "))

make_fp <- function(reduction) {
  lapply(goi_present, function(g) {
    FeaturePlot(srt, features = g, reduction = reduction, order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue") +
      ggtitle(g) +
      theme(plot.title = element_text(size = 10), legend.key.width = unit(0.3, "cm"))
  })
}

ggsave(file.path(out_dir, "1a_FeaturePlot_genes.png"),
       wrap_plots(make_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "Integrated — gene expression (no batch correction)"),
       width = length(goi_present) * 3, height = 3, dpi = 400)

ggsave(file.path(out_dir, "1b_FeaturePlot_genes_batch.png"),
       wrap_plots(make_fp("pearsonbatchumap"), nrow = 1) +
         plot_annotation(title = "Integrated — gene expression (batch corrected)"),
       width = length(goi_present) * 3, height = 3, dpi = 400)

# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(srt, reduction = "pearsonumap",
              group.by = "pearson_clusters", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (no batch)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

dp_batch <- DimPlot(srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (batch corrected)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2b_DimPlot_clusters_batch.png"), dp_batch,
       width = 6, height = 5, dpi = 400)

# ── Slide layout: UMAP (no batch) + UMAP (batch corrected) ────────────────
dp_s <- DimPlot(srt, reduction = "pearsonumap", group.by = "slide", label = FALSE) +
  ggtitle("UMAP (no batch)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

dp_s_batch <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "slide",
                      label = FALSE) +
  ggtitle("UMAP (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "2c_slide_layout.png"),
       (dp_s | dp_s_batch) + plot_annotation(title = "Integrated — cells by slide"),
       width = 12, height = 5, dpi = 400)

# ── Module score FeaturePlots ──────────────────────────────────────────────
make_mod_fp <- function(reduction) {
  lapply(seq_along(module_cols), function(j) {
    FeaturePlot(srt, features = module_cols[j], reduction = reduction, order = TRUE) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = 0, name = module_names[j]) +
      ggtitle(module_names[j]) +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))
  })
}

ggsave(file.path(out_dir, "3a_FeaturePlot_modules.png"),
       wrap_plots(make_mod_fp("pearsonumap"), ncol = 3) +
         plot_annotation(title = "Integrated — archetype scores (no batch)"),
       width = 12, height = 6, dpi = 400)

ggsave(file.path(out_dir, "3b_FeaturePlot_modules_batch.png"),
       wrap_plots(make_mod_fp("pearsonbatchumap"), ncol = 3) +
         plot_annotation(title = "Integrated — archetype scores (batch corrected)"),
       width = 12, height = 6, dpi = 400)

# ── Correlation: gene expression vs module scores ──────────────────────────
expr_df <- as.data.frame(t(as.matrix(
  GetAssayData(srt, layer = "data")[goi_present, , drop = FALSE]
)))
colnames(expr_df) <- paste0(goi_present, "_expr")
cor_df <- cbind(expr_df, srt@meta.data[, module_cols, drop = FALSE])
saveRDS(cor_df, file.path(out_dir, "cor_df.Rds"))

cor_df$cell <- rownames(cor_df)
melt_df <- cor_df %>%
  pivot_longer(cols = ends_with("_expr"), names_to = "gene", values_to = "expr") %>%
  mutate(gene = sub("_expr$", "", gene)) %>%
  pivot_longer(cols = all_of(module_cols), names_to = "module", values_to = "score") %>%
  filter(expr > 0) %>%
  mutate(gene   = factor(gene,   levels = goi_present),
         module = factor(module, levels = module_cols)) %>%
  select(cell, gene, expr, module, score)

label_df <- melt_df %>%
  group_by(gene, module) %>%
  summarise(ct = list(cor.test(expr, score, method = "spearman", exact = FALSE)),
            .groups = "drop") %>%
  mutate(rho   = sapply(ct, function(x) round(x$estimate, 3)),
         pv    = sapply(ct, function(x) format.pval(x$p.value, digits = 2, eps = 0.001)),
         label = paste0("ρ = ", rho, "\np = ", pv))

saveRDS(setNames(label_df$ct, paste(label_df$gene, label_df$module, sep = "_")),
        file.path(out_dir, "cor_test_results.Rds"))

cor_page <- ggplot(melt_df, aes(x = expr, y = score)) +
  geom_point(size = 0.15, alpha = 0.25, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.5) +
  geom_text(data = label_df, aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.1, size = 2.5, inherit.aes = FALSE) +
  facet_grid(gene ~ module, scales = "free") +
  labs(x = "Gene expression (log-norm)", y = "Module score",
       title = "Gene × module correlation — Integrated") +
  theme_classic(base_size = 18) +
  theme(strip.text = element_text(size = 18))
ggsave(file.path(out_dir, "4_correlation_genes_vs_modules.png"), cor_page,
       width = length(module_cols) * 3, height = length(goi_present) * 2.8 + 0.5,
       dpi = 400, limitsize = FALSE)

# ── Tumour / normal score FeaturePlots ─────────────────────────────────────
make_ucell_fp <- function(reduction) {
  lapply(seq_along(ucell_cols), function(j) {
    FeaturePlot(srt, features = ucell_cols[j], reduction = reduction, order = TRUE) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = 0, name = ucell_names[j]) +
      ggtitle(ucell_names[j]) +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))
  })
}

ggsave(file.path(out_dir, "5a_FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "Integrated — tumour/normal scores (no batch)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

ggsave(file.path(out_dir, "5b_FeaturePlot_tumour_normal_batch.png"),
       wrap_plots(make_ucell_fp("pearsonbatchumap"), nrow = 1) +
         plot_annotation(title = "Integrated — tumour/normal scores (batch corrected)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

# ── Export to AnnData ──────────────────────────────────────────────────────
h5ad_path <- file.path(out_dir, "integrated_cells.h5ad")
if (!file.exists(h5ad_path)) {
  cat("Exporting to AnnData:", h5ad_path, "\n")
  srt2anndata(srt,
              count_assay = "Spatial",
              data_assay  = "Spatial",
              save_name   = file.path(out_dir, "integrated_cells"),
              svg_path    = NULL)
} else {
  cat("AnnData already exists, skipping:", h5ad_path, "\n")
}

cat("Done. Outputs in", out_dir, "\n")
