library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(scPearsonPCA, lib.loc = "~/R_Library/4.5")
library(qs,          lib.loc = "~/R_Library/4.5")
library(qs2)
source("~/VisHD/functions.R")  # filter_artefacts_knn

# ── Paths ──────────────────────────────────────────────────────────────────
paths   <- system("realpath ~/VisHD/LUT-245-*/raw_srt.qs", intern = TRUE)
slides  <- basename(dirname(paths))
out_dir <- path.expand("~/VisHD/1.integrate_raw_cell")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pearson_path <- file.path(out_dir, "integrated_pearson_srt.qs2")

cat("Found", length(paths), "slides:", paste(slides, collapse = ", "), "\n")

# ── Archetype modules ──────────────────────────────────────────────────────
archetype_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(archetype_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
module_names <- names(archetype_module)
module_cols  <- module_names

# ── Normal / tumour marker gene sets ──────────────────────────────────────
source("~/VisHD/normal_markers.R")
tumour_markers <- c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")
ucell_names <- c("tumour_score", "normal_score")
ucell_cols  <- ucell_names

if (file.exists(pearson_path)) {
  cat("Loading precomputed integrated srt from:", pearson_path, "\n")
  srt <- qs_read(pearson_path)
} else {
  # ── Load, per-slide kNN artefact filter, then build merge-ready minis ────
  cat("Loading and merging", length(paths), "slides...\n")
  srt_list <- lapply(seq_along(paths), function(i) {
    cat("  Loading", slides[i], "\n")
    srt_full <- qread(paths[i])
    srt_full <- UpdateSeuratObject(srt_full)
    n_before <- ncol(srt_full)
    srt_full <- filter_artefacts_knn(srt_full, min_neighbours = 5)
    cat("    ", slides[i], ": ", n_before, " -> ", ncol(srt_full),
        " cells after kNN artefact filter\n", sep = "")

    counts   <- GetAssayData(srt_full, assay = "Spatial", layer = "counts")
    meta     <- srt_full@meta.data
    meta$slide <- slides[i]
    coords   <- GetTissueCoordinates(srt_full, which = "centroids")
    meta$x_centroid <- coords[rownames(meta), "x"]
    meta$y_centroid <- coords[rownames(meta), "y"]
    srt_mini <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
    rm(srt_full); gc()
    srt_mini
  })
  srt <- merge(srt_list[[1]], y = srt_list[-1], add.cell.ids = slides)
  rm(srt_list); gc()
  cat("Merged:", ncol(srt), "cells\n")

  DefaultAssay(srt) <- "Spatial"
  srt <- JoinLayers(srt, assay = "Spatial")

  # ── Post-merge QC filter ────────────────────────────────────────────────
  ca_lo <- quantile(srt$cell_area, 0.05, na.rm = TRUE)
  ca_hi <- quantile(srt$cell_area, 0.99, na.rm = TRUE)
  keep  <- srt$nFeature_Spatial > 20 &
           srt$nCount_Spatial   > 50 &
           srt$cell_area >  ca_lo &
           srt$cell_area <  ca_hi
  srt <- subset(srt, cells = colnames(srt)[which(keep)])
  cat("After QC filter (nFeature>20, nCount>50, cell_area in [Q5, Q99]):", ncol(srt), "cells\n")

  srt <- NormalizeData(srt)

  srt <- AddModuleScore(srt, features = archetype_module, name = "module_score")
  old_cols <- paste0("module_score", seq_along(module_names))
  colnames(srt@meta.data)[match(old_cols, colnames(srt@meta.data))] <- module_names

  # ── Pearson residual PCA ─────────────────────────────────────────────────
  counts_mat <- GetAssayData(srt, assay = "Spatial", layer = "counts")
  tc         <- Matrix::colSums(counts_mat)
  srt        <- FindVariableFeatures(srt, nfeatures = 8000)
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
      scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, NA)) +
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

# ── Primary cell-type signatures + Gavish meta-programs ────────────────────
library(readxl)
library(purrr)

# Seminal Vesicle Epithelial Cells
SVEC_marker <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                 "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")

# Gavish et al. meta-programs by cell type — drop "Malignant" sheet
meta_xlsx <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
meta_programs <- set_names(sheetname, sheetname) |>
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  map(~ map(.x, ~ as.character(na.omit(.x))))

# ── TME cell-type annotation pipeline ─────────────────────────────────────
# tme_markers: one vector per cell type — SVEC literature markers, Gavish
# meta-programs per non-malignant cell type, and Gavish Malignant programs as
# the tumour signature. Genes shared across ≥2 cell types are dropped from
# every entry so each marker list is exclusive.
tme_markers <- c(
  list(SVEC = SVEC_marker),
  lapply(meta_programs, function(progs) unique(unlist(progs))),
  list(Malignant = tumour_markers)
)

gene_counts  <- table(unlist(tme_markers))
shared_genes <- names(gene_counts[gene_counts > 1])
tme_markers  <- lapply(tme_markers, setdiff, shared_genes)
cat("tme_markers: dropped", length(shared_genes),
    "shared genes; per-cell-type sizes:",
    paste(names(tme_markers), lengths(tme_markers), sep = "=", collapse = ", "), "\n")

library(msigdbr)
library(stringr)

c8_data <- readRDS("~/VisHD/public_signature/c8_data_human.Rds")

# Map each tme_markers entry to a MSigDB C8 regex; unknown names fall back to
# the cell-type label uppercased so the pipeline still gets a non-empty hit set.
default_search <- function(nm) gsub("[^A-Z0-9]+", "_", toupper(nm))
known_search_terms <- list(
  SVEC          = "SEMINAL_VESICLE",
  `B cells`     = "B_CELL",
  Endothelial   = "ENDOTHELIAL",
  Epithelial    = "EPITHELIAL",
  Fibroblasts   = "FIBROBLAST",
  Macrophages   = "MACROPHAGE",
  `CD4 T cells` = "CD4.*T_CELL",
  `CD8 T cells` = "CD8.*T_CELL",
  Malignant     = "MALIGNANT|CANCER|TUMOR|PROSTATE_CANCER"
)
tme_search_terms <- setNames(
  lapply(names(tme_markers), function(nm)
    if (!is.null(known_search_terms[[nm]])) known_search_terms[[nm]] else default_search(nm)),
  names(tme_markers)
)
extract_c8_genes <- function(df, search_pattern) {
  df %>% filter(str_detect(gs_name, regex(search_pattern, ignore_case = TRUE))) %>%
    pull(gene_symbol) %>% unique()
}
c8_tme_markers <- map(tme_search_terms, ~extract_c8_genes(c8_data, .x))

source("~/VisHD/celltype_annotation_function.R")
if (!"celltype_annotation" %in% colnames(srt@meta.data)) {
  srt <- tme_cluster_annotation_pipeline(
    obj                 = srt,
    tme_markers         = tme_markers,
    secondary_genes     = c8_tme_markers,
    cluster_col         = "pearson_clusters",
    assay               = "Spatial",
    data_slot           = "data",
    expr_min_val        = 0,
    primary_expr_frac   = 0.05,
    secondary_expr_frac = 0.01,
    min_markers         = 3,
    conf_threshold      = 0.2,
    exclusivity_weight  = 0.30,
    detection_min       = 0.01,
    trim                = 0.10
  )
  qs_save(srt, pearson_path)
}

ggsave(file.path(out_dir, "10a_celltype_annotation_DimPlot.png"),
       DimPlot(srt, group.by = "celltype_annotation",
               cols = as.vector(pals::polychrome()), reduction = "pearsonumap") +
         ggtitle("TME cell-type annotation"),
       width = 7, height = 5, dpi = 350)

ggsave(file.path(out_dir, "10b_celltype_annotation_QC.png"),
       FeaturePlot(srt, "secondary_expr_frac", reduction = "pearsonumap") +
         VlnPlot(srt, "nFeature_Spatial", group.by = "celltype_annotation",
                 pt.size = 0) + theme(legend.position = "none"),
       width = 12, height = 5, dpi = 350)

# ── Export to AnnData ──────────────────────────────────────────────────────
h5ad_path <- file.path(out_dir, "integrated_cells.h5ad")
if (!file.exists(h5ad_path)) {
  srt2anndata(srt,
              data_assay = "Spatial",
              save_name  = sub("\\.h5ad$", "", h5ad_path),
              svg_path   = NULL)
} else {
  cat("AnnData already exists, skipping:", h5ad_path, "\n")
}

cat("Done. Outputs in", out_dir, "\n")
