library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(scPearsonPCA, lib.loc = "~/R_Library/4.5")
library(qs,          lib.loc = "~/R_Library/4.5")
library(qs2)
library(pals)

# ── Paths ──────────────────────────────────────────────────────────────────
paths   <- system("realpath ~/VisHD/LUT-245-*/normal/normal_srt.qs2", intern = TRUE)
slides  <- basename(gsub("/normal", "", dirname(paths)))
out_dir <- path.expand("~/VisHD/7.1.normal_cell_integration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pearson_path <- file.path(out_dir, "integrated_pearson_srt2.qs2")

cat("Found", length(paths), "slides:", paste(slides, collapse = ", "), "\n")

# ── Archetype modules ──────────────────────────────────────────────────────
archetype_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(archetype_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
module_names <- names(archetype_module)
module_cols  <- module_names

# ── Normal / tumour marker gene sets ──────────────────────────────────────
source("~/VisHD/normal_markers.R")
source("~/VisHD/functions.R")
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

  srt <- JoinLayers(srt, assay = "Spatial", layers = c("counts"),
                     new.layer.names = c("counts"))
  DefaultAssay(srt) <- "Spatial"
  removed_tumour_cells <- readRDS("~/VisHD/7.1.normal_cell_integration/removed_tumour_cells.Rds")
  srt <- subset(srt, cells = setdiff(colnames(srt),removed_tumour_cells ))
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 1e4)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 5000)

  srt <- AddModuleScore(srt, features = archetype_module, name = "module_score")
  old_cols <- paste0("module_score", seq_along(module_names))
  colnames(srt@meta.data)[match(old_cols, colnames(srt@meta.data))] <- module_names

  # ── Pearson residual PCA (no batch + batch-corrected by slide) ──────────
  srt <- do.pearson_pca(srt, resolution = 1.5)
  srt <- do.pearson_pca(srt, batch_variable = "slide", resolution = 2)

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
  # For each sheet (a tibble), split into a sublist keyed by column name,
  # dropping NA padding so each program is a clean character vector of genes.
  map(~ map(.x, ~ as.character(na.omit(.x))))

# Flatten meta_programs into one named list "<sheet>::<program>" → genes
flat_meta <- unlist(
  lapply(names(meta_programs), function(s)
    setNames(meta_programs[[s]],
             paste(s, names(meta_programs[[s]]), sep = "::"))),
  recursive = FALSE)

primary_signatures <- c(list(SVEC = SVEC_marker), flat_meta)
# Restrict to genes present, drop signatures with too few mappable genes
primary_signatures <- lapply(primary_signatures, intersect, rownames(srt))
primary_signatures <- primary_signatures[lengths(primary_signatures) >= 3]

sig_cols <- paste0("sig.", make.names(names(primary_signatures)))
if (!all(sig_cols %in% colnames(srt@meta.data))) {
  cat("Scoring", length(primary_signatures), "primary signatures...\n")
  srt <- AddModuleScore(srt, features = primary_signatures, name = "primary_sig_")
  added <- paste0("primary_sig_", seq_along(primary_signatures))
  colnames(srt@meta.data)[match(added, colnames(srt@meta.data))] <- sig_cols
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
  plot_annotation(title = "Normal-cell integration — QC metrics")
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
  plot_annotation(title = "Normal-cell integration — QC by cluster")
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
         plot_annotation(title = "Normal-cell integration — gene expression (no batch correction)"),
       width = length(goi_present) * 3, height = 3, dpi = 400)

ggsave(file.path(out_dir, "1b_FeaturePlot_genes_batch.png"),
       wrap_plots(make_fp("pearsonbatchumap"), nrow = 1) +
         plot_annotation(title = "Normal-cell integration — gene expression (batch corrected)"),
       width = length(goi_present) * 3, height = 3, dpi = 400)

# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(srt, reduction = "pearsonumap",
              group.by = "pearson_clusters", label = TRUE, label.size = 3,
              cols = as.vector(polychrome())) +
  ggtitle("Normal-cell integration — clusters (no batch)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

dp_batch <- DimPlot(srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE, label.size = 3,
                    cols = as.vector(polychrome())) +
  ggtitle("Normal-cell integration — clusters (batch corrected)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2b_DimPlot_clusters_batch.png"), dp_batch,
       width = 6, height = 5, dpi = 400)

# ── Slide layout: UMAP (no batch) + UMAP (batch corrected) ────────────────
slide_cols <- as.vector(brewer.set1(length(unique(srt$slide))))

dp_s <- DimPlot(srt, reduction = "pearsonumap", group.by = "slide", label = FALSE,
                cols = slide_cols) +
  ggtitle("UMAP (no batch)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

dp_s_batch <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "slide",
                      label = FALSE, cols = slide_cols) +
  ggtitle("UMAP (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "2c_slide_layout.png"),
       (dp_s | dp_s_batch) + plot_annotation(title = "Normal-cell integration — cells by slide"),
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
         plot_annotation(title = "Normal-cell integration — archetype scores (no batch)"),
       width = 12, height = 6, dpi = 400)

ggsave(file.path(out_dir, "3b_FeaturePlot_modules_batch.png"),
       wrap_plots(make_mod_fp("pearsonbatchumap"), ncol = 3) +
         plot_annotation(title = "Normal-cell integration — archetype scores (batch corrected)"),
       width = 12, height = 6, dpi = 400)

# ── Primary signature FeaturePlots ─────────────────────────────────────────
make_sig_panel <- function(reduction, titles, cols) {
  lapply(seq_along(titles), function(j) {
    FeaturePlot(srt, features = cols[j], reduction = reduction, order = TRUE) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = 0, name = titles[j]) +
      ggtitle(titles[j]) +
      theme(plot.title = element_text(size = 8), legend.key.width = unit(0.3, "cm"))
  })
}

# SVEC (single panel, both reductions side-by-side)
svec_plots <- lapply(c("pearsonumap", "pearsonbatchumap"), function(red) {
  FeaturePlot(srt, features = "sig.SVEC", reduction = red, order = TRUE) +
    scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                          midpoint = 0, name = "SVEC") +
    ggtitle(if (red == "pearsonumap") "no batch" else "batch corrected") +
    theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))
})
ggsave(file.path(out_dir, "6_FeaturePlot_SVEC.png"),
       wrap_plots(svec_plots, nrow = 1) +
         plot_annotation(title = "Seminal Vesicle Epithelial Cells signature"),
       width = 10, height = 4.5, dpi = 400)

# One PNG per Gavish cell-type sheet, with all programs in that sheet
for (sheet in names(meta_programs)) {
  full_keys  <- paste(sheet, names(meta_programs[[sheet]]), sep = "::")
  full_keys  <- intersect(full_keys, names(primary_signatures))
  if (length(full_keys) == 0) next
  short_keys <- sub(paste0("^", sheet, "::"), "", full_keys)
  cols_keys  <- paste0("sig.", make.names(full_keys))

  for (red in c("pearsonumap", "pearsonbatchumap")) {
    suffix <- if (red == "pearsonumap") "" else "_batch"
    ncol_g <- min(3, length(full_keys))
    nrow_g <- ceiling(length(full_keys) / ncol_g)
    ggsave(file.path(out_dir, sprintf("7_meta_%s%s.png", make.names(sheet), suffix)),
           wrap_plots(make_sig_panel(red, short_keys, cols_keys), ncol = ncol_g) +
             plot_annotation(title = sprintf("Meta-programs — %s%s", sheet,
                                             if (red == "pearsonumap") " (no batch)" else " (batch corrected)")),
           width = ncol_g * 3.5, height = nrow_g * 3, dpi = 400, limitsize = FALSE)
  }
}

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
         plot_annotation(title = "Normal-cell integration — tumour/normal scores (no batch)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

ggsave(file.path(out_dir, "5b_FeaturePlot_tumour_normal_batch.png"),
       wrap_plots(make_ucell_fp("pearsonbatchumap"), nrow = 1) +
         plot_annotation(title = "Normal-cell integration — tumour/normal scores (batch corrected)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

# ── DEG / GSEA enrichment on batch-corrected clusters ─────────────────────
deg_path    <- file.path(out_dir, "deg_spanorm.Rds")
enrich_path <- file.path(out_dir, "deg_enrich.Rds")

DefaultAssay(srt) <- "Spatial"
Idents(srt)       <- srt@meta.data$pearson_clusters_batch

if (!file.exists(deg_path)) {
  cat("Running FindAllMarkers on integrated normal cells...\n")
  DEG <- FindAllMarkers(srt, test.used = "MAST")
  saveRDS(DEG, deg_path)
} else {
  DEG <- readRDS(deg_path)
}

if (!file.exists(enrich_path) && nrow(DEG) > 0) {
  Hall <- readRDS("~/VisHD/Hall.Rds")
  C6   <- readRDS("~/VisHD/C6.Rds")
  C5   <- readRDS("~/VisHD/C5.Rds")
  DEG_filt <- DEG %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
  if (nrow(DEG_filt) > 0) {
    geneList <- lapply(split(DEG_filt, DEG_filt$cluster),
                       function(x) setNames(x$avg_log2FC, x$gene))
    enrichlist <- lapply(list(Hallmark = Hall, C6 = C6, C5 = C5), function(geneset) {
      enrich <- lapply(names(geneList), function(cn) {
        x <- geneList[[cn]]
        tryCatch(clusterProfiler::GSEA(x, TERM2GENE = geneset),
                 error   = function(e) { message("Skipping ", cn, ": ", conditionMessage(e)); NULL },
                 warning = function(w) { message("Warning ",  cn, ": ", conditionMessage(w)); NULL })
      })
      setNames(enrich, names(geneList))
    })
    saveRDS(enrichlist, enrich_path)
    cat("ENRICHMENT DONE\n")
  }
}

# ── Per-cell-type normal module scores (named all_marker entries) ─────────
if (!all(names(all_marker) %in% colnames(srt@meta.data))) {
  srt <- AddModuleScore(srt, features = all_marker)
  added <- paste0("Cluster", seq_along(all_marker))
  colnames(srt@meta.data)[match(added, colnames(srt@meta.data))] <- names(all_marker)
  qs_save(srt, pearson_path)
}

ggsave(file.path(out_dir, "8_FeaturePlot_normal_celltype_modules.png"),
       FeaturePlot(srt, names(all_marker), reduction = "pearsonbatchumap",
                   cols = c("white", "red"), order = TRUE),
       width = 18, height = 14, dpi = 350, limitsize = FALSE)

# ── Top-5 DEG FeaturePlot per cluster ─────────────────────────────────────
if (nrow(DEG) > 0) {
  top5_genes <- DEG %>%
    filter(p_val_adj < 0.05, abs(pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 5)
  if (nrow(top5_genes) > 0) {
    ggsave(file.path(out_dir, "9_FeaturePlot_DEG_top5.png"),
           FeaturePlot(srt, unique(top5_genes$gene), reduction = "pearsonbatchumap"),
           width = 15,
           height = max(6, 3 * ceiling(length(unique(top5_genes$gene)) / 5)),
           dpi = 350, limitsize = FALSE)
  }
}

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

# Columns tme_cluster_annotation_pipeline() writes on `obj` — must be
# preserved across the no-batch / batch double-run.
annot_cols <- c("celltype_annotation", "celltype_confidence",
                "celltype_score_raw", "celltype_runner_up",
                "annotation_source", "cluster_qc_low_frac",
                "primary_expr_frac", "secondary_expr_frac", "qc_label")

rename_annot <- function(srt, suffix) {
  for (col in annot_cols) {
    if (col %in% colnames(srt@meta.data)) {
      srt@meta.data[[paste0(col, "_", suffix)]] <- srt@meta.data[[col]]
      srt@meta.data[[col]] <- NULL
    }
  }
  srt
}

# Migrate a legacy un-suffixed celltype_annotation column (from before the
# nobatch/batch split) so the guards below don't re-run the pipeline.
if ("celltype_annotation" %in% colnames(srt@meta.data) &&
    !"celltype_annotation_nobatch" %in% colnames(srt@meta.data)) {
  srt <- rename_annot(srt, "nobatch")
  qs_save(srt, pearson_path)
}

# ── Annotation on pearson_clusters (no batch) ────────────────────────────
if (!"celltype_annotation_nobatch" %in% colnames(srt@meta.data)) {
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
  srt <- rename_annot(srt, "nobatch")
  qs_save(srt, pearson_path)
}

ggsave(file.path(out_dir, "10a_celltype_annotation_DimPlot.png"),
       DimPlot(srt, group.by = "celltype_annotation_nobatch",
               cols = as.vector(polychrome()), reduction = "pearsonumap") +
         ggtitle("TME cell-type annotation (no batch)"),
       width = 7, height = 5, dpi = 350)

ggsave(file.path(out_dir, "10b_celltype_annotation_QC.png"),
       FeaturePlot(srt, "secondary_expr_frac_nobatch", reduction = "pearsonumap") +
         VlnPlot(srt, "nFeature_Spatial", group.by = "celltype_annotation_nobatch",
                 pt.size = 0) + theme(legend.position = "none"),
       width = 12, height = 5, dpi = 350)

# ── Annotation on pearson_clusters_batch (batch corrected) ───────────────
if (!"celltype_annotation_batch" %in% colnames(srt@meta.data)) {
  srt <- tme_cluster_annotation_pipeline(
    obj                 = srt,
    tme_markers         = tme_markers,
    secondary_genes     = c8_tme_markers,
    cluster_col         = "pearson_clusters_batch",
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
  srt <- rename_annot(srt, "batch")
  qs_save(srt, pearson_path)
}

ggsave(file.path(out_dir, "10c_celltype_annotation_DimPlot_batch.png"),
       DimPlot(srt, group.by = "celltype_annotation_batch",
               cols = as.vector(polychrome()), reduction = "pearsonbatchumap") +
         ggtitle("TME cell-type annotation (batch corrected)"),
       width = 7, height = 5, dpi = 350)

ggsave(file.path(out_dir, "10d_celltype_annotation_QC_batch.png"),
       FeaturePlot(srt, "secondary_expr_frac_batch", reduction = "pearsonbatchumap") +
         VlnPlot(srt, "nFeature_Spatial", group.by = "celltype_annotation_batch",
                 pt.size = 0) + theme(legend.position = "none"),
       width = 12, height = 5, dpi = 350)

# ── Export to AnnData ──────────────────────────────────────────────────────
h5ad_path <- file.path(out_dir, "integrated_normal_cells2.h5ad")
if (!file.exists(h5ad_path)) {
  srt2anndata(srt,
              data_assay = "Spatial",
              save_name  = sub("\\.h5ad$", "", h5ad_path),
              svg_path   = file.path(out_dir, "SVGs.Rds"))
} else {
  cat("AnnData already exists, skipping:", h5ad_path, "\n")
}

cat("Done. Outputs in", out_dir, "\n")
