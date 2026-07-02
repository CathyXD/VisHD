library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)
library(qs,          lib.loc = "~/R_Library/4.5")
library(qs2)
library(pals)

# ── Paths ──────────────────────────────────────────────────────────────────
paths   <- system("realpath ~/VisHD/LUT-245-*/normal/normal_srt.qs2", intern = TRUE)
slides  <- basename(gsub("/normal", "", dirname(paths)))
out_dir <- path.expand("~/VisHD/8.2.normal_cell_integration_scvi")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
scvi_path <- file.path(out_dir, "integrated_scvi_srt.qs2")

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




if (file.exists(scvi_path)) {
  cat("Loading precomputed integrated srt from:", scvi_path, "\n")
  srt <- qs_read(scvi_path)
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

  # Follow the Seurat v5 integration vignette:
  # https://satijalab.org/seurat/articles/seurat5_integration
  # Layers stay split by slide through NormalizeData -> FindVariableFeatures
  # -> ScaleData -> RunPCA -> IntegrateLayers so scVI can model batch effects.
  # JoinLayers runs AFTER the graph/UMAP step (post-integration), so downstream
  # DEG / module scoring sees a single counts/data layer.
  DefaultAssay(srt) <- "Spatial"
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 1e4)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 5000)
  srt <- ScaleData(srt)
  srt <- RunPCA(srt, npcs = 30, verbose = FALSE)

  # ── scVI integration (batch corrected across slides via split layers) ──
  # scVIIntegration runs its own VAE on raw counts, so no orig.reduction is
  # passed (unlike CCA / RPCA / Harmony, which take orig.reduction = "pca").
  srt <- IntegrateLayers(
    object        = srt,
    method        = scVIIntegration,
    new.reduction = "integrated.scvi",
    conda_env     = "/scratch/pawsey1172/sweng/conda/envs/scvi-tools/",
    verbose       = FALSE
  )

  srt <- FindNeighbors(srt, reduction = "integrated.scvi", dims = 1:30,
                       graph.name = c("scvi_nn", "scvi_snn"))
  srt <- FindClusters(srt, resolution = 1, graph.name = "scvi_snn",
                      cluster.name = "scvi_clusters")
  srt <- RunUMAP(srt, reduction = "integrated.scvi", dims = 1:30,
                 reduction.name = "scviumap")

  # Re-join per-slide counts/data layers (vignette order: post-UMAP, before DEG)
  srt[["Spatial"]] <- JoinLayers(srt[["Spatial"]])

  # SCTransform on the joined counts as the normalization for downstream
  # AddModuleScore / cell annotation / DEG. SpaNorm is not used in this script.
  srt <- SCTransform(srt, assay = "Spatial", verbose = FALSE)
  DefaultAssay(srt) <- "SCT"

  srt <- AddModuleScore(srt, features = archetype_module, name = "module_score")
  old_cols <- paste0("module_score", seq_along(module_names))
  colnames(srt@meta.data)[match(old_cols, colnames(srt@meta.data))] <- module_names

  qs_save(srt, scvi_path)
}

# Ensure SCT is the default assay even when the object is loaded from cache
if ("SCT" %in% Assays(srt)) DefaultAssay(srt) <- "SCT"

# ── Tumour / normal module scores ──────────────────────────────────────────
if (!"tumour_score" %in% colnames(srt@meta.data)) {
  normal_modules <- intersect(unlist(all_marker), rownames(srt))
  srt <- AddModuleScore(srt,
    features = list(tumour_markers, normal_modules),
    name     = "marker_score")
  colnames(srt@meta.data)[match(c("marker_score1", "marker_score2"),
                                colnames(srt@meta.data))] <- ucell_names
  qs_save(srt, scvi_path)
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
  map(~ map(.x, ~ as.character(na.omit(.x))))

flat_meta <- unlist(
  lapply(names(meta_programs), function(s)
    setNames(meta_programs[[s]],
             paste(s, names(meta_programs[[s]]), sep = "::"))),
  recursive = FALSE)

primary_signatures <- c(list(SVEC = SVEC_marker), flat_meta)
primary_signatures <- lapply(primary_signatures, intersect, rownames(srt))
primary_signatures <- primary_signatures[lengths(primary_signatures) >= 3]

sig_cols <- paste0("sig.", make.names(names(primary_signatures)))
if (!all(sig_cols %in% colnames(srt@meta.data))) {
  cat("Scoring", length(primary_signatures), "primary signatures...\n")
  srt <- AddModuleScore(srt, features = primary_signatures, name = "primary_sig_")
  added <- paste0("primary_sig_", seq_along(primary_signatures))
  colnames(srt@meta.data)[match(added, colnames(srt@meta.data))] <- sig_cols
  qs_save(srt, scvi_path)
}

# ── QC plots ───────────────────────────────────────────────────────────────
qc_vars <- c("nCount_Spatial", "nFeature_Spatial")

qc_fp <- wrap_plots(
  lapply(qc_vars, function(v)
    FeaturePlot(srt, features = v, reduction = "scviumap", order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue") +
      ggtitle(v) +
      theme(plot.title = element_text(size = 10), legend.key.width = unit(0.3, "cm"))
  ), nrow = 1) +
  plot_annotation(title = "Normal-cell integration (scVI) — QC metrics")
ggsave(file.path(out_dir, "scvi_0a_FeaturePlot_QC.png"), qc_fp,
       width = 8, height = 3.5, dpi = 400)

n_clusters <- length(unique(srt@meta.data$scvi_clusters))
qc_vln <- wrap_plots(
  lapply(qc_vars, function(v)
    VlnPlot(srt, features = v, group.by = "scvi_clusters", pt.size = 0) +
      ggtitle(v) +
      theme(plot.title = element_text(size = 10), axis.title.x = element_blank(),
            legend.position = "none")
  ), ncol = 1) +
  plot_annotation(title = "Normal-cell integration (scVI) — QC by cluster")
ggsave(file.path(out_dir, "scvi_0b_VlnPlot_QC.png"), qc_vln,
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

ggsave(file.path(out_dir, "scvi_1_FeaturePlot_genes.png"),
       wrap_plots(make_fp("scviumap"), nrow = 1) +
         plot_annotation(title = "Normal-cell integration (scVI) — gene expression"),
       width = length(goi_present) * 3, height = 3, dpi = 400)

# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(srt, reduction = "scviumap",
              group.by = "scvi_clusters", label = TRUE, label.size = 3,
              cols = as.vector(polychrome())) +
  ggtitle("Normal-cell integration (scVI) — clusters") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "scvi_2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

# ── Slide layout ──────────────────────────────────────────────────────────
slide_cols <- as.vector(brewer.set1(length(unique(srt$slide))))

dp_s <- DimPlot(srt, reduction = "scviumap", group.by = "slide", label = FALSE,
                cols = slide_cols) +
  ggtitle("UMAP (scVI — batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "scvi_2c_slide_layout.png"),
       dp_s + plot_annotation(title = "Normal-cell integration (scVI) — cells by slide"),
       width = 6, height = 5, dpi = 400)

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

ggsave(file.path(out_dir, "scvi_3_FeaturePlot_modules.png"),
       wrap_plots(make_mod_fp("scviumap"), ncol = 3) +
         plot_annotation(title = "Normal-cell integration (scVI) — archetype scores"),
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

# SVEC (single reduction)
svec_plot <- FeaturePlot(srt, features = "sig.SVEC", reduction = "scviumap", order = TRUE) +
  scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                        midpoint = 0, name = "SVEC") +
  ggtitle("SVEC (scVI)") +
  theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))
ggsave(file.path(out_dir, "scvi_6_FeaturePlot_SVEC.png"),
       svec_plot + plot_annotation(title = "Seminal Vesicle Epithelial Cells signature (scVI)"),
       width = 5, height = 4.5, dpi = 400)

# One PNG per Gavish cell-type sheet
for (sheet in names(meta_programs)) {
  full_keys  <- paste(sheet, names(meta_programs[[sheet]]), sep = "::")
  full_keys  <- intersect(full_keys, names(primary_signatures))
  if (length(full_keys) == 0) next
  short_keys <- sub(paste0("^", sheet, "::"), "", full_keys)
  cols_keys  <- paste0("sig.", make.names(full_keys))

  ncol_g <- min(3, length(full_keys))
  nrow_g <- ceiling(length(full_keys) / ncol_g)
  ggsave(file.path(out_dir, sprintf("scvi_7_meta_%s.png", make.names(sheet))),
         wrap_plots(make_sig_panel("scviumap", short_keys, cols_keys), ncol = ncol_g) +
           plot_annotation(title = sprintf("Meta-programs — %s (scVI)", sheet)),
         width = ncol_g * 3.5, height = nrow_g * 3, dpi = 400, limitsize = FALSE)
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

ggsave(file.path(out_dir, "scvi_5_FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("scviumap"), nrow = 1) +
         plot_annotation(title = "Normal-cell integration (scVI) — tumour/normal scores"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

# ── DEG / GSEA enrichment on scVI clusters ────────────────────────────────
deg_path    <- file.path(out_dir, "deg_scvi_spanorm.Rds")
enrich_path <- file.path(out_dir, "deg_scvi_enrich.Rds")

DefaultAssay(srt) <- "SCT"
Idents(srt)       <- srt@meta.data$scvi_clusters

if (!file.exists(deg_path)) {
  cat("Running FindAllMarkers on scVI-integrated normal cells...\n")
  # PrepSCTFindMarkers recorrects counts so FindMarkers is valid on the SCT assay
  srt <- PrepSCTFindMarkers(srt)
  DEG <- FindAllMarkers(srt, assay = "SCT", test.used = "MAST")
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
  qs_save(srt, scvi_path)
}

ggsave(file.path(out_dir, "scvi_8_FeaturePlot_normal_celltype_modules.png"),
       FeaturePlot(srt, names(all_marker), reduction = "scviumap",
                   cols = c("white", "red"), order = TRUE),
       width = 18, height = 14, dpi = 350, limitsize = FALSE)

# ── Top-5 DEG FeaturePlot per cluster ─────────────────────────────────────
if (nrow(DEG) > 0) {
  top5_genes <- DEG %>%
    filter(p_val_adj < 0.05, abs(pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 5)
  if (nrow(top5_genes) > 0) {
    ggsave(file.path(out_dir, "scvi_9_FeaturePlot_DEG_top5.png"),
           FeaturePlot(srt, unique(top5_genes$gene), reduction = "scviumap"),
           width = 15,
           height = max(6, 3 * ceiling(length(unique(top5_genes$gene)) / 5)),
           dpi = 350, limitsize = FALSE)
  }
}

# ── TME cell-type annotation pipeline ─────────────────────────────────────
tme_markers <- c(
  list(SVEC = SVEC_marker),
  lapply(meta_programs, function(progs) unique(unlist(progs)))
)

library(msigdbr)
library(stringr)

c8_data <- msigdbr(species = "Homo sapiens", category = "C8")

default_search <- function(nm) gsub("[^A-Z0-9]+", "_", toupper(nm))
known_search_terms <- list(
  SVEC          = "EPITHELIAL|SEMINAL_VESICLE",
  `B cells`     = "B_CELL",
  Endothelial   = "ENDOTHELIAL",
  Epithelial    = "EPITHELIAL|MALIGNANT",
  Fibroblasts   = "FIBROBLAST",
  Macrophages   = "MACROPHAGE",
  `CD4 T cells` = "CD4.*T_CELL",
  `CD8 T cells` = "CD8.*T_CELL"
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
    cluster_col         = "scvi_clusters",
    assay               = "SCT",
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
  qs_save(srt, scvi_path)
}

ggsave(file.path(out_dir, "scvi_10a_celltype_annotation_DimPlot.png"),
       DimPlot(srt, group.by = "celltype_annotation",
               cols = as.vector(polychrome()), reduction = "scviumap") +
         ggtitle("TME cell-type annotation (scVI)"),
       width = 7, height = 5, dpi = 350)

ggsave(file.path(out_dir, "scvi_10b_celltype_annotation_QC.png"),
       FeaturePlot(srt, "secondary_expr_frac", reduction = "scviumap") +
         VlnPlot(srt, "nFeature_Spatial", group.by = "celltype_annotation",
                 pt.size = 0) + theme(legend.position = "none"),
       width = 12, height = 5, dpi = 350)

# ── Export to AnnData ──────────────────────────────────────────────────────
h5ad_path <- file.path(out_dir, "integrated_normal_cells_scvi.h5ad")
if (!file.exists(h5ad_path)) {
  srt2anndata(srt,
              data_assay = "SCT",
              save_name  = sub("\\.h5ad$", "", h5ad_path),
              svg_path   = file.path(out_dir, "SVGs_scvi.Rds"))
} else {
  cat("AnnData already exists, skipping:", h5ad_path, "\n")
}

cat("Done. Outputs in", out_dir, "\n")
