suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SpaNorm,    lib.loc = "~/R_Library/4.5")
  library(qs2)
  library(leidenbase, lib.loc = "~/R_Library/4.5")
  library(UCell,      lib.loc = "~/R_Library/4.5")
  library(ggplot2)
  library(readxl)
  library(purrr)
  library(patchwork)
  library(stringr)
  library(pals)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")   # all_marker

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
heat.shock_mod <- c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "NR4A1", "DUSP1", "ZFP36", "HSPA1A","HSPA1B", "IER3","IER2", "BTG2", "SOCS3")
ifp_save <- function(srt, feats, outfile, ncol = 3, dpi = 300, cells_highlight = NULL) {
  feats <- feats[feats %in% colnames(srt@meta.data)]
  if (length(feats) == 0) return(invisible(NULL))
  nrows <- ceiling(length(feats) / ncol)
  p <- dark_feature_plot(srt, features = feats, ncol = ncol, cells_highlight = cells_highlight)
  ggsave(outfile, p,
         width = ncol * 8, height = nrows * 7,
         dpi = dpi, limitsize = FALSE)
}
# gene-expression version of ifp_save: features are pulled from an assay
# (not meta.data). One figure per gene set, one dark_feature_plot panel per gene.
ifp_save_genes <- function(srt, genes, outfile, assay = "SpaNorm",
                           ncol = 3, dpi = 300, cells_highlight = NULL) {
  genes <- genes[genes %in% rownames(GetAssayData(srt, assay = assay, layer = "data"))]
  if (length(genes) == 0) return(invisible(NULL))
  nrows <- ceiling(length(genes) / ncol)
  p <- dark_feature_plot(srt, features = genes, assay = assay, ncol = ncol, cells_highlight = cells_highlight)
  ggsave(outfile, p,
         width = ncol * 8, height = nrows * 7,
         dpi = dpi, limitsize = FALSE)
}

# Key gene sets for spatial expression (one figure each)
ar_genes <- c("AR", "FOLH1", "KLK3", "NKX3-1", "TMPRSS2", "KLK4", "STEAP2", "STEAP1")
ne_genes <- c("CHGA", "CHGB", "SCG2", "SLC18A1", "SYNGR4", "NPB", "PTPN5")

# Prostate epithelial lineage signatures (scored per cell with UCell)
epi_markers <- list(
  Basal   = c("TP63","KRT5","KRT14","KRT15","KRT6A","KRT17","NGFR",
              "DST","BCAM","LGALS1","MMP7","SNAI2","COL17A1","LAMB3"),
  Luminal = c("KRT8","KRT18","KRT19","AR","NKX3-1","KLK3","KLK2","ACP3",
              "MSMB","TMPRSS2","FOLH1","SLC45A3","AZGP1","FOXA1","HOXB13","DPP4")
)

meta_xlsx     <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname     <- excel_sheets(meta_xlsx)

# groupgene: named list of the 3 curated tumour signatures (G1/G2/G3), each a
# character vector of gene symbols (from 6.3.archetype_module.r)
groupgene <- readRDS("~/VisHD/6.2archetype_downstream_tumour/archetype_module/groupgene.Rds")

# ── groupgene signature expression, aggregated across all samples ────────────
# Per sample: average SpaNorm expression per signature gene within each
# final_annotation group (tumour cells only). aggregate_signature_heatmap()
# then combines every sample's matrix into a single heatmap so the same group
# (e.g. G1) can be compared side-by-side across samples.
signature_group_avg <- function(srt, groupgene, annot_col = "final_annotation",
                                assay = "SpaNorm") {
  if ("compartment" %in% colnames(srt@meta.data))
    srt <- srt[, srt$compartment == "Tumour"]
  if (ncol(srt) < 10) return(NULL)

  genes <- intersect(unique(unlist(groupgene, use.names = FALSE)), rownames(srt))
  if (length(genes) == 0) return(NULL)

  expr  <- as.matrix(GetAssayData(srt, assay = assay, layer = "data")[genes, , drop = FALSE])
  group <- as.character(srt@meta.data[[annot_col]])
  group[is.na(group)] <- "NA"

  grp_idx <- split(seq_len(ncol(expr)), group)
  vapply(grp_idx, function(idx) rowMeans(expr[, idx, drop = FALSE]), numeric(nrow(expr)))
}

aggregate_signature_heatmap <- function(avg_list, groupgene, outfile, n_display = 40) {
  avg_list <- avg_list[!vapply(avg_list, is.null, logical(1))]
  if (length(avg_list) == 0) return(invisible(NULL))

  common_genes <- Reduce(intersect, lapply(avg_list, rownames))
  gene_list    <- lapply(groupgene, function(g) intersect(g, common_genes))
  gene_list    <- gene_list[lengths(gene_list) > 0]
  if (length(gene_list) == 0) return(invisible(NULL))

  # cap each signature to its top n_display genes by expression pooled across samples
  pooled_mean <- rowMeans(sapply(avg_list, function(m) rowMeans(m[common_genes, , drop = FALSE])))
  gene_list <- lapply(gene_list, function(g) {
    if (length(g) <= n_display) return(g)
    g[order(pooled_mean[g], decreasing = TRUE)][seq_len(n_display)]
  })

  all_genes    <- unique(unlist(gene_list, use.names = FALSE))
  gene_sig_map <- unlist(unname(mapply(
    function(genes, nm) setNames(rep(nm, length(genes)), genes),
    gene_list, names(gene_list), SIMPLIFY = FALSE)))
  gene_sig_map <- gene_sig_map[!duplicated(names(gene_sig_map))][all_genes]

  # one column per sample x group; z-score each gene across all such columns
  mat_list <- lapply(names(avg_list), function(s) {
    m <- avg_list[[s]][all_genes, , drop = FALSE]
    colnames(m) <- paste0(s, "__", colnames(m))
    m
  })
  mat <- do.call(cbind, mat_list)
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0

  sample_split <- sub("__.*$", "", colnames(mat))
  group_split  <- sub("^.*__", "", colnames(mat))

  sig_colors <- setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(length(gene_list)),
    names(gene_list))
  group_lvls   <- sort(unique(group_split))
  group_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(group_lvls)), group_lvls)
  sample_lvls   <- sort(unique(sample_split))
  sample_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Dark2"))(length(sample_lvls)), sample_lvls)
  col_z <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  gene_ha <- rowAnnotation(
    Signature = gene_sig_map,
    col       = list(Signature = sig_colors),
    show_annotation_name = FALSE)
  top_ha <- HeatmapAnnotation(
    Sample = sample_split,
    Group  = group_split,
    col    = list(Sample = sample_colors, Group = group_colors),
    show_annotation_name = TRUE)

  png(outfile, width = 1800, height = 900, res = 150)
  draw(Heatmap(
    mat,
    name              = "z-score",
    col               = col_z,
    left_annotation   = gene_ha,
    top_annotation    = top_ha,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    row_names_gp      = gpar(fontsize = 7),
    column_title      = "groupgene signature expression (avg per sample x group)",
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_split         = gene_sig_map,
    column_split      = group_split))
  dev.off()
}
Imagescore_meta_programs <- function(obj, sheetname, meta_cols, out_dir,
                                tag = "") {
  require(Seurat)
  require(ggplot2)

  for (sheet in sheetname) {
    # Columns for this sheet: ucell_score named them _meta<SheetName>.<ProgramName>_UCell
    sheet_pat  <- paste0("^", make.names(sheet))
    sheet_cols <- grep(sheet_pat, meta_cols, value = TRUE)
    sheet_cols <- sheet_cols[sheet_cols %in% colnames(obj@meta.data)]
    if (length(sheet_cols) == 0) next

    ncol_g <- min(3L, length(sheet_cols))
    nrow_g <- ceiling(length(sheet_cols) / ncol_g)
    base   <- make.names(sheet)
    plots_s <- lapply(sheet_cols, function(f) {
      v   <- obj@meta.data[[f]]
      mid <- (max(v, na.rm = TRUE) + min(v, na.rm = TRUE)) * 0.3
      tt <- gsub(paste0(sheet, "."), fix = T, "", f)
      tt <- gsub(".", " ",fix = T, tt)
      tt <- gsub("_meta_UCell", "", tt)
      ImageFeaturePlot(obj, features = f) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = mid)+
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                              midpoint = mid) + ggtitle(tt)
    })
    g <- patchwork::wrap_plots(plots_s, ncol = ncol_g)
    ggsave(file.path(out_dir, paste0("meta_", base, tag, ".png")),
           plot = g, width = ncol_g * 8, height = nrow_g * 7,
           dpi = 300, limitsize = FALSE)
  }
}


avg_expr_list <- list()
for (arg in 1:8){
  path  <- paths[arg]
  i     <- basename(path)
  setwd(path)
  srt            <- qs_read("tumour_normal_anno_srt.qs2")
  spatial_dir  <- file.path(path, "final_png", "spatial")
  srt <- AddModuleScore_UCell(srt, features = list(heat.shock_mod), name = "_heat_shock_mod", slot = "data", assay = "SpaNorm")
  heatmap_dir  <- file.path(path, "final_png", "heatmap")
  dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)
  srt <- AddModuleScore_UCell(srt, features = epi_markers, name = "_epi", slot = "data", assay = "SpaNorm")
  avg_expr_list[[i]] <- signature_group_avg(srt, groupgene)

  arch_mod_cols <- grep("_arch", colnames(srt@meta.data), value = TRUE)
  mod_score_cols<- grep("_gd", colnames(srt@meta.data), value = TRUE)
  gs23_cols <-  grep("_gs23", colnames(srt@meta.data), value = TRUE)
  meta_cols <-  grep("_meta", colnames(srt@meta.data), value = TRUE)
  tn_cols <- c("tumour_score_UCell", "normal_score_UCell")

  tumour_cells <- colnames(srt)[srt$compartment == "Tumour"]

  ifp_save(srt, feats = arch_mod_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_archetype_modules.png"), cells_highlight = tumour_cells)
  ifp_save(srt, feats = mod_score_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_G123_scores.png"), cells_highlight = tumour_cells)
  ifp_save(srt, feats = gs23_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_genesets2023.png"))
  ifp_save(srt, feats = tn_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_tumour_normal_scores.png"))
  Imagescore_meta_programs(srt, sheetname = sheetname, meta_cols = meta_cols, out_dir = spatial_dir)

  avg_expr_list[[i]] <- signature_group_avg(srt, groupgene)

  ifp_save_genes(srt, ar_genes, outfile = file.path(spatial_dir, "ImageFeaturePlot_AR_genes.png"))
  ifp_save_genes(srt, ne_genes, outfile = file.path(spatial_dir, "ImageFeaturePlot_NE_genes.png"))

  ifp_save(srt, feats = paste0(names(epi_markers), "_epi"),
           outfile = file.path(spatial_dir, "ImageFeaturePlot_epi_markers.png"))

  ggsave(file.path(spatial_dir, "ImageFeaturePlot_heat_shock_mod.png"),
         plot = ImageFeaturePlot(srt, features = "signature_1_heat_shock_mod") +
           scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                midpoint = 0.3) +
           scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                 midpoint = 0.3),
         width = 8, height = 7, dpi = 300, limitsize = FALSE)

  aggregate_signature_heatmap(avg_expr_list[i], groupgene,
                              outfile = file.path(heatmap_dir, "groupgene_signature_heatmap.png"))

   ifp_save(srt, feats = paste0(names(epi_markers), "_epi"),
           outfile = file.path(spatial_dir, "ImageFeaturePlot_epi_markers.png"))
  ifp_save_genes(srt, ar_genes, outfile = file.path(spatial_dir, "ImageFeaturePlot_AR_genes.png"))
  ifp_save_genes(srt, ne_genes, outfile = file.path(spatial_dir, "ImageFeaturePlot_NE_genes.png"))
}


agg_heatmap_dir <- "~/VisHD/9.2.aggregate_heatmap"
dir.create(agg_heatmap_dir, showWarnings = FALSE, recursive = TRUE)
aggregate_signature_heatmap(avg_expr_list, groupgene,
                            outfile = file.path(agg_heatmap_dir, "groupgene_signature_heatmap_aggregated.png"))