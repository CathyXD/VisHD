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
  library(tidyr)
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

# ── Output roots ─────────────────────────────────────────────────────────────
# Per-sample plots/metas: 9.2.per_sample_visual/<sample>/  (point 1)
# Cross-sample aggregates: 9.2.tumour_normal_visual/       (point 4, "another folder")
per_sample_root <- "~/VisHD/9.2.per_sample_visual"
agg_out_dir      <- "~/VisHD/9.2.tumour_normal_visual"
dir.create(per_sample_root, showWarnings = FALSE, recursive = TRUE)

# ── Single-plot-per-file helpers (point 2: split every combined/wrapped figure) ─
strip_suffix <- function(x) sub("_UCell$", "", x)

ifp_save <- function(srt, feats, outdir, dpi = 300, cells_highlight = NULL) {
  feats <- feats[feats %in% colnames(srt@meta.data)]
  if (length(feats) == 0) return(invisible(NULL))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (f in feats) {
    p <- dark_feature_plot(srt, features = f, ncol = 1, cells_highlight = cells_highlight)
    ggsave(file.path(outdir, paste0(strip_suffix(f), ".png")), p,
           width = 8, height = 7, dpi = dpi, limitsize = FALSE)
  }
}
# gene-expression version of ifp_save: features are pulled from an assay
# (not meta.data). One file per gene, one dark_feature_plot panel each.
ifp_save_genes <- function(srt, genes, outdir, assay = "SpaNorm",
                           dpi = 300, cells_highlight = NULL) {
  genes <- genes[genes %in% rownames(GetAssayData(srt, assay = assay, layer = "data"))]
  if (length(genes) == 0) return(invisible(NULL))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (g in genes) {
    p <- dark_feature_plot(srt, features = g, assay = assay, ncol = 1, cells_highlight = cells_highlight)
    ggsave(file.path(outdir, paste0(g, ".png")), p,
           width = 8, height = 7, dpi = dpi, limitsize = FALSE)
  }
}

# Key gene sets for spatial expression (one file each per gene)
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
groupgene <- readRDS("~/VisHD/6.3.archetype_module/groupgene.Rds")

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

# One ImageFeaturePlot file per meta-program column (point 2: split, no wrap_plots)
Imagescore_meta_programs <- function(obj, sheetname, meta_cols, out_dir) {
  require(Seurat)
  require(ggplot2)

  for (sheet in sheetname) {
    # Columns for this sheet: ucell_score named them _meta<SheetName>.<ProgramName>_UCell
    sheet_pat  <- paste0("^", make.names(sheet))
    sheet_cols <- grep(sheet_pat, meta_cols, value = TRUE)
    sheet_cols <- sheet_cols[sheet_cols %in% colnames(obj@meta.data)]
    if (length(sheet_cols) == 0) next

    sheet_dir <- file.path(out_dir, make.names(sheet))
    dir.create(sheet_dir, showWarnings = FALSE, recursive = TRUE)

    for (f in sheet_cols) {
      v   <- obj@meta.data[[f]]
      mid <- (max(v, na.rm = TRUE) + min(v, na.rm = TRUE)) * 0.3
      tt <- gsub(paste0(sheet, "."), fix = T, "", f)
      tt <- gsub(".", " ",fix = T, tt)
      tt <- gsub("_meta_UCell", "", tt)
      p <- ImageFeaturePlot(obj, features = f) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = mid) +
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                              midpoint = mid) +
        ggtitle(tt)
      ggsave(file.path(sheet_dir, paste0(make.names(tt), ".png")),
             plot = p, width = 8, height = 7, dpi = 300, limitsize = FALSE)
    }
  }
}


# ── Aggregated cross-sample plots (facet-by-slide), mirroring 4.6.tumour_visual.r
# but on tumour_normal_anno_srt.qs2 — compartment/module_anno/category/subclone
# on "umap"/"pearsonumap" + tissue centroids, and the 6 archetype module
# (_arch_UCell) scores as facet_grid(module ~ slide). Written at the end of
# this script into agg_out_dir (point 4, "another folder").
# module_anno palette reused verbatim from 9.3.aggreate_cell_composition_analysis.R
# so colours match every other cross-sample module_anno figure.
group_pal <- c("Neg" = "lightblue", "G1" = "red", "G2" = "gold", "G3" = "royalblue",
               "G1/G2" = "orange", "G1/G3" = "purple", "G2/G3" = "green",
               "G1/G2/G3" = "grey")
canon <- function(x) vapply(strsplit(x, "/"),
                            function(p) paste(sort(p), collapse = "/"), character(1))
build_mg_pal <- function(levs) {
  levs <- unique(as.character(levs))
  pal  <- setNames(group_pal[canon(levs)], levs)
  pal["Normal"] <- "lightpink"
  pal
}
get_centroids <- function(srt) {
  coords <- tryCatch(GetTissueCoordinates(srt, which = "centroids"), error = function(e) NULL)
  if (is.null(coords)) return(data.frame(x_centroid = NA_real_, y_centroid = NA_real_))
  idx <- match(colnames(srt), coords$cell)
  data.frame(x_centroid = coords$x[idx], y_centroid = coords$y[idx])
}
subclone_pal <- c("Normal" = "grey", "1" = "red", "2" = "royalblue")

agg_umap_dir <- file.path(agg_out_dir, "umap")
agg_spa_dir  <- file.path(agg_out_dir, "spatial")
agg_mod_dir  <- file.path(agg_out_dir, "modules")
agg_deg_dir  <- file.path(agg_out_dir, "deg")
agg_heat_dir <- file.path(agg_out_dir, "heatmap")
agg_tbl_dir  <- file.path(agg_out_dir, "combined_tables")
for (d in c(agg_umap_dir, agg_spa_dir, agg_mod_dir, agg_deg_dir, agg_heat_dir, agg_tbl_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

agg_long_list    <- list()
agg_deg_list     <- list()
agg_idp_cat_list <- list()   # per-sample ImageDimPlot(category), combined into one grid later
agg_idp_sub_list <- list()   # per-sample ImageDimPlot(subclone), combined into one grid later

avg_expr_list <- list()
for (arg in 1:8){
  path  <- paths[arg]
  i     <- basename(path)
  setwd(path)
  srt <- qs_read("tumour_normal_anno_srt.qs2")
  srt <- AddModuleScore_UCell(srt, features = list(heat.shock_mod), name = "_heat_shock_mod", slot = "data", assay = "SpaNorm")
  srt <- AddModuleScore_UCell(srt, features = epi_markers, name = "_epi", slot = "data", assay = "SpaNorm")

  # ── Per-sample output folder (point 1: sample-per-folder under 9.2) ────────
  sample_dir  <- file.path(per_sample_root, i)
  spatial_dir <- file.path(sample_dir, "spatial")
  dim_dir     <- file.path(sample_dir, "dimplot")
  heatmap_dir <- file.path(sample_dir, "heatmap")
  for (d in c(spatial_dir, dim_dir, heatmap_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  avg_expr_list[[i]] <- signature_group_avg(srt, groupgene)

  arch_mod_cols <- grep("_arch", colnames(srt@meta.data), value = TRUE)
  mod_score_cols<- grep("_gd", colnames(srt@meta.data), value = TRUE)
  gs23_cols <-  grep("_gs23", colnames(srt@meta.data), value = TRUE)
  meta_cols <-  grep("_meta", colnames(srt@meta.data), value = TRUE)
  tn_cols <- c("tumour_score_UCell", "normal_score_UCell")

  tumour_cells <- colnames(srt)[srt$compartment == "Tumour"]

  # Per-cell embeddings/centroids, reused for this sample's Dim/ImageDim plots,
  # the cross-sample aggregate table, and the full per-cell metadata export.
  have_reduc <- function(r) r %in% Reductions(srt)
  emb <- function(r) if (have_reduc(r)) as.data.frame(Embeddings(srt, r)[, 1:2]) else
                      data.frame(V1 = NA_real_, V2 = NA_real_)
  emb_umap <- emb("umap"); emb_pumap <- emb("pearsonumap")
  centroids <- get_centroids(srt)

  agg_df <- data.frame(
    slide       = i,
    compartment = as.character(srt$compartment),
    module_anno = as.character(srt$module_anno),
    category    = as.character(srt$category),
    subclone    = as.character(srt$subclone)
  )
  agg_df[c("umap1", "umap2")]           <- emb_umap
  agg_df[c("pumap1", "pumap2")]         <- emb_pumap
  agg_df[c("x_centroid", "y_centroid")] <- centroids
  for (mc in arch_mod_cols) agg_df[[sub("_arch_UCell$", "", mc)]] <- srt@meta.data[[mc]]
  agg_long_list[[i]] <- agg_df

  # ── Point 5: save ALL per-cell metadata + cell coordinates for this sample ──
  full_meta <- srt@meta.data
  full_meta$cell       <- rownames(full_meta)
  full_meta$x_centroid <- centroids$x_centroid
  full_meta$y_centroid <- centroids$y_centroid
  full_meta$umap1      <- emb_umap[, 1]
  full_meta$umap2      <- emb_umap[, 2]
  full_meta$pumap1     <- emb_pumap[, 1]
  full_meta$pumap2     <- emb_pumap[, 2]
  saveRDS(full_meta, file.path(sample_dir, "meta.Rds"))
  write.csv(full_meta, file.path(sample_dir, "meta.csv"), row.names = FALSE)

  # ── Point 3: category + subclone, in both DimPlot (UMAP) and ImageDimPlot ──
  ggsave(file.path(dim_dir, "DimPlot_umap_category.png"),
         DimPlot(srt, reduction = "umap", group.by = "category", raster = FALSE) +
           ggtitle(paste(i, "— category")),
         width = 7, height = 6, dpi = 200)
  ggsave(file.path(dim_dir, "DimPlot_umap_subclone.png"),
         DimPlot(srt, reduction = "umap", group.by = "subclone", cols = subclone_pal, raster = FALSE) +
           ggtitle(paste(i, "— subclone")),
         width = 7, height = 6, dpi = 200)
  ggsave(file.path(dim_dir, "DimPlot_pearsonumap_category.png"),
         DimPlot(srt, reduction = "pearsonumap", group.by = "category", raster = FALSE) +
           ggtitle(paste(i, "— category")),
         width = 7, height = 6, dpi = 200)
  ggsave(file.path(dim_dir, "DimPlot_pearsonumap_subclone.png"),
         DimPlot(srt, reduction = "pearsonumap", group.by = "subclone", cols = subclone_pal, raster = FALSE) +
           ggtitle(paste(i, "— subclone")),
         width = 7, height = 6, dpi = 200)

  idp_cat <- ImageDimPlot(srt, group.by = "category", border.color = "#00000000", size = 0.3) +
    ggtitle(paste(i, "— category"))
  idp_sub <- ImageDimPlot(srt, group.by = "subclone", cols = subclone_pal,
                          border.color = "#00000000", size = 0.3) +
    ggtitle(paste(i, "— subclone"))
  ggsave(file.path(spatial_dir, "ImageDimPlot_category.png"), idp_cat, width = 7, height = 6, dpi = 200)
  ggsave(file.path(spatial_dir, "ImageDimPlot_subclone.png"), idp_sub, width = 7, height = 6, dpi = 200)
  agg_idp_cat_list[[i]] <- idp_cat
  agg_idp_sub_list[[i]] <- idp_sub

  if (file.exists("tumour_normal_deg_spanorm.Rds")) {
    deg_top <- readRDS("tumour_normal_deg_spanorm.Rds") %>%
      dplyr::filter(abs(pct.1 - pct.2) > 0.2, p_val_adj < 0.05) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = avg_log2FC, n = 5) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(slide = i, gene, cluster = as.character(cluster), avg_log2FC)
    if (nrow(deg_top) > 0) agg_deg_list[[i]] <- deg_top
  }

  # ── Point 2: every gene-set panel split into one file per feature ──────────
  ifp_save(srt, feats = arch_mod_cols,  outdir = file.path(spatial_dir, "archetype_modules"),   cells_highlight = tumour_cells)
  ifp_save(srt, feats = mod_score_cols, outdir = file.path(spatial_dir, "G123_scores"),         cells_highlight = tumour_cells)
  ifp_save(srt, feats = gs23_cols,      outdir = file.path(spatial_dir, "genesets2023"))
  ifp_save(srt, feats = tn_cols,        outdir = file.path(spatial_dir, "tumour_normal_scores"))
  Imagescore_meta_programs(srt, sheetname = sheetname, meta_cols = meta_cols, out_dir = file.path(spatial_dir, "meta_programs"))

  ifp_save_genes(srt, ar_genes, outdir = file.path(spatial_dir, "AR_genes"))
  ifp_save_genes(srt, ne_genes, outdir = file.path(spatial_dir, "NE_genes"))
  ifp_save(srt, feats = paste0(names(epi_markers), "_epi"), outdir = file.path(spatial_dir, "epi_markers"))

  ggsave(file.path(spatial_dir, "heat_shock_mod.png"),
         plot = ImageFeaturePlot(srt, features = "signature_1_heat_shock_mod") +
           scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                midpoint = 0.3) +
           scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                 midpoint = 0.3),
         width = 8, height = 7, dpi = 300, limitsize = FALSE)

  aggregate_signature_heatmap(avg_expr_list[i], groupgene,
                              outfile = file.path(heatmap_dir, "groupgene_signature_heatmap.png"))

  rm(srt); gc()
}

aggregate_signature_heatmap(avg_expr_list, groupgene,
                            outfile = file.path(agg_heat_dir, "groupgene_signature_heatmap_aggregated.png"))

# ════════════════════════════════════════════════════════════════════════════════
# Aggregated cross-sample plots (facet-by-slide), mirroring 4.6.tumour_visual.r
# Written under agg_out_dir (point 4). Every previously-stacked pair is split
# into its own single-variable file (point 2).
# ════════════════════════════════════════════════════════════════════════════════

agg_big_df <- do.call(rbind, agg_long_list)
slide_levels <- names(agg_long_list)
agg_big_df$slide       <- factor(agg_big_df$slide, levels = slide_levels)
agg_big_df$compartment <- factor(agg_big_df$compartment)
agg_big_df$module_anno <- factor(agg_big_df$module_anno)
agg_big_df$category    <- factor(agg_big_df$category)
agg_big_df$subclone    <- factor(agg_big_df$subclone)
write.csv(agg_big_df, file.path(agg_tbl_dir, "per_cell_long.csv"), row.names = FALSE)

n_slides <- length(slide_levels)
ncf      <- min(4, n_slides)
nr       <- ceiling(n_slides / ncf)
mg_pal <- build_mg_pal(levels(agg_big_df$module_anno))

# One file per (reduction x colour variable) — no patchwork stacking.
save_agg_scatter <- function(xcol, ycol, colour_var, pal, title, file) {
  d <- agg_big_df[!is.na(agg_big_df[[xcol]]), ]
  p <- ggplot(d, aes(.data[[xcol]], .data[[ycol]], colour = .data[[colour_var]])) +
    geom_point(size = 0.3) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    labs(title = title, x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 8)
  if (!is.null(pal)) p <- p + scale_colour_manual(values = pal, name = colour_var)
  ggsave(file, p, width = ncf * 3.2, height = nr * 3.2, dpi = 300, limitsize = FALSE)
}

# SpaNorm UMAP — compartment, module_anno, category, subclone
save_agg_scatter("umap1", "umap2", "compartment", NULL, "SpaNorm UMAP — compartment",
                 file.path(agg_umap_dir, "spanorm_compartment_by_slide.png"))
save_agg_scatter("umap1", "umap2", "module_anno", mg_pal, "SpaNorm UMAP — module group",
                 file.path(agg_umap_dir, "spanorm_module_by_slide.png"))
save_agg_scatter("umap1", "umap2", "category", NULL, "SpaNorm UMAP — category",
                 file.path(agg_umap_dir, "spanorm_category_by_slide.png"))
save_agg_scatter("umap1", "umap2", "subclone", subclone_pal, "SpaNorm UMAP — subclone",
                 file.path(agg_umap_dir, "spanorm_subclone_by_slide.png"))

# Pearson UMAP — compartment, module_anno, category, subclone
save_agg_scatter("pumap1", "pumap2", "compartment", NULL, "Pearson UMAP — compartment",
                 file.path(agg_umap_dir, "pearson_compartment_by_slide.png"))
save_agg_scatter("pumap1", "pumap2", "module_anno", mg_pal, "Pearson UMAP — module group",
                 file.path(agg_umap_dir, "pearson_module_by_slide.png"))
save_agg_scatter("pumap1", "pumap2", "category", NULL, "Pearson UMAP — category",
                 file.path(agg_umap_dir, "pearson_category_by_slide.png"))
save_agg_scatter("pumap1", "pumap2", "subclone", subclone_pal, "Pearson UMAP — subclone",
                 file.path(agg_umap_dir, "pearson_subclone_by_slide.png"))

# ── Spatial (tissue centroid) compartment/module group, faceted by slide ───
d_sp <- agg_big_df[!is.na(agg_big_df$x_centroid), ]
if (nrow(d_sp) > 0) {
  save_agg_scatter_sp <- function(colour_var, pal, title, file) {
    p <- ggplot(d_sp, aes(x_centroid, y_centroid, colour = .data[[colour_var]])) +
      geom_point(size = 0.2) +
      facet_wrap(~ slide, ncol = ncf, scales = "free") +
      coord_fixed() + scale_y_reverse() +
      labs(title = title) +
      theme_bw(base_size = 8)
    if (!is.null(pal)) p <- p + scale_colour_manual(values = pal, name = colour_var)
    ggsave(file, p, width = ncf * 3.2, height = nr * 3.2, dpi = 300, limitsize = FALSE)
  }
  save_agg_scatter_sp("compartment", NULL, "Tissue centroids — compartment",
                      file.path(agg_spa_dir, "spatial_compartment_by_slide.png"))
  save_agg_scatter_sp("module_anno", mg_pal, "Tissue centroids — module group",
                      file.path(agg_spa_dir, "spatial_module_by_slide.png"))
} else {
  message("No tissue centroids recovered for any slide — skipping spatial plot")
}

# ── ImageDimPlot category/subclone, one grid per variable, one panel per slide
# (genuine tissue image, not the centroid scatter above) ───────────────────
ggsave(file.path(agg_spa_dir, "ImageDimPlot_category_by_slide.png"),
       patchwork::wrap_plots(agg_idp_cat_list, ncol = ncf),
       width = ncf * 4.5, height = nr * 3.5, dpi = 200, limitsize = FALSE)
ggsave(file.path(agg_spa_dir, "ImageDimPlot_subclone_by_slide.png"),
       patchwork::wrap_plots(agg_idp_sub_list, ncol = ncf),
       width = ncf * 4.5, height = nr * 3.5, dpi = 200, limitsize = FALSE)

# ── Archetype module scores, facet_grid(module ~ slide) ─────────────────────
module_names <- setdiff(names(agg_long_list[[1]]),
  c("slide","compartment","module_anno","category","subclone",
    "umap1","umap2","pumap1","pumap2","x_centroid","y_centroid"))
if (length(module_names) > 0) {
  make_module_grid <- function(xcol, ycol, title) {
    d <- agg_big_df[!is.na(agg_big_df[[xcol]]), c("slide", xcol, ycol, module_names)]
    d_long <- d %>%
      pivot_longer(cols = all_of(module_names), names_to = "module", values_to = "score") %>%
      mutate(module = factor(module, levels = module_names))
    ggplot(d_long, aes(.data[[xcol]], .data[[ycol]], colour = score)) +
      geom_point(size = 0.25) +
      facet_grid(module ~ slide, scales = "free") +
      scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = 0, name = "score") +
      labs(title = title, x = "UMAP_1", y = "UMAP_2") +
      theme_bw(base_size = 7)
  }
  ggsave(file.path(agg_mod_dir, "archetype_module_umap_by_slide.png"),
         make_module_grid("umap1", "umap2", "Archetype modules — SpaNorm UMAP"),
         width = n_slides * 2.2, height = length(module_names) * 2, dpi = 300, limitsize = FALSE)
  ggsave(file.path(agg_mod_dir, "archetype_module_pearsonumap_by_slide.png"),
         make_module_grid("pumap1", "pumap2", "Archetype modules — Pearson UMAP"),
         width = n_slides * 2.2, height = length(module_names) * 2, dpi = 300, limitsize = FALSE)
}

# ── Top-DEG heatmap across slides (BANKSY clusters, per-sample) ────────────
if (length(agg_deg_list) == 0) {
  message("No tumour_normal_deg_spanorm.Rds found for any slide — skipping DEG heatmap")
} else {
  agg_deg_df <- do.call(rbind, agg_deg_list)
  write.csv(agg_deg_df, file.path(agg_tbl_dir, "top_deg_per_slide.csv"), row.names = FALSE)

  genes <- agg_deg_df %>% dplyr::count(gene) %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(gene)
  deg_mat <- matrix(NA_real_, nrow = length(genes), ncol = n_slides,
                    dimnames = list(genes, slide_levels))
  deg_agg <- agg_deg_df %>% dplyr::group_by(gene, slide) %>%
    dplyr::summarise(avg_log2FC = max(avg_log2FC), .groups = "drop")
  deg_mat[cbind(deg_agg$gene, as.character(deg_agg$slide))] <- deg_agg$avg_log2FC

  ht <- ComplexHeatmap::Heatmap(
    deg_mat, name = "avg_log2FC", na_col = "grey95",
    cluster_columns = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_title = "Top per-sample DEGs (own BANKSY clustering) x slide"
  )
  png(file.path(agg_deg_dir, "top_deg_by_slide_heatmap.png"),
      width = max(6, n_slides * 0.6 + 3), height = max(5, length(genes) * 0.15 + 2),
      units = "in", res = 200)
  ComplexHeatmap::draw(ht)
  dev.off()
}

cat("\n==================== 9.2 aggregated cross-sample plots done ====================\n")
cat("Per-sample outputs under", per_sample_root, "\n")
cat("Aggregated outputs under", agg_out_dir, "\n")
