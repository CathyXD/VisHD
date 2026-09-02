#!/usr/bin/env Rscript
# 4.6.tumour_visual.r   (run-once, all 8 samples)
# Re-creates the per-sample diagnostic plots from 4.1.tumour_split.R
# (LUT-245-XX/tumour/tumour_srt.qs2) as combined, facet-by-slide figures
# instead of 8 separate per-sample files:
#
#   umap/     category + subclone, faceted by slide, on the SpaNorm "umap"
#             and on "pearsonumap" (recaptures spanorm/pearson_category_subclone.png)
#   spatial/  same two labels, faceted by slide, on tissue centroid coordinates
#             (recaptures the ImageDimPlot half of spanorm_category_subclone.png)
#   modules/  the 6 archetype module scores (AR/Inflammation/NE1/NE2/Cycling/
#             Glycolysis), facet_grid(module ~ slide) on "umap" and
#             "pearsonumap" (recaptures archetype_module_exp.pdf /
#             pearson_archetype_module_exp.pdf)
#   svec/     SVEC marker panel, faceted by slide, on "pearsonumap" and
#             "banksy0.2.umap" (recaptures SVEC_marker.png / pearson_SVEC_marker.png).
#             13 genes x 8 slides is too many panels to show individually, so
#             this uses one composite score (mean SpaNorm expression of the
#             present markers) per cell instead of per-gene FeaturePlots.
#   deg/      cross-sample heatmap of each sample's own top-DEG genes (from
#             deg_spanorm.Rds, if present) x slide, filled by avg_log2FC.
#             Skipped with a message if no sample has written that file yet.
#   combined_tables/   the per-cell long table backing the umap/module/svec plots.
#
# `subclone`/cluster labels are per-sample (own clustering run) and are NOT
# biologically aligned across slides — colours/order are per-panel only.
#
#   Rscript 4.6.tumour_visual.r

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pals)
})

setwd("~/VisHD")

out_dir  <- "~/VisHD/4.6.tumour_visual"
umap_dir <- file.path(out_dir, "umap")
spa_dir  <- file.path(out_dir, "spatial")
mod_dir  <- file.path(out_dir, "modules")
svec_dir <- file.path(out_dir, "svec")
deg_dir  <- file.path(out_dir, "deg")
tbl_dir  <- file.path(out_dir, "combined_tables")
for (d in c(umap_dir, spa_dir, mod_dir, svec_dir, deg_dir, tbl_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

module_names <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")
mod_cols     <- paste0(module_names, "_Module")

SVEC_marker <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                 "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")

# ── Helpers ──────────────────────────────────────────────────────────────────
qual_pal <- function(levs) {
  levs <- as.character(levs)
  n    <- length(levs)
  base <- as.vector(pals::polychrome(36))
  setNames(if (n <= 36) base[seq_len(n)] else grDevices::colorRampPalette(base)(n), levs)
}

# Tissue centroids matched by cell id (barcode ordering silently mismatches — see 4.4).
get_centroids <- function(srt) {
  coords <- tryCatch(GetTissueCoordinates(srt, which = "centroids"),
                     error = function(e) NULL)
  if (is.null(coords)) return(data.frame(x_centroid = NA_real_, y_centroid = NA_real_))
  idx <- match(colnames(srt), coords$cell)
  data.frame(x_centroid = coords$x[idx], y_centroid = coords$y[idx])
}

# ── Discover samples ─────────────────────────────────────────────────────────
paths  <- system("realpath ~/VisHD/LUT-245-*/tumour/tumour_srt.qs2", intern = TRUE)
slides <- basename(dirname(dirname(paths)))
cat("Found", length(paths), "slides:", paste(slides, collapse = ", "), "\n")

slide_levels <- slides
slide_pal    <- qual_pal(slide_levels)

# ── 1. Per-sample load -> combined long table + per-sample top-DEG table ────
long_list <- list()
deg_list  <- list()

for (i in seq_along(paths)) {
  slide <- slides[i]
  cat("Loading", slide, "\n")
  srt <- qs_read(paths[i])

  clone_col <- if ("subclone" %in% colnames(srt@meta.data)) "subclone" else
               if ("ATAClone_cluster" %in% colnames(srt@meta.data)) "ATAClone_cluster" else "category"

  have_reduc <- function(r) r %in% Reductions(srt)
  emb <- function(r) if (have_reduc(r)) as.data.frame(Embeddings(srt, r)[, 1:2]) else
                      data.frame(V1 = NA_real_, V2 = NA_real_)

  df <- data.frame(
    slide    = slide,
    category = as.character(srt@meta.data$category),
    clone    = as.character(srt@meta.data[[clone_col]])
  )
  df[c("umap1", "umap2")]         <- emb("umap")
  df[c("pumap1", "pumap2")]       <- emb("pearsonumap")
  df[c("bumap1", "bumap2")]       <- emb("banksy0.2.umap")
  df[c("x_centroid", "y_centroid")] <- get_centroids(srt)

  missing_mod <- setdiff(mod_cols, colnames(srt@meta.data))
  if (length(missing_mod) > 0)
    message("  ", slide, ": missing module cols ", paste(missing_mod, collapse = ", "))
  for (mc in mod_cols)
    df[[mc]] <- if (mc %in% colnames(srt@meta.data)) srt@meta.data[[mc]] else NA_real_

  DefaultAssay(srt) <- "SpaNorm"
  svec_present <- intersect(SVEC_marker, rownames(srt))
  df$SVEC_score <- if (length(svec_present) > 0) {
    Matrix::colMeans(GetAssayData(srt, assay = "SpaNorm", layer = "data")[svec_present, , drop = FALSE])
  } else NA_real_

  long_list[[slide]] <- df

  deg_path <- file.path(dirname(paths[i]), "deg_spanorm.Rds")
  if (file.exists(deg_path)) {
    deg <- readRDS(deg_path)
    top <- deg %>%
      dplyr::filter(abs(pct.1 - pct.2) > 0.2, p_val_adj < 0.05) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = avg_log2FC, n = 5) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(slide = slide, gene, cluster = as.character(cluster), avg_log2FC, pct.1, pct.2)
    if (nrow(top) > 0) deg_list[[slide]] <- top
  }

  rm(srt); gc()
}

big_df <- do.call(rbind, long_list)
big_df$slide    <- factor(big_df$slide, levels = slide_levels)
big_df$category <- factor(big_df$category)
big_df$clone    <- factor(big_df$subclone)
write.csv(big_df, file.path(tbl_dir, "per_cell_long.csv"), row.names = FALSE)

n_slides <- length(slide_levels)
ncf      <- min(4, n_slides)
category_pal <- qual_pal(levels(big_df$category))
clone_pal    <- qual_pal(levels(big_df$clone))

# ── 2. UMAP category/subclone, faceted by slide ─────────────────────────────
make_umap_pair <- function(xcol, ycol, title_prefix) {
  d <- big_df[!is.na(big_df[[xcol]]), ]
  p_cat <- ggplot(d, aes(.data[[xcol]], .data[[ycol]], colour = category)) +
    geom_point(size = 0.3) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    scale_colour_manual(values = category_pal, name = "category") +
    labs(title = paste(title_prefix, "— category"), x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 8)
  p_clone <- ggplot(d, aes(.data[[xcol]], .data[[ycol]], colour = clone)) +
    geom_point(size = 0.3) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    scale_colour_manual(values = clone_pal, name = "clone") +
    labs(title = paste(title_prefix, "— subclone"), x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 8)
  p_cat / p_clone
}

nr <- ceiling(n_slides / ncf)
ggsave(file.path(umap_dir, "spanorm_category_subclone_by_slide.png"),
       make_umap_pair("umap1", "umap2", "SpaNorm UMAP"),
       width = ncf * 3.2, height = nr * 3.2 * 2, dpi = 300, limitsize = FALSE)

ggsave(file.path(umap_dir, "pearson_category_subclone_by_slide.png"),
       make_umap_pair("pumap1", "pumap2", "Pearson UMAP"),
       width = ncf * 3.2, height = nr * 3.2 * 2, dpi = 300, limitsize = FALSE)

# ── 3. Spatial (tissue centroid) category/subclone, faceted by slide ───────
d_sp <- big_df[!is.na(big_df$x_centroid), ]
if (nrow(d_sp) > 0) {
  p_sp_cat <- ggplot(d_sp, aes(x_centroid, y_centroid, colour = category)) +
    geom_point(size = 0.2) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    scale_colour_manual(values = category_pal, name = "category") +
    coord_fixed() + scale_y_reverse() +
    labs(title = "Tissue centroids — category") +
    theme_bw(base_size = 8)
  p_sp_clone <- ggplot(d_sp, aes(x_centroid, y_centroid, colour = clone)) +
    geom_point(size = 0.2) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    scale_colour_manual(values = clone_pal, name = "clone") +
    coord_fixed() + scale_y_reverse() +
    labs(title = "Tissue centroids — subclone") +
    theme_bw(base_size = 8)
  ggsave(file.path(spa_dir, "spatial_category_subclone_by_slide.png"),
         p_sp_cat / p_sp_clone,
         width = ncf * 3.2, height = nr * 3.2 * 2, dpi = 300, limitsize = FALSE)
} else {
  message("No tissue centroids recovered for any slide — skipping spatial plot")
}

# ── 4. Archetype module scores, facet_grid(module ~ slide) ──────────────────
make_module_grid <- function(xcol, ycol, title) {
  d <- big_df[!is.na(big_df[[xcol]]), c("slide", xcol, ycol, mod_cols)]
  d_long <- d %>%
    pivot_longer(cols = all_of(mod_cols), names_to = "module", values_to = "score") %>%
    mutate(module = factor(sub("_Module$", "", module), levels = module_names))
  ggplot(d_long, aes(.data[[xcol]], .data[[ycol]], colour = score)) +
    geom_point(size = 0.25) +
    facet_grid(module ~ slide, scales = "free") +
    scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                           midpoint = 0, name = "score") +
    labs(title = title, x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 7)
}

ggsave(file.path(mod_dir, "archetype_module_umap_by_slide.png"),
       make_module_grid("umap1", "umap2", "Archetype modules — SpaNorm UMAP"),
       width = n_slides * 2.2, height = length(module_names) * 2, dpi = 300, limitsize = FALSE)

ggsave(file.path(mod_dir, "archetype_module_pearsonumap_by_slide.png"),
       make_module_grid("pumap1", "pumap2", "Archetype modules — Pearson UMAP"),
       width = n_slides * 2.2, height = length(module_names) * 2, dpi = 300, limitsize = FALSE)

# ── 5. SVEC composite score, faceted by slide ───────────────────────────────
make_svec_plot <- function(xcol, ycol, title) {
  d <- big_df[!is.na(big_df[[xcol]]) & !is.na(big_df$SVEC_score), ]
  ggplot(d, aes(.data[[xcol]], .data[[ycol]], colour = SVEC_score)) +
    geom_point(size = 0.3) +
    facet_wrap(~ slide, ncol = ncf, scales = "free") +
    scale_colour_gradient(low = "grey90", high = "darkblue", name = "SVEC score") +
    labs(title = title, x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 8)
}

ggsave(file.path(svec_dir, "SVEC_score_pearsonumap_by_slide.png"),
       make_svec_plot("pumap1", "pumap2", "SVEC composite score — Pearson UMAP"),
       width = ncf * 3.2, height = nr * 3.2, dpi = 300, limitsize = FALSE)

ggsave(file.path(svec_dir, "SVEC_score_banksyumap_by_slide.png"),
       make_svec_plot("bumap1", "bumap2", "SVEC composite score — BANKSY (λ=0.2) UMAP"),
       width = ncf * 3.2, height = nr * 3.2, dpi = 300, limitsize = FALSE)

# ── 6. Top-DEG heatmap across slides ─────────────────────────────────────────
if (length(deg_list) == 0) {
  message("No deg_spanorm.Rds found for any slide — skipping DEG heatmap")
} else {
  deg_df <- do.call(rbind, deg_list)
  write.csv(deg_df, file.path(tbl_dir, "top_deg_per_slide.csv"), row.names = FALSE)

  genes <- deg_df %>% dplyr::count(gene) %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(gene)
  mat <- matrix(NA_real_, nrow = length(genes), ncol = n_slides,
               dimnames = list(genes, slide_levels))
  agg <- deg_df %>% dplyr::group_by(gene, slide) %>%
    dplyr::summarise(avg_log2FC = max(avg_log2FC), .groups = "drop")
  mat[cbind(agg$gene, as.character(agg$slide))] <- agg$avg_log2FC

  library(ComplexHeatmap)
  ht <- ComplexHeatmap::Heatmap(
    mat, name = "avg_log2FC", na_col = "grey95",
    cluster_columns = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_title = "Top per-sample DEGs (own clustering) x slide"
  )
  png(file.path(deg_dir, "top_deg_by_slide_heatmap.png"),
      width = max(6, n_slides * 0.6 + 3), height = max(5, length(genes) * 0.15 + 2),
      units = "in", res = 200)
  ComplexHeatmap::draw(ht)
  dev.off()
}

cat("\n==================== 4.6.tumour_visual done ====================\n")
cat("Outputs under", out_dir, "\n")
