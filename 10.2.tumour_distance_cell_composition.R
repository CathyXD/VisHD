#!/usr/bin/env Rscript
# 10.2.tumour_distance_cell_composition.R   (run-once, all 8 samples)
# Radius-based neighbour-composition analysis on the per-sample
# LUT-245-XX/tumour_normal_anno_srt.qs2 objects. Neighbour counts are
# inclusive (a tumour cell's own module group / tumour cells in general are
# not excluded from the composition — they are just labelled coarsely).
#
#   Section 1 — for EVERY cell (tumour + normal), distance (µm) to its
#               nearest TUMOUR cell (of any subclone/module group); at each
#               radius r, tabulate the composition of ALL cells within r of
#               any tumour cell. Tumour neighbours are collapsed to a single
#               "Tumour" identity; normal neighbours keep their
#               final_annotation cell-type name.
#   Section 2 — same idea, but the reference population is each tumour
#               module group in turn (G1 / G2 / G3, inclusive of combo
#               labels e.g. "G1/G2"). The composition of ALL cells within r
#               (including the focal group's own cells) is tabulated, with
#               neighbours labelled by normal final_annotation cell type OR,
#               for tumour cells, their exact module_anno value ("Neg",
#               "G1", ..., "G1/G2/G3" — no collapsing).
#
# Both sections get the same per-cell visualization: for every cell in the
# reference/self population, its own neighbour composition at each radius is
# saved as a table (columns = individual self-population cells) and rendered
# as a column-scaled heatmap. ALL tumour-identity rows ("Tumour" in section 1;
# every module_anno value in section 2) are pulled out of the heatmap body —
# so the body holds only normal cell types — and shown instead as a stacked
# top annotation (alongside pearson_clusters); columns are ranked by the
# self population's own proportion.
#
# Coordinates come from GetTissueCoordinates(srt, which="centroids") and are
# in FOV pixel units (no Space Ranger scale-factor file exists for this
# GeoJSON-segmentation pipeline). Converted to µm using the same ~0.29 µm/px
# approximation already used in 10.1.Statial.R.
#
# Outputs under ~/VisHD/10.2.tumour_distance_cell_composition/:
#   section1_any_tumour/    per_sample_tables/ per_sample_barplots/
#                           combined_tables/   aggregate_barplots/
#                           per_cell_tables/    <slide>_r<radius>.csv —
#                             columns = individual tumour cells, rows =
#                             "Tumour" + normal final_annotation
#                           heatmaps/            <slide>_r<radius>.png
#   section2_module_group/  per_sample_tables/ per_sample_barplots/
#                           combined_tables/   aggregate_barplots/
#                           per_cell_tables/     <slide>_<group>_r<radius>.csv —
#                             columns = individual cells of that module group,
#                             rows = normal final_annotation + exact module_anno
#                           heatmaps/             <slide>_<group>_r<radius>.png
#
#   Rscript 10.2.tumour_distance_cell_composition.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(pals)
  library(qs2)
  library(FNN)
  library(dbscan)
  library(ComplexHeatmap)
  library(circlize)
})

setwd("~/VisHD")

UM_PER_PX     <- 0.29                      # same approximation as 10.1.Statial.R
RADII_UM      <- c(10, 20, 30, 50, 100)
RADII_PX      <- RADII_UM / UM_PER_PX
MODULE_GROUPS <- c("G1", "G2", "G3")

out_dir <- "~/VisHD/10.2.tumour_distance_cell_composition"
s1_dir  <- file.path(out_dir, "section1_any_tumour")
s2_dir  <- file.path(out_dir, "section2_module_group")
for (d in c(s1_dir, s2_dir))
  for (sd in c("per_sample_tables", "per_sample_barplots", "combined_tables",
               "aggregate_barplots", "per_cell_tables", "heatmaps"))
    dir.create(file.path(d, sd), recursive = TRUE, showWarnings = FALSE)

# ── Palette helper ─────────────────────────────────────────────────────────
# module_anno colours reused verbatim from 9.1/9.3, covering every value the
# column can take (Neg, the 3 pure groups, and all combos).
group_pal <- c("Neg" = "lightblue", "G1" = "red", "G2" = "gold", "G3" = "royalblue",
              "G2/G1" = "orange", "G3/G1" = "purple", "G3/G2" = "green",
              "G3/G2/G1" = "grey")

# Normal cell-type names get a Set3 palette; anything named in `special_pal`
# (e.g. Tumour=red, or G1/G2/G3/Other) keeps its fixed colour instead.
build_composition_pal <- function(levs, special_pal) {
  levs <- unique(as.character(levs))
  oth  <- setdiff(levs, names(special_pal))
  pal  <- if (length(oth))
    setNames(as.vector(pals::brewer.set3(max(3, length(oth))))[seq_along(oth)], oth)
  else character(0)
  matched <- intersect(names(special_pal), levs)
  if (length(matched)) pal[matched] <- special_pal[matched]
  pal
}

get_centroids <- function(srt) {
  coords <- GetTissueCoordinates(srt, which = "centroids")
  idx    <- match(colnames(srt), coords$cell)
  stopifnot(!anyNA(idx))
  as.matrix(coords[idx, c("x", "y")])
}

# ── Per-cell composition table + heatmap, shared by section 1 and section 2 ──
# For every cell in `self_idx` (the reference/self population defining the
# heatmap's columns), counts ALL cells (including itself) within each radius,
# labelled by `identity_vec`. Every value in `tumour_levels` is pulled out of
# the heatmap body (which then holds only normal cell types) and shown
# instead as a stacked top annotation named `tumour_anno_name`, coloured by
# `tumour_pal`. `self_label` (one of `tumour_levels`) is used to rank columns
# by that self population's own proportion. A second top annotation,
# `ref_module_anno` (coloured by `ref_module_pal`), shows each reference/self
# cell's own `module_anno` identity (from `module_anno_all`) as a single
# categorical strip — distinct from the neighbour-composition barplot above.
build_percell_heatmaps <- function(xy, cell_ids, self_idx, identity_vec, self_label,
                                   tumour_levels, tumour_pal, tumour_anno_name,
                                   module_anno_all, ref_module_pal,
                                   pearson_clusters_all, pc_pal, out_dir, tag, sample_name) {
  self_xy  <- xy[self_idx, , drop = FALSE]
  self_ids <- cell_ids[self_idx]
  if (length(self_ids) == 0) {
    cat("  ", tag, "skipped: no self cells\n")
    return(invisible(NULL))
  }

  fa_levels       <- sort(unique(identity_vec))
  frnn            <- dbscan::frNN(x = xy, query = self_xy, eps = max(RADII_PX), sort = FALSE)
  pearson_self    <- as.character(pearson_clusters_all[self_idx])
  ref_module_self <- as.character(module_anno_all[self_idx])

  for (ridx in seq_along(RADII_UM)) {
    r_px <- RADII_PX[ridx]
    mat  <- matrix(0L, nrow = length(fa_levels), ncol = length(self_ids),
                   dimnames = list(fa_levels, self_ids))
    for (t in seq_along(self_ids)) {
      ids <- frnn$id[[t]]
      if (length(ids) == 0) next
      keep <- ids[frnn$dist[[t]] <= r_px]
      if (length(keep) == 0) next
      tab <- table(identity_vec[keep])
      mat[names(tab), t] <- as.integer(tab)
    }
    out_tbl <- data.frame(identity = rownames(mat), mat, check.names = FALSE)
    write.csv(out_tbl,
              file.path(out_dir, "per_cell_tables", sprintf("%s_%s_r%d.csv", sample_name, tag, RADII_UM[ridx])),
              row.names = FALSE)

    # Column-scaled (proportion of each cell's own neighbour total).
    mat_prop <- apply(mat, 2, function(col) {
      s <- sum(col)
      if (s > 0) col / s else col
    })
    self_prop  <- if (self_label %in% rownames(mat_prop)) mat_prop[self_label, ] else rep(0, ncol(mat_prop))
    tum_rows   <- intersect(tumour_levels, rownames(mat_prop))
    mat_body   <- mat_prop[setdiff(rownames(mat_prop), tumour_levels), , drop = FALSE]
    mat_tumour <- mat_prop[tum_rows, , drop = FALSE]

    # No clustering — just rank columns by the self population's own proportion.
    ord        <- order(self_prop, decreasing = TRUE)
    mat_body   <- mat_body[, ord, drop = FALSE]
    mat_tumour <- mat_tumour[, ord, drop = FALSE]
    pc_ord     <- pearson_self[ord]
    ref_ord    <- ref_module_self[ord]

    anno_args   <- list(pearson_cluster = pc_ord, ref_module_anno = ref_ord)
    anno_cols   <- list(pearson_cluster = pc_pal, ref_module_anno = ref_module_pal)
    legend_list <- list()
    if (nrow(mat_tumour) > 0) {
      anno_args[[tumour_anno_name]] <- ComplexHeatmap::anno_barplot(
        t(mat_tumour), bar_width = 1, border = FALSE,
        gp = grid::gpar(fill = tumour_pal[rownames(mat_tumour)], col = NA),
        axis_param = list(gp = grid::gpar(fontsize = 6))
      )
      legend_list <- list(ComplexHeatmap::Legend(
        labels = rownames(mat_tumour), title = tumour_anno_name,
        legend_gp = grid::gpar(fill = tumour_pal[rownames(mat_tumour)])
      ))
    }
    top_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
                        c(anno_args, list(col = anno_cols, annotation_name_gp = grid::gpar(fontsize = 8))))

    ht <- ComplexHeatmap::Heatmap(
      mat_body,
      name              = "proportion",
      col               = circlize::colorRamp2(c(0, 0.5, 1), c("white", "orange", "indianred")),
      cluster_rows      = TRUE,
      cluster_columns   = FALSE,
      show_column_names = FALSE,
      row_names_gp      = grid::gpar(fontsize = 8),
      top_annotation    = top_anno,
      column_title      = sprintf("%s %s — per-cell composition (r = %d µm)", sample_name, tag, RADII_UM[ridx])
    )
    png(file.path(out_dir, "heatmaps", sprintf("%s_%s_r%d.png", sample_name, tag, RADII_UM[ridx])),
        width = 10, height = max(5, nrow(mat_body) * 0.3 + 2) + 2, units = "in", res = 150)
    ComplexHeatmap::draw(ht, annotation_legend_list = legend_list)
    dev.off()
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Per-sample loop
# ══════════════════════════════════════════════════════════════════════════════
paths   <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
samples <- basename(paths)

section1_tables <- list()
section2_tables <- list()

for (k in seq_along(paths)) {
  path <- paths[k]; i <- samples[k]
  cat("==================== ", i, " ====================\n")
  srt      <- qs_read(file.path(path, "tumour_normal_anno_srt.qs2"))
  xy       <- get_centroids(srt)
  cell_ids <- colnames(srt)
  meta     <- srt@meta.data
  is_tumour <- meta$compartment == "Tumour"
  cat("  tumour cells:", sum(is_tumour), " normal cells:", sum(!is_tumour), "\n")

  pearson_clusters_all <- as.character(meta$pearson_clusters)
  pc_levels <- sort(unique(pearson_clusters_all))
  pc_pal    <- setNames(as.vector(pals::polychrome())[seq_along(pc_levels)], pc_levels)
  module_anno_chr <- as.character(meta$module_anno)   # reference-cell identity annotation (both sections)

  # ── Section 1: any-tumour reference, ALL cells counted, tumour collapsed ──
  tumour_xy <- xy[is_tumour, , drop = FALSE]
  identity1 <- ifelse(is_tumour, "Tumour", as.character(meta$final_annotation))

  if (nrow(tumour_xy) == 0) {
    cat("  section 1 skipped: no tumour cells\n")
  } else {
    dist_um <- FNN::get.knnx(data = tumour_xy, query = xy, k = 1)$nn.dist[, 1] * UM_PER_PX
    df1 <- do.call(rbind, lapply(RADII_UM, function(r) {
      near <- dist_um <= r
      if (!any(near)) return(NULL)
      tab <- table(final_annotation = identity1[near])
      data.frame(radius = r, final_annotation = names(tab), Freq = as.numeric(tab),
                stringsAsFactors = FALSE)
    }))
    write.csv(df1, file.path(s1_dir, "per_sample_tables", paste0(i, ".csv")), row.names = FALSE)
    section1_tables[[i]] <- df1

    pal1 <- build_composition_pal(df1$final_annotation, c(Tumour = "red"))
    p1 <- ggplot(df1, aes(factor(radius), Freq, fill = final_annotation)) +
      geom_col(position = "fill") +
      scale_fill_manual(values = pal1, na.value = "grey80", name = "final_annotation") +
      scale_y_continuous(labels = percent) +
      labs(title = paste(i, "— composition around any tumour cell"),
           x = "radius (µm)", y = "proportion") +
      theme_bw(base_size = 10)
    ggsave(file.path(s1_dir, "per_sample_barplots", paste0(i, ".png")),
           p1, width = 8, height = 6, dpi = 150)

    build_percell_heatmaps(xy, cell_ids, is_tumour, identity1, self_label = "Tumour",
                           tumour_levels = "Tumour", tumour_pal = c(Tumour = "red"),
                           tumour_anno_name = "Tumour",
                           module_anno_all = module_anno_chr, ref_module_pal = group_pal,
                           pearson_clusters_all = pearson_clusters_all, pc_pal = pc_pal,
                           out_dir = s1_dir, tag = "Tumour", sample_name = i)
  }

  # ── Section 2: per module-group reference, ALL cells counted, labelled by
  #    normal final_annotation OR (for tumour cells) their exact module_anno ──
  final_annotation_chr <- as.character(meta$final_annotation)
  identity2      <- ifelse(!is_tumour, final_annotation_chr, module_anno_chr)
  tumour_levels2 <- sort(unique(module_anno_chr[is_tumour]))

  df2_parts <- lapply(MODULE_GROUPS, function(g) {
    self_idx <- is_tumour & grepl(g, module_anno_chr)
    if (sum(self_idx) == 0) {
      cat("  section 2 module group", g, "skipped: 0 self cells\n")
      return(NULL)
    }
    self_xy <- xy[self_idx, , drop = FALSE]
    dist_um <- FNN::get.knnx(data = self_xy, query = xy, k = 1)$nn.dist[, 1] * UM_PER_PX
    tbl <- do.call(rbind, lapply(RADII_UM, function(r) {
      near <- dist_um <= r
      if (!any(near)) return(NULL)
      tab <- table(neighbour_group = identity2[near])
      data.frame(module_group = g, radius = r, neighbour_group = names(tab),
                Freq = as.numeric(tab), stringsAsFactors = FALSE)
    }))

    build_percell_heatmaps(xy, cell_ids, self_idx, identity2, self_label = g,
                           tumour_levels = tumour_levels2, tumour_pal = group_pal,
                           tumour_anno_name = "module_anno",
                           module_anno_all = module_anno_chr, ref_module_pal = group_pal,
                           pearson_clusters_all = pearson_clusters_all, pc_pal = pc_pal,
                           out_dir = s2_dir, tag = g, sample_name = i)
    tbl
  })
  df2 <- do.call(rbind, df2_parts)

  if (!is.null(df2)) {
    write.csv(df2, file.path(s2_dir, "per_sample_tables", paste0(i, ".csv")), row.names = FALSE)
    section2_tables[[i]] <- df2

    pal2 <- build_composition_pal(df2$neighbour_group, group_pal)
    p2 <- ggplot(df2, aes(factor(radius), Freq, fill = neighbour_group)) +
      geom_col(position = "fill") +
      facet_wrap(~ module_group) +
      scale_fill_manual(values = pal2, na.value = "grey80", name = "neighbour_group") +
      scale_y_continuous(labels = percent) +
      labs(title = paste(i, "— composition around tumour module groups"),
           x = "radius (µm)", y = "proportion") +
      theme_bw(base_size = 10)
    ggsave(file.path(s2_dir, "per_sample_barplots", paste0(i, ".png")),
           p2, width = 12, height = 5, dpi = 150)
  }

  rm(srt); gc()
}

# ══════════════════════════════════════════════════════════════════════════════
# Aggregation across all 8 samples
# ══════════════════════════════════════════════════════════════════════════════
cat("\n==================== aggregating across samples ====================\n")

# ── Section 1 ──────────────────────────────────────────────────────────────
if (length(section1_tables) > 0) {
  agg1 <- do.call(rbind, Map(function(nm, d) cbind(slide = nm, d),
                             names(section1_tables), section1_tables))
  agg1 <- agg1 %>%
    group_by(slide, radius) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ungroup()
  write.csv(agg1, file.path(s1_dir, "combined_tables", "all_samples_composition.csv"),
            row.names = FALSE)

  pal1_all <- build_composition_pal(agg1$final_annotation, c(Tumour = "red"))
  n_s1 <- length(unique(agg1$slide))

  p_facet1 <- ggplot(agg1, aes(factor(radius), prop, fill = final_annotation)) +
    geom_col() +
    facet_wrap(~ slide, ncol = min(4, n_s1)) +
    scale_fill_manual(values = pal1_all, na.value = "grey80", name = "final_annotation") +
    scale_y_continuous(labels = percent) +
    labs(title = "Composition around any tumour cell, by sample",
         x = "radius (µm)", y = "proportion") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(s1_dir, "aggregate_barplots", "faceted_by_slide.png"),
         p_facet1, width = 14, height = 10, dpi = 150, limitsize = FALSE)

  mean1 <- agg1 %>%
    group_by(radius, final_annotation) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")
  p_mean1 <- ggplot(mean1, aes(factor(radius), prop, fill = final_annotation)) +
    geom_col() +
    scale_fill_manual(values = pal1_all, na.value = "grey80", name = "final_annotation") +
    scale_y_continuous(labels = percent) +
    labs(title = "Mean composition around any tumour cell, across samples",
         x = "radius (µm)", y = "mean proportion") +
    theme_bw(base_size = 10)
  ggsave(file.path(s1_dir, "aggregate_barplots", "mean_across_samples.png"),
         p_mean1, width = 8, height = 6, dpi = 150)
}

# ── Section 2 ──────────────────────────────────────────────────────────────
if (length(section2_tables) > 0) {
  agg2 <- do.call(rbind, Map(function(nm, d) cbind(slide = nm, d),
                             names(section2_tables), section2_tables))
  agg2 <- agg2 %>%
    group_by(slide, module_group, radius) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ungroup()
  write.csv(agg2, file.path(s2_dir, "combined_tables", "all_samples_composition.csv"),
            row.names = FALSE)

  pal2_all <- build_composition_pal(agg2$neighbour_group, group_pal)
  n_s2 <- length(unique(agg2$slide))

  p_facet2 <- ggplot(agg2, aes(factor(radius), prop, fill = neighbour_group)) +
    geom_col() +
    facet_grid(module_group ~ slide) +
    scale_fill_manual(values = pal2_all, na.value = "grey80", name = "neighbour_group") +
    scale_y_continuous(labels = percent) +
    labs(title = "Composition around tumour module groups, by sample",
         x = "radius (µm)", y = "proportion") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(s2_dir, "aggregate_barplots", "faceted_by_slide.png"),
         p_facet2, width = max(14, n_s2 * 2.5), height = 8, dpi = 150, limitsize = FALSE)

  mean2 <- agg2 %>%
    group_by(module_group, radius, neighbour_group) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")
  p_mean2 <- ggplot(mean2, aes(factor(radius), prop, fill = neighbour_group)) +
    geom_col() +
    facet_wrap(~ module_group) +
    scale_fill_manual(values = pal2_all, na.value = "grey80", name = "neighbour_group") +
    scale_y_continuous(labels = percent) +
    labs(title = "Mean composition around tumour module groups, across samples",
         x = "radius (µm)", y = "mean proportion") +
    theme_bw(base_size = 10)
  ggsave(file.path(s2_dir, "aggregate_barplots", "mean_across_samples.png"),
         p_mean2, width = 12, height = 5, dpi = 150)
}

cat("\n==================== 10.2.tumour_distance_cell_composition done ====================\n")
cat("Outputs under", out_dir, "\n")
