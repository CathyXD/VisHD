#!/usr/bin/env Rscript
# 6.4.merged_module_analysis.R   (run-once, all 8 samples)
# Merged-module counterpart of the 6.4.DT_signature_* scripts, consolidated into
# one file (no Slurm array). Same scoring / binarisation / DEG / aggregate /
# pseudobulk-DESeq2 workflow, but sourced from the merged CB+DT groupdeg
# (6.3.archetype_module_Jaccard.r's G1/G2/G3 signatures, derived from the
# Jaccard-clustered CB+DT archetype modules) instead of a single-arm groupdeg.
#
# Section 1 mirrors 6.4.DT_signature_analysis.r (per-sample scoring,
#   binarisation, Module_group label, DEG + pathway enrichment, spatial plots),
#   looped over all 8 samples in-process instead of via SLURM_ARRAY_TASK_ID.
#   Skips a sample already processed (per-sample meta.Rds present).
# Section 2 mirrors 6.4.DT_signature_analysis_aggregate.R (cross-sample
#   boxplots, positivity, Module_group composition, per-cell Wilcoxon,
#   spatial grids, recurrent pathway enrichment).
# Section 3 mirrors 6.4.DT_signature_DEseq2.R (pseudobulk DESeq2 on
#   Module_group vs "Neg", one contrast per non-Neg level).
#
#   Rscript 6.4.merged_module_analysis.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(patchwork)
  library(ggpubr)
  library(RColorBrewer)
  library(data.table)
  library(DESeq2, lib.loc = "~/R_Library/4.5")
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
})

# ── Config ────────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir <- "~/VisHD"
source(file.path(base_dir, "functions.R"))        # binarise_expression(), run_gsea_panel(), pathwayenrich_plot()

groupdeg <- readRDS(file.path(base_dir,
  "6.3.archetype_module_Jaccard",
  "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))
module_names <- paste0("module_", names(groupdeg))
banksy_col   <- "BANKSY_0.2_snn_res.1"            # only BANKSY clustering in tumour_srt.qs2

# gene sets for Module_group DEG pathway enrichment (run_gsea_panel(), functions.R)
Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

outdir <- file.path(base_dir, "6.4.merged_module_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pos_pal <- c(neg = "grey85", pos = "indianred")

# Fixed Module_group levels/palette (Neg + every signature combination) so the
# colours are consistent across samples in the aggregate spatial grid. Palette
# is generated generically (not hardcoded per group name) since the merged
# groupdeg group count/names are not fixed.
group_combos  <- unlist(lapply(seq_along(module_names), function(k)
  combn(sub("^module_", "", module_names), k, FUN = function(x) paste(x, collapse = "/"))))
group_levels  <- c("Neg", group_combos)
single_combos <- group_combos[!grepl("/", group_combos)]        # one signature positive
multi_combos  <- group_combos[grepl("/", group_combos)]         # >= 2 signatures positive
group_pal     <- c(Neg = "lightgrey",
                    setNames(brewer.pal(max(3, length(single_combos)), "Set1")[seq_along(single_combos)],
                             single_combos),
                    setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(multi_combos)),
                             multi_combos))[group_levels]

## ===============================================================================
## Section 1 — per-sample scoring, binarisation, DEG + enrichment, spatial plots
## (mirrors 6.4.DT_signature_analysis.r, looped over all 8 samples)
## ===============================================================================

for (s in samples) {

  sample_dir <- file.path(outdir, s)
  meta_f     <- file.path(sample_dir, "meta.Rds")

  bin_dir_s <- file.path(sample_dir, "binarisation")
  deg_dir_s <- file.path(sample_dir, "module_group_deg")
  png_dir_s <- file.path(sample_dir, "png")
  dir.create(bin_dir_s, showWarnings = FALSE, recursive = TRUE)
  dir.create(deg_dir_s, showWarnings = FALSE, recursive = TRUE)
  dir.create(png_dir_s, showWarnings = FALSE, recursive = TRUE)

  srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(srt_f)) { message("Skipping ", s, " — missing tumour_srt.qs2"); next }

  message("Loading ", s, " ...")
  srt <- qs_read(srt_f)

  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  if (file.exists(meta_f)) {
    message(s, " already processed (", meta_f, " found) — reusing saved metadata")
    saved_meta <- readRDS(meta_f)
    srt  <- AddMetaData(srt, saved_meta)
    keep <- module_names %in% colnames(saved_meta)
    if (!any(keep)) stop("no scored module columns found in ", meta_f)
  } else {
    # restrict each signature to genes present in this object (need >= 2 for scoring)
    feats <- lapply(groupdeg, function(g) intersect(g, rownames(srt)))
    keep  <- lengths(feats) >= 2
    if (!any(keep)) { message("no scorable signatures for ", s, " — skipping"); rm(srt); gc(); next }

    srt <- AddModuleScore(srt, features = feats[keep], name = "GMmod_", assay = score_assay)
    new_cols <- paste0("GMmod_", seq_len(sum(keep)))                  # appended in feats[keep] order
    names(srt@meta.data)[match(new_cols, names(srt@meta.data))] <- module_names[keep]

    # binarise each per-cell module score (GMM background threshold) → pos/neg factor
    for (m in module_names[keep]) {
      sc  <- setNames(srt@meta.data[[m]], colnames(srt))
      message("  binarising ", s, " ", m)
      bin <- binarise_expression(sc, verbose = TRUE,
               plot_out = file.path(bin_dir_s, sprintf("%s_binarisation.png", m)))
      srt@meta.data[[paste0(m, "_pos")]] <-
        factor(ifelse(bin == 1L, "pos", "neg"), levels = c("neg", "pos"))
    }

    # combine the per-signature pos calls into one label, e.g. "G1", "G1/G2", "Neg"
    kept_mods <- module_names[keep]
    labs      <- sub("^module_", "", kept_mods)
    pos_mat   <- do.call(cbind, lapply(kept_mods,
                    function(m) srt@meta.data[[paste0(m, "_pos")]] == "pos"))
    srt@meta.data$Module_group <- factor(apply(pos_mat, 1, function(r) {
      hit <- labs[which(r)]                           # drops NA + neg, keeps module order
      if (!length(hit)) "Neg" else paste(hit, collapse = "/")
    }), levels = group_levels)
  }

  # spatial map of the combined Module_group label
  sp_group_s <- ImageDimPlot(srt, group.by = "Module_group", cols = group_pal,
                     border.color = "#00000000", size = 0.4) +
    ggtitle(s) +
    theme(plot.title = element_text(size = 8), legend.text = element_text(size = 6))
  ggsave(file.path(png_dir_s, "sp_group.png"), sp_group_s,
         width = 6, height = 5, dpi = 150)

  # spatial map of positive cells, one ImageDimPlot per signature
  sp_pos_s <- list()
  for (m in paste0(module_names[keep], "_pos")) {
    npos <- sum(srt@meta.data[[m]] == "pos", na.rm = TRUE)
    mod  <- sub("_pos$", "", m)
    p    <- ImageDimPlot(srt, group.by = m, cols = pos_pal,
                   border.color = "#00000000", size = 0.4) +
      ggtitle(sprintf("%s  %s+ (%d)", s, mod, npos)) +
      theme(plot.title = element_text(size = 8), legend.position = "none")
    sp_pos_s[[mod]] <- p
    ggsave(file.path(png_dir_s, sprintf("sp_pos_%s.png", mod)), p,
           width = 6, height = 5, dpi = 150)
  }

  # ── DEG by Module_group + pathway enrichment on those DEGs ─────────────────
  srt$Module_group <- droplevels(srt$Module_group)
  deg_mg_f <- file.path(deg_dir_s, "deg_modulegroup.Rds")
  if (file.exists(deg_mg_f)) {
    message("  ", s, " Module_group DEG already detected (", deg_mg_f, ") — skipping")
  } else if (dplyr::n_distinct(srt$Module_group) >= 2) {
    Idents(srt) <- "Module_group"
    DEG_mg <- tryCatch(
      FindAllMarkers(srt, assay = score_assay, test.use = "MAST",
                     only.pos = TRUE, verbose = FALSE),
      error = function(e) { message("  FindAllMarkers failed for ", s, ": ",
                                    conditionMessage(e)); NULL })
    if (!is.null(DEG_mg) && nrow(DEG_mg) > 0) {
      saveRDS(DEG_mg, deg_mg_f)
      write.csv(DEG_mg, file.path(deg_dir_s, "deg_modulegroup.csv"),
                row.names = FALSE)
      run_gsea_panel(DEG_mg, gene_sets,
                     file.path(deg_dir_s, "enrich_modulegroup.Rds"))
    }
  }

  # ── DEG per-signature (pos vs neg, one direction: genes up in pos only) ─────
  for (m in module_names[keep]) {
    lab <- sub("^module_", "", m)
    col <- paste0(m, "_pos")
    deg_pos_f <- file.path(deg_dir_s, sprintf("deg_%s_pos.Rds", lab))
    if (file.exists(deg_pos_f)) {
      message("  ", s, " ", lab, " pos/neg DEG already detected (", deg_pos_f, ") — skipping")
    } else if (dplyr::n_distinct(srt@meta.data[[col]]) >= 2) {
      Idents(srt) <- col
      DEG_pos <- tryCatch(
        FindMarkers(srt, ident.1 = "pos", ident.2 = "neg", assay = score_assay,
                    test.use = "MAST", only.pos = TRUE, verbose = FALSE),
        error = function(e) { message("  FindMarkers failed for ", s, " ", lab, ": ",
                                      conditionMessage(e)); NULL })
      if (!is.null(DEG_pos) && nrow(DEG_pos) > 0) {
        DEG_pos <- tibble::rownames_to_column(DEG_pos, "gene")
        DEG_pos$cluster <- lab
        saveRDS(DEG_pos, deg_pos_f)
        write.csv(DEG_pos, file.path(deg_dir_s, sprintf("deg_%s_pos.csv", lab)),
                  row.names = FALSE)
        run_gsea_panel(DEG_pos, gene_sets,
                       file.path(deg_dir_s, sprintf("enrich_%s_pos.Rds", lab)))
      }
    }
  }

  # ── Save per-sample outputs for the aggregate step ────────────────────────────
  saveRDS(srt@meta.data, meta_f)
  qs_save(sp_group_s, file.path(sample_dir, "sp_group.qs2"))
  qs_save(sp_pos_s, file.path(sample_dir, "sp_pos.qs2"))

  rm(srt); gc()
  message("Done processing ", s)
}

## ===============================================================================
## Section 2 — cross-sample aggregate (mirrors 6.4.DT_signature_analysis_aggregate.R)
## ===============================================================================

box_dir  <- file.path(outdir, "boxplots")
sts_dir  <- file.path(outdir, "stats")
spt_dir  <- file.path(outdir, "spatial")
path_dir <- file.path(outdir, "pathway_enrichment")
corr_dir <- file.path(outdir, "correlation")
for (d in c(box_dir, sts_dir, spt_dir, path_dir, corr_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ── Load per-sample outputs ────────────────────────────────────────────────────
meta_list <- list()
sp_pos    <- setNames(vector("list", length(module_names)), module_names)  # ImageDimPlots
sp_group  <- list()                                                        # Module_group maps
for (s in samples) {
  meta_f <- file.path(outdir, s, "meta.Rds")
  if (!file.exists(meta_f)) { message("Skipping ", s, " — no meta.Rds (not yet processed)"); next }
  meta_list[[s]] <- readRDS(meta_f)
  sp_group_f <- file.path(outdir, s, "sp_group.qs2")
  if (file.exists(sp_group_f)) sp_group[[s]] <- qs_read(sp_group_f)
  sp_pos_f <- file.path(outdir, s, "sp_pos.qs2")
  if (file.exists(sp_pos_f)) {
    sp_pos_s <- qs_read(sp_pos_f)
    for (m in names(sp_pos_s)) sp_pos[[m]][[s]] <- sp_pos_s[[m]]
  }
}
if (!length(meta_list)) stop("No per-sample outputs found — Section 1 must process at least one sample first")

# ── Save per-cell binarisation results ────────────────────────────────────────
# idcol="slide" stamps each cell with its sample (per-sample objects lack one);
# fill=TRUE tolerates samples whose scorable signatures (hence columns) differ.
metas <- data.table::rbindlist(meta_list, idcol = "slide", fill = TRUE)
write.csv(metas, file.path(outdir, "metas.csv"), row.names = FALSE)
saveRDS(metas, file.path(outdir, "metas.Rds"))

# ── Long-format per-cell tables ───────────────────────────────────────────────
metas_df <- as.data.frame(metas)
# collapse category to broad CB vs DT (drop CB 0/CB 1 sub-levels) for all category comparisons
if ("category" %in% names(metas_df))
  metas_df$category <- ifelse(grepl("^CB", metas_df$category), "CB", metas_df$category)
present_modules <- intersect(module_names, names(metas_df))
present_pos     <- intersect(paste0(present_modules, "_pos"), names(metas_df))
present_modules <- sub("_pos$", "", present_pos)                  # keep score/pos paired
group_cols      <- intersect(c("category", "subclone", banksy_col, "Module_group"), names(metas_df))

score_long <- metas_df %>%
  select(slide, all_of(group_cols), all_of(present_modules)) %>%
  pivot_longer(all_of(present_modules), names_to = "module", values_to = "score")

pos_long <- metas_df %>%
  select(slide, all_of(group_cols), all_of(present_pos)) %>%
  pivot_longer(all_of(present_pos), names_to = "module", values_to = "pos") %>%
  mutate(module = sub("_pos$", "", module), pos = as.character(pos))

# ── Module score vs AR / FOLH1 expression correlation ────────────────────────
# Per-cell SpaNorm expression for AR/FOLH1 (5.4.tumour_expression_proportion.R
# output) joined to this script's per-cell module scores by (slide, cell).
# Scatter + Spearman correlation across all cells/samples pooled, one panel
# per module, plus a module-score boxplot split by AR/FOLH1-negative
# (expression == 0, i.e. zero counts) vs -positive cells.
expr_f <- file.path(base_dir, "5.4.expression_proportion", "spanorm_expression_per_cell.Rds")
if (!file.exists(expr_f)) {
  message("Skipping module vs AR/FOLH1 correlation — missing ", expr_f)
} else {
  ar_folh1 <- readRDS(expr_f) %>%
    filter(gene %in% c("AR", "FOLH1")) %>%
    transmute(slide = as.character(sample), cell = as.character(cell), gene, expression) %>%
    pivot_wider(names_from = gene, values_from = expression, names_prefix = "expr_")

  mod_wide <- metas_df %>%
    mutate(cell = as.character(cell)) %>%
    select(slide, cell, all_of(present_modules))

  corr_long <- mod_wide %>%
    inner_join(ar_folh1, by = c("slide", "cell")) %>%
    pivot_longer(all_of(present_modules), names_to = "module", values_to = "score")
  cat("Module x AR/FOLH1 correlation table:", nrow(corr_long), "rows (",
      dplyr::n_distinct(corr_long$cell), "cells matched)\n")

  # Scatter + Spearman correlation, one panel per module, pooled across samples
  make_module_scatter <- function(gene_col, gene_name, file) {
    d <- corr_long %>% filter(!is.na(score), !is.na(.data[[gene_col]]))
    if (nrow(d) == 0) return(invisible())
    p <- ggplot(d, aes(x = .data[[gene_col]], y = score)) +
      geom_point(size = 0.3, alpha = 0.15) +
      geom_smooth(method = "lm", se = FALSE, colour = "firebrick", linewidth = 0.5) +
      stat_cor(method = "spearman", size = 2.5, colour = "black") +
      facet_wrap(~ module, scales = "free") +
      labs(x = paste(gene_name, "expression (SpaNorm)"), y = "module score",
           title = sprintf("Module score vs %s expression (per cell, all samples)", gene_name)) +
      theme_bw(base_size = 9)
    ggsave(file, p, width = length(present_modules) * 3 + 1, height = 4,
           dpi = 200, limitsize = FALSE)
  }
  make_module_scatter("expr_AR",    "AR",    file.path(corr_dir, "scatter_module_vs_AR.png"))
  make_module_scatter("expr_FOLH1", "FOLH1", file.path(corr_dir, "scatter_module_vs_FOLH1.png"))

  # Boxplot: module score in gene-negative (expression == 0) vs -positive cells
  make_module_negpos_box <- function(gene_col, gene_name, file) {
    d <- corr_long %>%
      filter(!is.na(score), !is.na(.data[[gene_col]])) %>%
      mutate(grp = factor(ifelse(.data[[gene_col]] == 0, "neg", "pos"), levels = c("neg", "pos")))
    if (nrow(d) == 0 || dplyr::n_distinct(d$grp) < 2) return(invisible())
    p <- ggplot(d, aes(grp, score, fill = grp)) +
      geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.15, lwd = 0.3) +
      facet_wrap(~ module, scales = "free_y") +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.6) +
      scale_fill_manual(values = pos_pal, guide = "none") +
      labs(x = sprintf("%s status (counts = 0 -> neg)", gene_name), y = "module score",
           title = sprintf("Module score by %s-negative vs %s-positive cells", gene_name, gene_name)) +
      theme_bw(base_size = 9)
    ggsave(file, p, width = length(present_modules) * 2.5 + 1, height = 4.5,
           dpi = 200, limitsize = FALSE)
  }
  make_module_negpos_box("expr_AR",    "AR",    file.path(corr_dir, "boxplot_module_by_AR_status.png"))
  make_module_negpos_box("expr_FOLH1", "FOLH1", file.path(corr_dir, "boxplot_module_by_FOLH1_status.png"))
}

# ── Per-cell boxplots: module score by group, facet module (rows) × sample (cols)
# Each point is one cell; a per-panel Wilcoxon/Kruskal p-value is annotated.
make_cell_box <- function(group_col, drop_vals, title, file) {
  d <- score_long %>%
    filter(!is.na(.data[[group_col]]), !.data[[group_col]] %in% drop_vals, !is.na(score)) %>%
    mutate(grp = factor(.data[[group_col]]))
  if (nrow(d) == 0 || dplyr::n_distinct(d$grp) < 2) return(invisible())
  p <- ggplot(d, aes(grp, score, fill = grp)) +
    geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.15, lwd = 0.3) +
    facet_grid(module ~ slide, scales = "free", space = "free_x") +
    stat_compare_means(label = "p.format", size = 2,
                       method = if (nlevels(d$grp) > 2) "kruskal.test" else "wilcox.test") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "module score", title = title) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  ggsave(file, p, width = length(samples) * 2.2 + 1,
         height = length(present_modules) * 2.5 + 1, dpi = 200, limitsize = FALSE)
}
make_cell_box("category", character(0),
  "Per-cell module score by category", file.path(box_dir, "boxplot_category.png"))
make_cell_box("subclone", c("Normal", "Removed"),
  "Per-cell module score by subclone", file.path(box_dir, "boxplot_subclone.png"))
make_cell_box("Module_group", character(0),
  "Per-cell module score by Module_group", file.path(box_dir, "boxplot_module_group.png"))

# ── Positivity bar plots: proportion positive per group, facet module × sample ─
make_pos_bar <- function(group_col, drop_vals, palette, title, file) {
  d <- pos_long %>%
    filter(!is.na(.data[[group_col]]), !.data[[group_col]] %in% drop_vals, !is.na(pos)) %>%
    mutate(grp = factor(.data[[group_col]])) %>%
    group_by(slide, module, grp) %>%
    summarise(prop_pos = mean(pos == "pos"), n = dplyr::n(), .groups = "drop")
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(grp, prop_pos, fill = grp)) +
      geom_col() +
      facet_grid(module ~ slide, scales = "free_x", space = "free_x") +
      scale_fill_brewer(palette = palette, guide = "none") +
      labs(x = NULL, y = "proportion positive", title = title) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
    ggsave(file, p, width = length(samples) * 2.6 + 1,
           height = length(present_modules) * 2.2 + 1, dpi = 200, limitsize = FALSE)
  }
  d %>% mutate(grouping = group_col)
}
bar_cat <- make_pos_bar("category", character(0), "Set2",
  "Signature positivity by category", file.path(box_dir, "positivity_barplot_category.png"))
bar_sub <- make_pos_bar("subclone", c("Normal", "Removed"), "Set1",
  "Signature positivity by subclone", file.path(box_dir, "positivity_barplot_subclone.png"))
bar_bk  <- make_pos_bar(banksy_col, character(0), "Set3",
  "Signature positivity by BANKSY_0.2 cluster", file.path(box_dir, "positivity_barplot_banksy.png"))

write.csv(bind_rows(bar_cat, bar_sub, bar_bk),
          file.path(sts_dir, "positivity_proportions.csv"), row.names = FALSE)

# ── Positivity boxplot by category: one point per sample, Wilcoxon significance
make_pos_box <- function(d, title, file) {
  d <- d %>% mutate(grp = droplevels(factor(grp)))
  if (nrow(d) == 0 || dplyr::n_distinct(d$grp) < 2) return(invisible())
  comparisons <- combn(levels(d$grp), 2, simplify = FALSE)
  p <- ggplot(d, aes(grp, prop_pos)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92") +
    geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.8) +
    facet_wrap(~ module, scales = "free_y") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", size = 2.5) +
    labs(x = "category", y = "proportion positive (per sample)", title = title) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file, p, width = length(present_modules) * 2.5 + 1, height = 4.5,
         dpi = 200, limitsize = FALSE)
}
make_pos_box(bar_cat, "Signature positivity by category (points = samples)",
  file.path(box_dir, "boxplot_positivity_by_category.png"))
make_pos_box(bar_sub, "Signature positivity by subclone (points = samples)",
  file.path(box_dir, "boxplot_positivity_by_subclone.png"))
make_pos_box(bar_bk, "Signature positivity by BANKSY_0.2 cluster (points = samples)",
  file.path(box_dir, "boxplot_positivity_by_banksy.png"))

# ── Module_group composition: proportion of each label per group, facet sample ─
make_group_comp <- function(group_col, drop_vals, title, file) {
  d <- metas_df %>%
    filter(!is.na(.data[[group_col]]), !.data[[group_col]] %in% drop_vals,
           !is.na(Module_group)) %>%
    mutate(grp = droplevels(factor(.data[[group_col]])),    # drop x levels absent after filtering
           Module_group = factor(Module_group, levels = group_levels)) %>%
    group_by(slide, grp, Module_group, .drop = T) %>%
    summarise(n = dplyr::n(), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  if (nrow(d) == 0) return(invisible())
  p <- ggplot(d, aes(grp, prop, fill = Module_group)) +
    geom_col() +
    facet_grid(.~slide, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = group_pal, drop = TRUE) +
    labs(x = NULL, y = "proportion of cells", fill = "Module_group", title = title) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  ggsave(file, p, width = length(samples) * 2.6 + 1, height = 4.5, dpi = 200, limitsize = FALSE)
  d %>% mutate(grouping = group_col)
}
comp_cat <- make_group_comp("category", character(0),
  "Module_group composition by category", file.path(box_dir, "module_group_composition_category.png"))
comp_sub <- make_group_comp("subclone", c("Normal", "Removed"),
  "Module_group composition by subclone", file.path(box_dir, "module_group_composition_subclone.png"))
comp_bk  <- make_group_comp(banksy_col, character(0),
  "Module_group composition by BANKSY_0.2 cluster", file.path(box_dir, "module_group_composition_banksy.png"))
write.csv(bind_rows(comp_cat, comp_sub, comp_bk),
          file.path(sts_dir, "module_group_composition.csv"), row.names = FALSE)

# ── Per-cell Wilcoxon: pairwise module-score test per (sample, module, grouping)
cell_wilcox <- function(group_col, drop_vals) {
  d <- score_long %>%
    filter(!is.na(.data[[group_col]]), !.data[[group_col]] %in% drop_vals, !is.na(score))
  out <- list()
  for (s in unique(d$slide)) for (m in unique(d$module)) {
    sub  <- d %>% filter(slide == s, module == m)
    grps <- sort(unique(as.character(sub[[group_col]])))
    if (length(grps) < 2) next
    for (pr in combn(grps, 2, simplify = FALSE)) {
      a <- sub$score[sub[[group_col]] == pr[1]]
      b <- sub$score[sub[[group_col]] == pr[2]]
      if (length(a) < 3 || length(b) < 3) next
      p <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_)
      out[[length(out) + 1]] <- tibble(slide = s, module = m, group1 = pr[1], group2 = pr[2],
                                       n1 = length(a), n2 = length(b), p = p)
    }
  }
  res <- bind_rows(out)
  if (nrow(res)) res$p_adj <- p.adjust(res$p, method = "BH")
  res
}
write.csv(bind_rows(mutate(cell_wilcox("category", character(0)), grouping = "category"),
                    mutate(cell_wilcox("subclone", c("Normal", "Removed")), grouping = "subclone")),
          file.path(sts_dir, "wilcoxon_percell_modulescore.csv"), row.names = FALSE)

# ── Spatial grid: signatures (rows) × samples (columns), positive cells red ───
build_grid <- function(plotmat) {
  cells <- list()
  for (m in module_names) for (s in samples) {
    p <- plotmat[[m]][[s]]
    cells[[length(cells) + 1]] <- if (is.null(p)) patchwork::plot_spacer() else p
  }
  wrap_plots(cells, nrow = length(module_names), ncol = length(samples), byrow = TRUE)
}
if (any(lengths(sp_pos) > 0))
  ggsave(file.path(spt_dir, "positive_cells_imagedimplot_all.png"), build_grid(sp_pos),
         width = length(samples) * 4, height = length(module_names) * 4,
         dpi = 150, limitsize = FALSE)

# ── Spatial grid: Module_group label, one panel per sample ────────────────────
if (length(sp_group) > 0)
  ggsave(file.path(spt_dir, "module_group_imagedimplot_all.png"),
         wrap_plots(sp_group[samples[samples %in% names(sp_group)]], nrow = 2) &
           scale_fill_manual(values = group_pal, drop = FALSE),
         width = length(sp_group) * 3, height = 10, dpi = 150, limitsize = FALSE)

# ── Aggregate Module_group pathway enrichment across samples ──────────────────
# Recurrent = significant (p.adjust < 0.05) in >= min_recur samples for a given
# (Module_group, geneset, pathway); top_n most-recurrent pathways per Module_group
# are plotted as a bar chart.
enrich_files <- list.files(outdir, pattern = "^enrich_modulegroup\\.Rds$",
                           recursive = TRUE, full.names = TRUE)
if (length(enrich_files) > 0) {
  enrich_long <- bind_rows(lapply(enrich_files, function(f) {
    s  <- basename(dirname(dirname(f)))          # .../<sample>/module_group_deg/enrich_modulegroup.Rds
    el <- readRDS(f)
    bind_rows(lapply(names(el), function(gs) bind_rows(lapply(names(el[[gs]]), function(m) {
      res <- el[[gs]][[m]]
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      as.data.frame(res@result) %>% mutate(slide = s, geneset = gs, Module_group = m)
    }))))
  }))
  write.csv(enrich_long, file.path(path_dir, "all_pathway_enrichment_long.csv"), row.names = FALSE)

  min_recur <- 2   # a pathway must be significant in >=2 samples to count as "recurrent"
  top_n     <- 10
  recur <- enrich_long %>%
    filter(p.adjust < 0.05) %>%
    group_by(Module_group, geneset, ID) %>%
    summarise(n_recur = dplyr::n_distinct(slide), mean_NES = mean(NES), .groups = "drop") %>%
    filter(n_recur >= min_recur)
  write.csv(recur, file.path(path_dir, "recurrent_pathways.csv"), row.names = FALSE)

  for (m in unique(recur$Module_group)) {
    d <- recur %>% filter(Module_group == m) %>%
      arrange(desc(n_recur), desc(abs(mean_NES))) %>% slice_head(n = top_n) %>%
      mutate(ID = factor(ID, levels = rev(ID)))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(ID, n_recur, fill = geneset)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = sprintf("# samples significant (of %d)", length(samples)),
           fill = "gene set", title = sprintf("Top recurrent enriched pathways — %s", m)) +
      theme_bw(base_size = 9)
    ggsave(file.path(path_dir, sprintf("top_pathways_%s.png", gsub("/", "-", m))), p,
           width = 8, height = max(3, nrow(d) * 0.35 + 1), dpi = 150, limitsize = FALSE)
  }
}

# ── Aggregate per-signature (pos vs neg) pathway enrichment across samples ────
# Recurrent = significant (p.adjust < 0.05) in >= min_recur samples for a given
# (signature, geneset, pathway); top_n most-recurrent pathways per signature
# are plotted as a bar chart.
enrich_pos_files <- list.files(outdir, pattern = "^enrich_.*_pos\\.Rds$",
                               recursive = TRUE, full.names = TRUE)
if (length(enrich_pos_files) > 0) {
  enrich_pos_long <- bind_rows(lapply(enrich_pos_files, function(f) {
    s   <- basename(dirname(dirname(f)))          # .../<sample>/module_group_deg/enrich_<sig>_pos.Rds
    lab <- sub("^enrich_(.*)_pos\\.Rds$", "\\1", basename(f))
    el  <- readRDS(f)
    bind_rows(lapply(names(el), function(gs) bind_rows(lapply(names(el[[gs]]), function(m) {
      res <- el[[gs]][[m]]
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      as.data.frame(res@result) %>% mutate(slide = s, geneset = gs, signature = lab)
    }))))
  }))
  write.csv(enrich_pos_long, file.path(path_dir, "all_pathway_enrichment_signature_pos_long.csv"),
            row.names = FALSE)

  min_recur_pos <- 1   # a pathway must be significant in >=1 sample to count as "recurrent"
  top_n_pos     <- 15
  recur_pos <- enrich_pos_long %>%
    filter(p.adjust < 0.05) %>%
    filter(NES > 0) %>%
    group_by(signature, geneset, ID) %>%
    summarise(n_recur = dplyr::n_distinct(slide), mean_NES = mean(NES), .groups = "drop") %>%
    filter(n_recur >= min_recur_pos)
  write.csv(recur_pos, file.path(path_dir, "recurrent_pathways_signature_pos.csv"), row.names = FALSE)

  for (lab in unique(recur_pos$signature)) {
    d <- recur_pos %>% filter(signature == lab) %>%
      arrange(desc(n_recur), desc(abs(mean_NES))) %>% slice_head(n = top_n_pos) %>%
      mutate(ID = str_wrap(gsub("_", " ", ID), width = 20),
             ID = factor(ID, levels = rev(ID)))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(ID, mean_NES, fill = geneset)) +
      geom_col() +
      geom_text(aes(label = n_recur), hjust = -0.2, size = 3) +
      coord_flip() +
      labs(x = NULL, y = "Average NES score",
           fill = "gene set", title = sprintf("Top recurrent enriched pathways — %s+ (pos vs neg)", lab)) +
      theme_bw(base_size = 9)
    ggsave(file.path(path_dir, sprintf("top_pathways_%s_pos.png", gsub("/", "-", lab))), p,
           width = 8, height = max(3, nrow(d) * 0.35 + 1), dpi = 150, limitsize = FALSE)
  }
}

message("\nDone. 6.4.merged_module_analysis aggregate written to ", outdir)

## ===============================================================================
## Section 3 — pseudobulk DESeq2 on Module_group vs "Neg" (mirrors 6.4.DT_signature_DEseq2.R)
## ===============================================================================

deseq_outdir <- file.path(outdir, "DESeq2")
png_dir      <- file.path(deseq_outdir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

min_cells       <- 10   # minimum cells to keep a pseudobulk sample
min_samples_grp <- 2    # a Module_group level needs pseudobulks in >=2 samples to be tested

# Gene sets for per-contrast GSEA (same collections as 5.3.integrate_DESeq2_deg.R)
archetype_module <- readRDS(file.path(base_dir, "public_signature/clean_module.Rds"))
meta_programs <- readxl::read_excel(
  file.path(base_dir, "public_signature/meta_programs_2025-01-29.xlsx"),
  sheet = "Malignant", col_names = TRUE
) |>
  purrr::map(~ as.character(na.omit(.x)))

list_to_term2gene <- function(lst) {
  data.frame(
    term = rep(names(lst), lengths(lst)),
    gene = unlist(lst, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}
Archetype <- list_to_term2gene(archetype_module)
MetaProgM <- list_to_term2gene(meta_programs)
gene_sets_deseq <- list(Hallmark = Hall, C6 = C6, C5 = C5,
                        Archetype = Archetype, MetaProgMalignant = MetaProgM)

# ── 1. Load all tumour objects + join per-cell Module_group ──────────────────
srt_list <- list()
for (s in samples) {
  srt_f  <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  meta_f <- file.path(outdir, s, "meta.Rds")
  if (!file.exists(srt_f))  { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
  if (!file.exists(meta_f)) { message("Skipping ", s, " — missing 6.4.merged_module_analysis meta.Rds"); next }
  message("Loading ", s, " ...")
  obj  <- qs_read(srt_f)
  meta <- readRDS(meta_f)
  obj$Module_group <- as.character(meta$Module_group[match(colnames(obj), rownames(meta))])
  obj$sample_id    <- s
  srt_list[[s]] <- obj
}
cat(length(srt_list), "samples loaded\n")

# ── 2. Merge and save integrated object ───────────────────────────────────────
merged_f <- file.path(deseq_outdir, "merged_tumour_srt.qs2")
if (file.exists(merged_f)) {
  cat("Loading existing merged object ...\n")
  srt_merged <- qs_read(merged_f)
} else {
  cat("Merging ...\n")
  srt_merged <- merge(
    srt_list[[1]],
    y            = srt_list[-1],
    add.cell.ids = names(srt_list),
    merge.data   = FALSE
  )
  srt_merged <- JoinLayers(srt_merged, assay = "Spatial")
  qs_save(srt_merged, merged_f)
  cat("Integrated object saved.\n")
}
cat("Merged:", ncol(srt_merged), "cells,", nrow(srt_merged), "genes\n")

# ── 3. Pseudobulk aggregation (sample x Module_group) ─────────────────────────
meta        <- srt_merged@meta.data %>% filter(!is.na(Module_group))
valid_cells <- rownames(meta)

counts_raw <- GetAssayData(srt_merged, assay = "Spatial", layer = "counts")[, valid_cells]

group_ids     <- setNames(paste0(meta$sample_id, "__", meta$Module_group), valid_cells)
unique_groups <- unique(group_ids)
cat(length(unique_groups), "pseudobulk groups before filtering\n")

# Sum counts per group (sparse-friendly)
pb_counts <- vapply(unique_groups, function(g) {
  cells_g <- names(group_ids)[group_ids == g]
  Matrix::rowSums(counts_raw[, cells_g, drop = FALSE])
}, numeric(nrow(counts_raw)))

pb_meta <- data.frame(
  group        = unique_groups,
  sample_id    = sub("__.*", "", unique_groups),
  Module_group = sub(".*__", "", unique_groups),
  n_cells      = as.integer(table(group_ids)[unique_groups]),
  row.names    = unique_groups,
  stringsAsFactors = FALSE
)
pb_meta   <- pb_meta %>% filter(n_cells >= min_cells)
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
cat(nrow(pb_meta), "pseudobulks retained (>=", min_cells, "cells)\n")

# Keep only Module_group levels seen in >= min_samples_grp distinct samples —
# a level present in only one sample would alias with that sample's blocking
# term and make the ~sample_id + Module_group design matrix rank-deficient.
grp_n_samples <- pb_meta %>%
  dplyr::distinct(sample_id, Module_group) %>%
  dplyr::count(Module_group, name = "n_samples")
keep_groups <- grp_n_samples %>% dplyr::filter(n_samples >= min_samples_grp) %>% dplyr::pull(Module_group)
if (!"Neg" %in% keep_groups) stop("'Neg' Module_group not present in enough samples — cannot run contrasts.")

pb_meta   <- pb_meta %>% filter(Module_group %in% keep_groups)
pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
cat(length(keep_groups), "Module_group levels retained:", paste(keep_groups, collapse = ", "), "\n")

if (nrow(pb_meta) < 4) stop("Too few pseudobulk samples for DESeq2 — check cell numbers.")

# Filter lowly expressed genes: ≥10 counts in ≥(n_pseudobulk / 5) samples
keep_genes <- rowSums(pb_counts >= 10) >= max(2, floor(nrow(pb_meta) / 5))
pb_counts  <- pb_counts[keep_genes, ]
cat(sum(keep_genes), "genes retained after count filtering\n")

# ── 4. DESeq2 (single fit; contrasts pulled per Module_group level below) ─────
pb_meta$Module_group <- factor(pb_meta$Module_group, levels = c("Neg", setdiff(keep_groups, "Neg")))
pb_meta$sample_id    <- factor(pb_meta$sample_id)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_counts),
  colData   = pb_meta,
  design    = ~sample_id + Module_group
)
dds <- DESeq(dds, parallel = FALSE)
cat("DESeq2 done\n")

saveRDS(dds, file.path(deseq_outdir, "deseq2_dds.Rds"))
write.csv(pb_meta, file.path(deseq_outdir, "pseudobulk_meta.csv"), row.names = FALSE)

# PCA of pseudobulk samples
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("Module_group", "sample_id"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = Module_group, shape = sample_id)) +
  geom_point(size = 3) +
  labs(title = "PCA — pseudobulk samples (VST)",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       colour = "Module_group", shape = "Sample") +
  theme_classic()
ggsave(file.path(png_dir, "0_PCA_pseudobulk.png"), p_pca, width = 9, height = 6, dpi = 200)

# ── 5. Module_group vs Neg contrasts ──────────────────────────────────────────
contrast_levels <- setdiff(levels(pb_meta$Module_group), "Neg")
summary_list <- list()

for (lvl in contrast_levels) {
  lvl_safe <- gsub("[^A-Za-z0-9]+", "_", lvl)
  cat("\nContrast:", lvl, "vs Neg\n")

  res    <- results(dds, contrast = c("Module_group", lvl, "Neg"), alpha = 0.05)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj)

  saveRDS(res_df, file.path(deseq_outdir, sprintf("deseq2_res_%s_vs_Neg.Rds", lvl_safe)))
  write.csv(res_df, file.path(deseq_outdir, sprintf("deseq2_res_%s_vs_Neg.csv", lvl_safe)),
            row.names = FALSE)

  n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5)
  cat("  significant DEGs (padj<0.05, |lfc|>0.5):", n_sig, "\n")
  summary_list[[lvl]] <- data.frame(Module_group = lvl, n_sig_DEGs = n_sig)

  # Volcano plot
  p_vol <- res_df %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      sig       = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.25,
      direction = case_when(
        sig & log2FoldChange > 0 ~ "Up",
        sig & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(sig, gene, NA_character_)
    ) %>%
    ggplot(aes(log2FoldChange, -log10(pvalue + 1e-300), colour = direction, label = label)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(size = 2.5, max.overlaps = 25, na.rm = TRUE) +
    scale_colour_manual(values = c(NS = "grey70", Up = "firebrick", Down = "royalblue")) +
    geom_vline(xintercept = c(-1.25, 1.25), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = sprintf("%s vs Neg — pseudobulk DESeq2", lvl),
         subtitle = paste0(n_sig, " DEGs  |  positive = ", lvl),
         x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() + theme(legend.position = "none")
  ggsave(file.path(png_dir, sprintf("volcano_%s_vs_Neg.png", lvl_safe)),
         p_vol, width = 7, height = 5, dpi = 200)

  # Heatmap of top DEGs (VST z-scores), restricted to this contrast's Neg + lvl samples
  sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% mutate(score = abs(log2FoldChange) / padj)
  top_up    <- sig_genes %>% filter(log2FoldChange > 0) %>% slice_max(score, n = 30)
  top_down  <- sig_genes %>% filter(log2FoldChange < 0) %>% slice_max(score, n = 30)
  top_genes <- bind_rows(top_up, top_down) %>% pull(gene)

  if (length(top_genes) >= 5) {
    hm_samples <- rownames(pb_meta)[pb_meta$Module_group %in% c("Neg", lvl)]
    mat_hm <- assay(vsd)[top_genes, hm_samples, drop = FALSE]
    mat_hm <- t(scale(t(mat_hm)))
    mat_hm[mat_hm >  3] <-  3
    mat_hm[mat_hm < -3] <- -3

    col_fun   <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
    grp_col   <- setNames(c("grey70", "firebrick"), c("Neg", lvl))
    samp_lvls <- levels(pb_meta$sample_id)
    samp_col  <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(samp_lvls)), samp_lvls)

    ha <- HeatmapAnnotation(
      Module_group = pb_meta[hm_samples, "Module_group"],
      Sample       = pb_meta[hm_samples, "sample_id"],
      col          = list(Module_group = grp_col, Sample = samp_col),
      annotation_name_side = "left"
    )

    png(file.path(png_dir, sprintf("heatmap_%s_vs_Neg.png", lvl_safe)),
        width = 1600, height = 2000, res = 200)
    draw(Heatmap(
      mat_hm,
      name              = "z-score",
      col               = col_fun,
      top_annotation    = ha,
      show_row_names    = TRUE,
      show_column_names = TRUE,
      row_names_gp      = gpar(fontsize = 7),
      column_names_gp   = gpar(fontsize = 7),
      column_title      = sprintf("Top %d DEGs — %s vs Neg", length(top_genes), lvl),
      clustering_method_columns = "ward.D2",
      column_split      = droplevels(pb_meta[hm_samples, "Module_group"])
    ))
    dev.off()
    cat("  heatmap saved\n")
  }

  # GSEA enrichment on this contrast's significant genes
  sig_df <- res_df %>% filter(!is.na(padj), padj < 0.05) %>% arrange(desc(log2FoldChange))
  if (nrow(sig_df) > 0) {
    gene_list   <- setNames(sig_df$log2FoldChange, sig_df$gene)
    enrich_list <- lapply(gene_sets_deseq, function(gs) {
      tryCatch(
        clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
        error   = function(e) NULL,
        warning = function(w) NULL
      )
    })
    saveRDS(enrich_list, file.path(deseq_outdir, sprintf("enrich_%s_vs_Neg.Rds", lvl_safe)))

    for (nm in names(enrich_list)) {
      res_e <- enrich_list[[nm]]
      if (is.null(res_e) || nrow(res_e@result) == 0) next
      sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
      if (sig_n == 0) next
      p_e <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e)
      ggsave(file.path(png_dir, sprintf("GSEA_%s_vs_Neg_%s.pdf", lvl_safe, nm)),
             p_e, width = 6, height = 10)
    }
    cat("  GSEA done\n")
  } else {
    cat("  no significant DEGs — skipping GSEA\n")
  }
}

summary_df <- bind_rows(summary_list) %>% arrange(desc(n_sig_DEGs))
write.csv(summary_df, file.path(deseq_outdir, "DEG_summary.csv"), row.names = FALSE)
print(summary_df)

cat("\nDone. 6.4.merged_module_analysis complete. Outputs in", outdir, "\n")
