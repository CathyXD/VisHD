# 6.2.3 — groupdeg signature analysis across tumour cells
# Score the final G1/G2/G3 signatures (from 6.2.2 groupdeg.rds) on every sample's
# tumour cells, binarise each per-cell module score, then compare per-cell module
# scores between `category` and tumour `subclone` (boxplots, facet module × sample),
# show signature positivity as bar plots grouped by category / subclone /
# BANKSY_0.2 cluster, and map the positive cells spatially. Saves test results,
# binarisation calls, and positivity proportions.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(patchwork)
  library(ggpubr)
  library(RColorBrewer)
})

# ── Config ──────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir <- "~/VisHD"
source(file.path(base_dir, "functions.R"))        # binarise_expression()

groupdeg <- readRDS(file.path(base_dir,
  "6.3.archetype_module",
  "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))
module_names <- paste0("module_", names(groupdeg))
banksy_col   <- "BANKSY_0.2_snn_res.1"            # only BANKSY clustering in tumour_srt.qs2

# gene sets for Module_group DEG pathway enrichment (run_gsea_panel(), functions.R)
Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

outdir   <- file.path(base_dir, "6.4.signature_analysis")
bin_dir  <- file.path(outdir, "binarisation")
box_dir  <- file.path(outdir, "boxplots")
sts_dir  <- file.path(outdir, "stats")
spt_dir  <- file.path(outdir, "spatial")
deg_dir  <- file.path(outdir, "module_group_deg")
path_dir <- file.path(outdir, "pathway_enrichment")
corr_dir <- file.path(outdir, "correlation")
for (d in c(bin_dir, box_dir, sts_dir, spt_dir, deg_dir, path_dir, corr_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

pos_pal <- c(neg = "grey85", pos = "indianred")

# ── Per-sample scoring + binarisation ─────────────────────────────────────────
# Fixed Module_group levels/palette (Neg + every signature combination) so the
# colours are consistent across samples in the spatial grid.
group_combos <- unlist(lapply(seq_along(module_names), function(k)
  combn(sub("^module_", "", module_names), k, FUN = function(x) paste(x, collapse = "/"))))
group_levels <- c("Neg", group_combos)
group_pal    <- c("Neg" = "lightgrey", "G1" = "red", "G2" = "gold", "G3" = "royalblue",
                  "G2/G1" = "orange", "G3/G1" = "purple", "G3/G2" = "green",
                  "G3/G2/G1" = "darkgrey")[group_levels]

sp_pos   <- setNames(vector("list", length(module_names)), module_names)  # ImageDimPlots
sp_group <- list()                                                        # Module_group maps
meta_list <- list()

# cached scoring/binarisation from a previous run — split by slide so each
# sample's cached columns can be merged back onto its freshly loaded srt
metas_f <- file.path(bin_dir, "metas.Rds")
cached_by_sample <- NULL
if (file.exists(metas_f)) {
  message("Found cached ", metas_f, " — reusing scoring/binarisation, skipping recompute")
  cached_metas     <- as.data.frame(readRDS(metas_f))
  metas_df <- cached_metas
  cached_by_sample <- split(cached_metas, cached_metas$slide)
  meta_list <- cached_by_sample
}

for (s in samples) {
  srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(srt_f)) { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
  message("Scoring groupdeg signatures for ", s, " ...")
  srt <- qs_read(srt_f)

  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  if (!is.null(cached_by_sample) && s %in% names(cached_by_sample)) {
    message("  using cached scoring/binarisation for ", s)
    meta_s <- cached_by_sample[[s]]
    srt  <- AddMetaData(srt, meta_s)
    keep <- module_names %in% intersect(module_names, names(meta_s))
    if (!any(keep)) { message("  no scorable signatures for ", s); rm(srt); gc(); next }
  } else {
    # restrict each signature to genes present in this object (need >= 2 for scoring)
    feats <- lapply(groupdeg, function(g) intersect(g, rownames(srt)))
    keep  <- lengths(feats) >= 2
    if (!any(keep)) { message("  no scorable signatures for ", s); rm(srt); gc(); next }

    srt <- AddModuleScore(srt, features = feats[keep], name = "GDmod_", assay = score_assay)
    new_cols <- paste0("GDmod_", seq_len(sum(keep)))                  # appended in feats[keep] order
    names(srt@meta.data)[match(new_cols, names(srt@meta.data))] <- module_names[keep]

    # binarise each per-cell module score (GMM background threshold) → pos/neg factor
    for (m in module_names[keep]) {
      sc  <- setNames(srt@meta.data[[m]], colnames(srt))
      message("  binarising ", s, " ", m)
      bin <- binarise_expression(sc, verbose = TRUE,
               plot_out = file.path(bin_dir, sprintf("%s_%s_binarisation.png", s, m)))
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
      # collect per-cell scores + pos calls + grouping metadata
  meta_list[[s]] <- tibble::rownames_to_column(srt@meta.data, "cell")
  }



  # spatial map of the combined Module_group label
  sp_group[[s]] <- ImageDimPlot(srt, group.by = "Module_group", cols = group_pal,
                     border.color = "#00000000", size = 0.4) +
    ggtitle(s) +
    theme(plot.title = element_text(size = 8), legend.text = element_text(size = 6))

  # spatial map of positive cells, one ImageDimPlot per signature
  for (m in paste0(module_names[keep], "_pos")) {
    npos <- sum(srt@meta.data[[m]] == "pos", na.rm = TRUE)
    sp_pos[[sub("_pos$", "", m)]][[s]] <-
      ImageDimPlot(srt, group.by = m, cols = pos_pal,
                   border.color = "#00000000", size = 0.4) +
      ggtitle(sprintf("%s  %s+ (%d)", s, sub("_pos$", "", m), npos)) +
      theme(plot.title = element_text(size = 8), legend.position = "none")
  }

  # ── DEG by Module_group + pathway enrichment on those DEGs ─────────────────
  srt$Module_group <- droplevels(srt$Module_group)
  deg_mg_f <- file.path(deg_dir, sprintf("%s_deg_modulegroup.Rds", s))
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
      write.csv(DEG_mg, file.path(deg_dir, sprintf("%s_deg_modulegroup.csv", s)),
                row.names = FALSE)
      run_gsea_panel(DEG_mg, gene_sets,
                     file.path(deg_dir, sprintf("%s_enrich_modulegroup.Rds", s)))
    }
  }

  # ── DEG per-signature (pos vs neg, one direction: genes up in pos only) ────
  for (m in module_names[keep]) {
    lab <- sub("^module_", "", m)
    col <- paste0(m, "_pos")
    deg_pos_f <- file.path(deg_dir, sprintf("%s_deg_%s_pos.Rds", s, lab))
    if (file.exists(deg_pos_f)) {
      message("  ", s, " ", lab, " pos DEG already detected (", deg_pos_f, ") — skipping")
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
        write.csv(DEG_pos, file.path(deg_dir, sprintf("%s_deg_%s_pos.csv", s, lab)),
                  row.names = FALSE)
        run_gsea_panel(DEG_pos, gene_sets,
                       file.path(deg_dir, sprintf("%s_enrich_%s_pos.Rds", s, lab)))
      }
    }
  }

  rm(srt); gc()
}


# ── Save per-cell binarisation results ────────────────────────────────────────
# idcol="slide" stamps each cell with its sample (per-sample objects lack one);
# fill=TRUE tolerates samples whose scorable signatures (hence columns) differ.
metas <- data.table::rbindlist(meta_list, idcol = "slide", fill = TRUE)
write.csv(metas, file.path(bin_dir, "metas.csv"), row.names = FALSE)
saveRDS(metas, file.path(bin_dir, "metas.Rds"))

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
      theme_bw(base_size = 15) +
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
    theme_bw(base_size = 15) +
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

# ── Module_group composition across all samples: facet_grid(category ~ slide) ─
comp_cat_slide <- metas_df %>%
  filter(!is.na(category), !is.na(Module_group)) %>%
  mutate(Module_group = factor(Module_group, levels = group_levels)) %>%
  group_by(slide, category, Module_group, .drop = TRUE) %>%
  summarise(n = dplyr::n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write.csv(comp_cat_slide, file.path(sts_dir, "module_group_composition_category_slide.csv"),
          row.names = FALSE)

if (nrow(comp_cat_slide) > 0) {
  p <- ggplot(comp_cat_slide, aes(x = "", y = prop, fill = Module_group)) +
    geom_col() +
    facet_grid(category ~ slide, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = group_pal, drop = TRUE) +
    labs(x = NULL, y = "proportion of cells", fill = "Module_group",
         title = "Module_group composition by category and slide") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(box_dir, "module_group_composition_category_slide.png"), p,
         width = length(samples) * 1.5 + 1, height = 4.5, dpi = 200, limitsize = FALSE)
}

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
  wrap_plots(cells, nrow = length(module_names), ncol = length(samples), byrow = TRUE, guides = "collect") 
}
if (any(lengths(sp_pos) > 0))
  ggsave(file.path(spt_dir, "positive_cells_imagedimplot_all.png"), build_grid(sp_pos),
         width = length(samples) * 4, height = length(module_names) * 4,
         dpi = 300, limitsize = FALSE)

# ── Spatial grid: Module_group label, one panel per sample ────────────────────
if (length(sp_group) > 0)
  ggsave(file.path(spt_dir, "module_group_imagedimplot_all.png"),
         wrap_plots(sp_group[samples[samples %in% names(sp_group)]], nrow = 2, guides = "collect") &
           scale_fill_manual(values = group_pal, drop = FALSE),
         width = length(sp_group) * 3, height = 10, dpi = 300, limitsize = FALSE)

# ── Aggregate Module_group pathway enrichment across samples ──────────────────
# Recurrent = significant (p.adjust < 0.05) in >= min_recur samples for a given
# (Module_group, geneset, pathway); top_n most-recurrent pathways per Module_group
# are plotted as a bar chart.
enrich_files <- list.files(deg_dir, pattern = "_enrich_modulegroup\\.Rds$", full.names = TRUE)
if (length(enrich_files) > 0) {
  enrich_long <- bind_rows(lapply(enrich_files, function(f) {
    s  <- sub("_enrich_modulegroup\\.Rds$", "", basename(f))
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

# ── Aggregate per-signature (pos vs neg) DEG genes across samples ─────────────
# Recurrent = significant (p_val_adj < 0.05) in >= min_recur samples for a given
# (signature, gene); top_n most-recurrent genes per signature are plotted as a
# bar chart (bar length = mean avg_log2FC, label = # samples recurrent).
deg_pos_files <- list.files(deg_dir, pattern = "_deg_.*_pos\\.Rds$", full.names = TRUE)
if (length(deg_pos_files) > 0) {
  deg_pos_long <- bind_rows(lapply(deg_pos_files, function(f) {
    s <- sub("^(.*)_deg_.*_pos\\.Rds$", "\\1", basename(f))
    readRDS(f) %>% mutate(slide = s)
  }))
  write.csv(deg_pos_long, file.path(path_dir, "all_deg_signature_pos_long.csv"), row.names = FALSE)

  min_recur_deg <- 2   # a gene must be significant in >=2 samples to count as "recurrent"
  top_n_deg     <- 15
  recur_deg <- deg_pos_long %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster, gene) %>%
    summarise(n_recur = dplyr::n_distinct(slide), mean_log2FC = mean(avg_log2FC), .groups = "drop") %>%
    filter(n_recur >= min_recur_deg)
  write.csv(recur_deg, file.path(path_dir, "recurrent_deg_genes_signature_pos.csv"), row.names = FALSE)

  for (lab in unique(recur_deg$cluster)) {
    d <- recur_deg %>% filter(cluster == lab) %>%
      arrange(desc(n_recur), desc(mean_log2FC)) %>% slice_head(n = top_n_deg) %>%
      mutate(gene = str_wrap(gsub("_", " ", gene), width = 20),
             gene = factor(gene, levels = rev(gene)))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(gene, mean_log2FC)) +
      geom_col(fill = "steelblue") +
      geom_text(aes(label = n_recur), hjust = -0.2, size = 3) +
      coord_flip() +
      labs(x = NULL, y = "Average log2FC",
           title = sprintf("Top recurrent DEG genes — %s+ (pos vs neg)", lab)) +
      theme_bw(base_size = 9)
    ggsave(file.path(path_dir, sprintf("top_deg_genes_%s_pos.png", gsub("/", "-", lab))), p,
           width = 8, height = max(3, nrow(d) * 0.35 + 1), dpi = 150, limitsize = FALSE)
  }
}

# ── Aggregate per-signature (pos vs neg) pathway enrichment across samples ────
# Recurrent = significant (p.adjust < 0.05) in >= min_recur samples for a given
# (signature, geneset, pathway); top_n most-recurrent pathways per signature
# are plotted as a bar chart.
enrich_pos_files <- list.files(deg_dir, pattern = "_enrich_.*_pos\\.Rds$", full.names = TRUE)
if (length(enrich_pos_files) > 0) {
  enrich_pos_long <- bind_rows(lapply(enrich_pos_files, function(f) {
    s   <- sub("_enrich_.*_pos\\.Rds$", "", basename(f))
    lab <- sub("^.*_enrich_(.*)_pos\\.Rds$", "\\1", basename(f))
    el  <- readRDS(f)
    bind_rows(lapply(names(el), function(gs) bind_rows(lapply(names(el[[gs]]), function(cl) {
      res <- el[[gs]][[cl]]
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      as.data.frame(res@result) %>% mutate(slide = s, geneset = gs, signature = lab)
    }))))
  }))
  write.csv(enrich_pos_long, file.path(path_dir, "all_pathway_enrichment_signature_pos_long.csv"),
            row.names = FALSE)

  min_recur_pos <- 2   # a pathway must be significant in >=2 samples to count as "recurrent"
  top_n_pos     <- 10
  recur_pos <- enrich_pos_long %>%
    filter(p.adjust < 0.05) %>%
    group_by(signature, geneset, ID) %>%
    summarise(n_recur = dplyr::n_distinct(slide), mean_NES = mean(NES), .groups = "drop") %>%
    filter(n_recur >= min_recur_pos)
  write.csv(recur_pos, file.path(path_dir, "recurrent_pathways_signature_pos.csv"), row.names = FALSE)

  for (lab in unique(recur_pos$signature)) {
    d <- recur_pos %>% filter(signature == lab) %>%
      arrange(desc(n_recur), desc(abs(mean_NES))) %>% slice_head(n = top_n_pos) %>%
      mutate(ID = factor(ID, levels = rev(ID)))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(ID, n_recur, fill = geneset)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = sprintf("# samples significant (of %d)", length(samples)),
           fill = "gene set", title = sprintf("Top recurrent enriched pathways — %s+ (pos vs neg)", lab)) +
      theme_bw(base_size = 9)
    ggsave(file.path(path_dir, sprintf("top_pathways_%s_pos.png", gsub("/", "-", lab))), p,
           width = 8, height = max(3, nrow(d) * 0.35 + 1), dpi = 150, limitsize = FALSE)
  }
}

message("\nDone. 6.2.3 signature analysis written to ", outdir)
