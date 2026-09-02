# 6.4.DT_signature_analysis_aggregate.R  (run-once)
# Aggregates the per-sample outputs written by 6.4.DT_signature_analysis.r
# (array 1-8: meta.Rds, sp_group.qs2, sp_pos.qs2, module_group_deg/enrich_modulegroup.Rds
# per sample) into cross-sample boxplots, positivity bar plots, Module_group
# composition, per-cell Wilcoxon tests, spatial grids, and recurrent pathway
# enrichment. Run after all 8 array tasks have completed.

suppressPackageStartupMessages({
  library(tidyverse)
  library(qs2)
  library(patchwork)
  library(ggpubr)
  library(RColorBrewer)
  library(data.table)
})

samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir <- "~/VisHD"
groupdeg <- readRDS(file.path(base_dir,
  "6.3.DT_archetype_module",
  "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))
module_names <- paste0("module_", names(groupdeg))
banksy_col   <- "BANKSY_0.2_snn_res.1"

outdir   <- file.path(base_dir, "6.4.DT_signature_analysis")
box_dir  <- file.path(outdir, "boxplots")
sts_dir  <- file.path(outdir, "stats")
spt_dir  <- file.path(outdir, "spatial")
path_dir <- file.path(outdir, "pathway_enrichment")
corr_dir <- file.path(outdir, "correlation")
for (d in c(box_dir, sts_dir, spt_dir, path_dir, corr_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# same fixed Module_group levels/palette as the per-sample script
group_combos <- unlist(lapply(seq_along(module_names), function(k)
  combn(sub("^module_", "", module_names), k, FUN = function(x) paste(x, collapse = "/"))))
group_levels <- c("Neg", group_combos)
single_combos <- group_combos[!grepl("/", group_combos)]        # one signature positive
multi_combos  <- group_combos[grepl("/", group_combos)]         # >= 2 signatures positive
group_pal     <- c(Neg = "lightgrey",
                    setNames(brewer.pal(max(3, length(single_combos)), "Set1")[seq_along(single_combos)],
                             single_combos),
                    setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(multi_combos)),
                             multi_combos))[group_levels]


# ── Load per-sample outputs ────────────────────────────────────────────────────
meta_list <- list()
sp_pos    <- setNames(vector("list", length(module_names)), module_names)  # ImageDimPlots
sp_group  <- list()                                                        # Module_group maps
for (s in samples) {
  cat(s, "\n")
  meta_f <- file.path(outdir, s, "meta.Rds")
  # if (!file.exists(meta_f)) { message("Skipping ", s, " — no meta.Rds (not yet processed)"); next }
  meta_list[[s]] <- readRDS(meta_f)
  cat("read in meta\n")
  sp_group_f <- file.path(outdir, s, "sp_group.qs2")
  if (file.exists(sp_group_f)) sp_group[[s]] <- qs_read(sp_group_f)
  cat("read in sp_group\n")
  sp_pos_f <- file.path(outdir, s, "sp_pos.qs2")
  if (file.exists(sp_pos_f)) {
    sp_pos_s <- qs_read(sp_pos_f)
    for (m in names(sp_pos_s)) sp_pos[[m]][[s]] <- sp_pos_s[[m]]
  }
  cat("read in sp_pos\n")
}
if (!length(meta_list)) stop("No per-sample outputs found — run 6.4.DT_signature_analysis.r (array 1-8) first")

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
  pos_pal <- c(neg = "grey85", pos = "indianred")

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

  min_recur_pos <- 1   # a pathway must be significant in >=2 samples to count as "recurrent"
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

message("\nDone. 6.4.DT_signature_analysis_aggregate written to ", outdir)
