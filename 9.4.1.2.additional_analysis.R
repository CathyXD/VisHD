#!/usr/bin/env Rscript
# 9.4.1.2.additional_analysis.R
# Extra analysis on the final normal annotation. Mirrors the descriptive /
# visual half of 9.2.scimilarity_check.R, but keyed on `final_annotation`
# (joined from 9.4.1's final_annotation.csv by cell_ID) instead of the
# SCimilarity hint, and on the 9.4.1 re-embedded object.
#
# Loads:
#   9.4.1.final_clear_normal_integration/normal_srt_final_anno.qs2
#   9.4.1.final_clear_normal_integration/final_annotation.csv   (cell_ID -> final_annotation)
# Writes (to 9.4.1 out_dir/final_annotation_analysis/):
#   1_final_annotation_DimPlot.png      DimPlot on pearsonbatchumap
#   2_composition_bar.png               per-slide x category composition
#   3_final_annotation_boxplot.png      per-slide proportion of each annotation
#   3b_DT_vs_CB_boxplot.png             DT vs CB per cell type (paired Wilcoxon)
#   composition.csv
#   meta_program_scores/<celltype>/<program>.png        (AddModuleScore)
#   meta_program_scores/<celltype>_UCell/<program>.png  (UCell)
#
# Skipped vs 9.2 (annotation is already final): TME annotation pipeline, the
# DT-vs-CB paired Wilcoxon, and the annotated qs2 / AnnData re-export.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(purrr)
  library(stringr)
  library(pals)
  library(ggpubr)
  library(qs2)
  library(UCell, lib.loc = "~/R_Library/4.5")
})

out_dir <- path.expand("~/VisHD/9.4.1.final_clear_normal_integration")
in_srt  <- file.path(out_dir, "normal_srt_final_anno.qs2")
in_csv  <- file.path(out_dir, "final_annotation.csv")
viz_dir <- file.path(out_dir, "final_annotation_analysis")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
cat("  ", ncol(srt), "cells\n")

# ── Join final_annotation from the CSV by cell_ID (mirror of 9.2's hint join) ─
anno <- read.csv(in_csv, check.names = FALSE, stringsAsFactors = FALSE)
rownames(anno) <- anno$cell_ID
cat("Loaded", nrow(anno), "annotated cells from", in_csv, "\n")

srt$final_annotation <- anno[srt$cell_ID, "final_annotation"]
srt$source_cluster   <- anno[srt$cell_ID, "source_cluster"]
na_frac <- mean(is.na(srt$final_annotation))
cat(sprintf("NA fraction after join: %.4f (%d / %d)\n",
            na_frac, sum(is.na(srt$final_annotation)), ncol(srt)))

# ── 1. DimPlot of the final annotation on the batch-corrected UMAP ────────────
dp <- DimPlot(srt, reduction = "pearsonbatchumap",
              group.by = "final_annotation", label = TRUE, label.size = 3,
              repel = TRUE, cols = as.vector(polychrome())) +
  ggtitle("Final normal annotation (batch corrected)") +
  theme(legend.text = element_text(size = 7))
ggsave(file.path(viz_dir, "1_final_annotation_DimPlot.png"),
       dp, width = 9, height = 7, dpi = 400)

# ── 2. Per-slide x category composition (stacked bar faceted by slide) ────────
has_category <- "category" %in% colnames(srt@meta.data)
if (has_category) {
  comp <- srt@meta.data %>%
    mutate(cat2 = ifelse(startsWith(as.character(category), "CB"), "CB", "DT")) %>%
    count(slide, cat2, final_annotation) %>%
    group_by(slide, cat2) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(slide, cat2, final_annotation, fill = list(n = 0, prop = 0))
  bar <- ggplot(comp, aes(x = cat2, y = prop, fill = final_annotation)) +
    geom_col() +
    facet_grid(~ slide) +
    labs(x = "Category", y = "Proportion",
         title = "Final annotation composition per slide x category")
} else {
  comp <- srt@meta.data %>%
    count(slide, final_annotation) %>%
    group_by(slide) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    complete(slide, final_annotation, fill = list(n = 0, prop = 0))
  bar <- ggplot(comp, aes(x = slide, y = prop, fill = final_annotation)) +
    geom_col() +
    labs(x = "Slide", y = "Proportion",
         title = "Final annotation composition per slide")
}
bar <- bar +
  scale_fill_manual(values = as.vector(polychrome())) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 9))
ggsave(file.path(viz_dir, "2_composition_bar.png"), bar, width = 16, height = 9, dpi = 400)
write.csv(comp, file.path(viz_dir, "composition.csv"), row.names = FALSE)

# ── 3. Boxplot: each annotation's proportion across the 8 slides ──────────────
comp_slide <- srt@meta.data %>%
  count(slide, final_annotation) %>%
  group_by(slide) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(slide, final_annotation, fill = list(n = 0, prop = 0))

box <- ggplot(comp_slide,
              aes(x = reorder(final_annotation, -prop, FUN = median), y = prop)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(aes(color = slide), width = 0.15, alpha = 0.6, size = 1) +
  labs(x = NULL, y = "Proportion (per slide)",
       title = "Per-slide variability of each final annotation") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.key.size = unit(0.3, "cm"))
ggsave(file.path(viz_dir, "3_final_annotation_boxplot.png"),
       box, width = 12, height = 6, dpi = 400)

# ── 3b. Paired DT vs CB boxplot per cell type (slide is the pairing unit) ─────
# Uses the per-slide x cat2 proportions from section 2 (CB = CB 0 + CB 1).
if (has_category) {
  comp_paired <- comp %>%
    group_by(final_annotation, slide) %>%
    filter(n() == 2) %>%                       # keep slides with both DT and CB
    ungroup() %>%
    mutate(group = factor(cat2, levels = c("DT", "CB")))

  box_pair <- ggpaired(comp_paired,
                       x = "group", y = "prop", id = "slide",
                       facet.by = "final_annotation",
                       fill = "group", palette = c("#E69F00", "#56B4E9"),
                       line.color = "grey60", line.size = 0.3,
                       point.size = 1, short.panel.labs = TRUE) +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       label = "p.format", size = 2.6) +
    labs(x = NULL, y = "Proportion (per slide x category)",
         title = "Final annotation — DT vs CB per cell type (paired Wilcoxon)") +
    theme(axis.text.x = element_text(size = 8),
          strip.text = element_text(size = 7))
  ggsave(file.path(viz_dir, "3b_DT_vs_CB_boxplot.png"),
         box_pair, width = 9, height = 8, dpi = 400)
}

# ── 4. Meta-program module scores (FeaturePlot per program, per cell-type) ────
meta_xlsx <- path.expand("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx")
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
meta_programs <- set_names(sheetname, sheetname) |>
  purrr::map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  purrr::map(~ purrr::map(.x, ~ as.character(na.omit(.x))))

mp_root <- file.path(viz_dir, "meta_program_scores")
dir.create(mp_root, showWarnings = FALSE, recursive = TRUE)
genes_in_srt <- rownames(srt)

# AddModuleScore (one program at a time, blue gradient)
for (ct in names(meta_programs)) {
  ct_dir <- file.path(mp_root, gsub("[^A-Za-z0-9]+", "_", ct))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  for (prog in names(meta_programs[[ct]])) {
    genes <- intersect(meta_programs[[ct]][[prog]], genes_in_srt)
    if (length(genes) < 3) next
    srt <- AddModuleScore(srt, features = list(genes), name = "mp_score", assay = "Spatial")
    fp <- FeaturePlot(srt, features = "mp_score1",
                      reduction = "pearsonbatchumap", order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, NA)) +
      ggtitle(sprintf("%s - %s (%d genes)", ct, prog, length(genes)))
    ggsave(file.path(ct_dir, sprintf("%s.png", gsub("[^A-Za-z0-9]+", "_", prog))),
           fp, width = 6, height = 5, dpi = 300)
    srt[["mp_score1"]] <- NULL
  }
}

# UCell (whole sheet at once, red gradient)
for (ct in names(meta_programs)) {
  ct_dir <- file.path(mp_root, paste(gsub("[^A-Za-z0-9]+", "_", ct), "UCell", sep = "_"))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  srt <- AddModuleScore_UCell(srt, features = meta_programs[[ct]], missing_genes = "skip")
  for (prog in names(meta_programs[[ct]])) {
    genes <- intersect(meta_programs[[ct]][[prog]], genes_in_srt)
    fp <- FeaturePlot(srt, features = paste0(prog, "_UCell"),
                      reduction = "pearsonbatchumap", order = TRUE) +
      scale_color_gradient(low = "snow", high = "red", limits = c(0, NA)) +
      ggtitle(sprintf("%s - %s (%d genes)", ct, prog, length(genes)))
    ggsave(file.path(ct_dir, sprintf("%s.png", gsub("[^A-Za-z0-9]+", "_", prog))),
           fp, width = 6, height = 5, dpi = 300)
  }
}



cat("\nDone. Outputs in", viz_dir, "\n")
