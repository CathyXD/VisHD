library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readxl)
library(purrr)
library(stringr)
library(msigdbr)
library(pals)
library(ggpubr)
library(qs2)

in_srt  <- path.expand("~/VisHD/1.integrate_raw_cell/normal_srt.qs2")
in_csv  <- path.expand("~/VisHD/9.normalcell_annotation/celltype_hint_per_cell.csv")
out_dir <- path.expand("~/VisHD/9.2.scimilarity_check")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading", in_srt, "\n")
srt  <- qs_read(in_srt)
cat("  ", ncol(srt), "cells\n")

hint <- read.csv(in_csv, row.names = 1, check.names = FALSE)
cat("Loaded", nrow(hint), "annotated cells from", in_csv, "\n")

srt$celltype_hint     <- hint[srt$cell_ID, "celltype_hint"]
srt$celltype_hint_raw <- hint[srt$cell_ID, "celltype_hint_raw"]
srt$min_dist          <- hint[srt$cell_ID, "min_dist"]

na_frac <- mean(is.na(srt$celltype_hint))
cat(sprintf("NA fraction after join: %.4f (%d / %d)\n",
            na_frac, sum(is.na(srt$celltype_hint)), ncol(srt)))

# 1 ── DimPlot on both reductions
dp_nb <- DimPlot(srt, reduction = "pearsonumap",
                 group.by = "celltype_hint", label = TRUE, label.size = 2, repel = TRUE,
                 cols = "polychrome") +
  ggtitle("SCimilarity hint (no batch correction)") +
  theme(legend.position = "none")
dp_b  <- DimPlot(srt, reduction = "pearsonbatchumap",
                 group.by = "celltype_hint", label = TRUE, label.size = 2, repel = TRUE,
                 cols = "polychrome") +
  ggtitle("SCimilarity hint (batch corrected)") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "1_celltype_DimPlot.png"),
       dp_nb | dp_b, width = 16, height = 7, dpi = 400)

# 2 ── Per-slide × category composition (stacked bar faceted by slide)
comp <- srt@meta.data %>%
  count(slide, category, celltype_hint) %>%
  group_by(slide, category) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(slide, category, celltype_hint, fill = list(n = 0, prop = 0))

bar <- ggplot(comp, aes(x = category, y = prop, fill = celltype_hint)) +
  geom_col() +
  facet_wrap(~ slide, nrow = 2) +
  labs(x = "Category", y = "Proportion",
       title = "Cell-type composition per slide × category") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 9))
ggsave(file.path(out_dir, "2_composition_bar.png"),
       bar, width = 16, height = 9, dpi = 400)

# 3 ── Boxplot: each cell type's proportion across the 8 slides (raw labels, no Unknown collapse)
comp_raw <- srt@meta.data %>%
  count(slide, celltype_hint_raw) %>%
  group_by(slide) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(slide, celltype_hint_raw, fill = list(n = 0, prop = 0))

box <- ggplot(comp_raw,
              aes(x = reorder(celltype_hint_raw, -prop, FUN = median), y = prop)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(aes(color = slide), width = 0.15, alpha = 0.6, size = 1) +
  labs(x = NULL, y = "Proportion (per slide)",
       title = "Per-slide variability of each cell type (celltype_hint_raw)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.key.size = unit(0.3, "cm"))
ggsave(file.path(out_dir, "3_celltype_boxplot.png"),
       box, width = 18, height = 6, dpi = 400)

write.csv(comp, file.path(out_dir, "composition.csv"), row.names = FALSE)

# ── TME cell-type annotation pipeline (mirror of 8.1) ─────────────────────
source("~/VisHD/normal_markers.R")
source("~/VisHD/functions.R")
source("~/VisHD/celltype_annotation_function.R")

tumour_markers <- c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")
SVEC_marker    <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                    "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")

meta_xlsx <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
meta_programs <- set_names(sheetname, sheetname) |>
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  map(~ map(.x, ~ as.character(na.omit(.x))))

# tme_markers: per-cell-type vectors with shared genes removed
tme_markers <- c(
  list(SVEC = SVEC_marker),
  lapply(meta_programs, function(progs) unique(unlist(progs)))
)
gene_counts  <- table(unlist(tme_markers))
shared_genes <- names(gene_counts[gene_counts > 1])
tme_markers  <- lapply(tme_markers, setdiff, shared_genes)
cat("tme_markers: dropped", length(shared_genes),
    "shared genes; per-cell-type sizes:",
    paste(names(tme_markers), lengths(tme_markers), sep = "=", collapse = ", "), "\n")

# C8 fallback markers — MSigDB regex per cell-type label
c8_data <- readRDS("~/VisHD/public_signature/c8_data_human.Rds")
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

# Annotation on pearson_clusters (no batch)
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

# Annotation on pearson_clusters_batch (batch corrected)
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
}

ggsave(file.path(out_dir, "10c_celltype_annotation_DimPlot_batch.png"),
       DimPlot(srt, group.by = "celltype_annotation_batch",
               cols = as.vector(polychrome()), reduction = "pearsonbatchumap") +
         ggtitle("TME cell-type annotation (batch corrected)"),
       width = 7, height = 5, dpi = 350)

# 11 ── Per-slide × category composition for TME annotation (batch)
comp_tme <- srt@meta.data %>%
  count(slide, category, celltype_annotation_batch) %>%
  group_by(slide, category) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(nesting(slide, category), celltype_annotation_batch, fill = list(n = 0, prop = 0))

n_tme   <- length(unique(comp_tme$celltype_annotation_batch))
tme_pal <- as.vector(polychrome())[seq_len(n_tme)]

bar_tme <- ggplot(comp_tme,
                  aes(x = category, y = prop, fill = celltype_annotation_batch)) +
  geom_col() +
  facet_grid(.~ slide,  scales = "free", space = "free") +
  scale_fill_manual(values = tme_pal, na.value = "grey80") +
  labs(x = "Category", y = "Proportion",
       title = "TME annotation (batch) — composition per slide × category") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        strip.text = element_text(size = 9))
ggsave(file.path(out_dir, "11_TME_composition_bar.png"),
       bar_tme, width = 16, height = 9, dpi = 400)
write.csv(comp_tme, file.path(out_dir, "TME_composition.csv"), row.names = FALSE)

# 11b ── Paired DT vs CB boxplot per cell type
comp_pair <- comp_tme %>%
  filter(category %in% c("DT", "CB 0", "CB 1", "CB")) %>%
  mutate(group = ifelse(grepl("^CB", category), "CB", "DT")) %>%
  group_by(slide, celltype_annotation_batch, group) %>%
  summarise(prop = sum(prop), .groups = "drop")

# paired wilcoxon per cell type (slide as the pairing unit)
comp_paired <- comp_pair %>%
  group_by(celltype_annotation_batch, slide) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("DT", "CB")))

box_pair <- ggpaired(comp_paired,
                     x = "group", y = "prop", id = "slide",
                     facet.by = "celltype_annotation_batch",
                     fill = "group", palette = c("#E69F00", "#56B4E9"),
                     line.color = "grey60", line.size = 0.3,
                     point.size = 1, short.panel.labs = TRUE) +
  stat_compare_means(method = "wilcox.test", paired = TRUE,
                     label = "p.format", size = 2.6) +
  labs(x = NULL, y = "Proportion (per slide × category)",
       title = "TME annotation (batch) — DT vs CB per cell type (paired Wilcoxon)") +
  theme(axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 7))
ggsave(file.path(out_dir, "11b_TME_DT_vs_CB_boxplot.png"),
       box_pair, width = 8, height = 8, dpi = 400)

# 12 ── Meta-program module scores (FeaturePlot per program, per cell-type subdir)
mp_root <- file.path(out_dir, "meta_program_scores")
dir.create(mp_root, showWarnings = FALSE, recursive = TRUE)

genes_in_srt <- rownames(srt)
for (ct in names(meta_programs)) {
  ct_dir <- file.path(mp_root, gsub("[^A-Za-z0-9]+", "_", ct))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  for (prog in names(meta_programs[[ct]])) {
    genes <- intersect(meta_programs[[ct]][[prog]], genes_in_srt)
    if (length(genes) < 3) next
    score_name <- "mp_score"
    srt <- AddModuleScore(srt, features = list(genes), name = score_name, assay = "Spatial")
    fp <- FeaturePlot(srt, features = paste0(score_name, "1"),
                      reduction = "pearsonbatchumap", order = TRUE) +
      scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, NA)) +
      ggtitle(sprintf("%s — %s (%d genes)", ct, prog, length(genes)))
    fname <- sprintf("%s.png", gsub("[^A-Za-z0-9]+", "_", prog))
    ggsave(file.path(ct_dir, fname), fp, width = 6, height = 5, dpi = 300)
    srt[[paste0(score_name, "1")]] <- NULL
  }
}

qs_save(srt, file.path(out_dir, "normal_srt_annotated.qs2"))

setwd(out_dir)
srt2anndata(srt, data_assay = "Spatial",
            save_name = "normal_srt_annotated", svg_path = NULL)

cat("Done. Outputs in", out_dir, "\n")
