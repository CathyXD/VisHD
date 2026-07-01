#!/usr/bin/env Rscript
# 10.0.per_sample_tumour_normal.R   (per-sample, array 1-8)
# Combine clean tumour cells (tumour/tumour_srt.qs2) and finally-annotated
# normal cells (normal/normal_anno_srt.qs2) into one per-sample object, run
# SpaNorm → Pearson PCA → BANKSY clustering, add module scores via UCell for
# archetype / G1-G2-G3 / genesets2023 / tumour+normal signatures, then
# visualize across three embedding spaces plus bar-plot summaries.
#
# Output subdirs under LUT-245-XX/final_png/:
#   spanorm/    DimPlot + FeaturePlot on SpaNorm UMAP
#   pearsonpca/ DimPlot + FeaturePlot on Pearson UMAP + Gavish meta-programs
#   banksy/     DimPlot + FeaturePlot on BANKSY UMAP + composition + DEG
#   spatial/    ImageDimPlot + ImageFeaturePlot (tissue-image space)
#   barplots/   Composition bar plots grouped by cluster / category / subclone
#
#   Rscript 10.0.per_sample_tumour_normal.R <sample-index 1-8>

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
})

source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")   # all_marker

# ── Reference data ─────────────────────────────────────────────────────────────
Hall <- readRDS("~/VisHD/Hall.Rds")
C6   <- readRDS("~/VisHD/C6.Rds")
C5   <- readRDS("~/VisHD/C5.Rds")
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

clean_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

genesets2023 <- readRDS("~/VisHD/public_signature/genesets2023.Rds")
genesets2023 <- unlist(genesets2023, recursive = FALSE)


meta_xlsx     <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname     <- excel_sheets(meta_xlsx)
meta_programs <- set_names(sheetname, sheetname) %>%
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) %>%
  map(~ map(.x, ~ as.character(na.omit(.x))))
meta_programs_unlist <- unlist(meta_programs, recursive = F)
names(meta_programs_unlist) <- make.names(names(meta_programs_unlist), unique = T)

metas    <- readRDS("~/VisHD/6.2.3.signature_analysis/binarisation/metas.Rds")
groupdeg <- readRDS(paste0("~/VisHD/6.2archetype_downstream_tumour/archetype_module/",
                           "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))

tumour_markers <- c("AR", "FOLH1", "KLK2", "KLK3", "KLK4", "TMPRSS2",
                    "NKX3-1", "HOXB13", "TRPM8")
source("~/VisHD/normal_markers.R")
normal_modules <- unlist(all_marker)
# ── CLI arg ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 10.0.per_sample_tumour_normal.R <sample-index 1-8>")
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[arg]
i     <- basename(path)
cat("Working at", path, "\n")
setwd(path)

spanorm_dir  <- file.path(path, "final_png", "spanorm")
pearson_dir  <- file.path(path, "final_png", "pearsonpca")
banksy_dir   <- file.path(path, "final_png", "banksy")
spatial_dir  <- file.path(path, "final_png", "spatial")
barplots_dir <- file.path(path, "final_png", "barplots")
for (d in c(spanorm_dir, pearson_dir, banksy_dir, spatial_dir, barplots_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ── Palette helpers (needed in both load and processing paths) ─────────────────
group_pal <- c("Neg"      = "lightgrey",
               "G1"       = "red",
               "G2"       = "gold",
               "G3"       = "royalblue",
               "G1/G2"    = "orange",
               "G1/G3"    = "purple",
               "G2/G3"    = "green",
               "G1/G2/G3" = "black")
canon <- function(x) vapply(strsplit(x, "/"),
                            function(p) paste(sort(p), collapse = "/"), character(1))

if (file.exists("tumour_normal_anno_srt.qs2")) {
  cat("Found existing tumour_normal_anno_srt.qs2 — loading directly\n")
  srt            <- qs_read("tumour_normal_anno_srt.qs2")
  DEG            <- readRDS("tumour_normal_deg_spanorm.Rds")
} else {

# ════════════════════════════════════════════════════════════════════════════════
# SECTION A: PROCESSING
# ════════════════════════════════════════════════════════════════════════════════

# ── A1. Load and assemble combined object ──────────────────────────────────────
tumour_srt <- qs_read(file.path(path, "tumour", "tumour_srt.qs2"))
normal_srt <- qs_read(file.path(path, "normal", "normal_anno_srt.qs2"))
tum_cells  <- colnames(tumour_srt)
nor_cells  <- colnames(normal_srt)
cat("  tumour cells:", length(tum_cells), " normal cells:", length(nor_cells), "\n")

tum_label <- paste0("Tumour ", as.character(tumour_srt$subclone))
names(tum_label) <- tum_cells
nor_label <- as.character(normal_srt$final_annotation)
names(nor_label) <- nor_cells

parent <- qs_read(file.path(path, "tumour_subclone_srt.qs2"))
srt    <- subset(parent, cells = c(tum_cells, nor_cells))
rm(parent); gc()

srt$cell_type        <- unname(c(tum_label, nor_label)[colnames(srt)])
srt$final_annotation <- unname(nor_label[colnames(srt)])
srt$compartment      <- ifelse(colnames(srt) %in% tum_cells, "Tumour", "Normal")
cat("\nCombined object:", ncol(srt), "cells\n")
print(table(srt$cell_type, useNA = "ifany"))

# ── A2. SpaNorm + Pearson PCA ──────────────────────────────────────────────────
srt <- do.spanorm(srt)
srt <- do.pearson_pca(srt)
srt <- FindClusters(srt, resolution = 0.8, algorithm = 4)
qs_save(srt, "tumour_normal_anno_srt.qs2")
cat("SpaNorm + Pearson PCA done\n")

# ── A3. BANKSY clustering (lambda 0.2) ────────────────────────────────────────
srt <- SeuratWrappers::RunBanksy(srt, lambda = 0.2, verbose = TRUE,
                                 use_agf = TRUE, assay = "SpaNorm", slot = "data",
                                 k_geom = c(15), assay_name = "BANKSY_0.2")
srt <- RunPCA(srt, npcs = 30, features = rownames(srt), reduction.name = "banksy0.2.pca")
srt <- RunUMAP(srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name = "banksy0.2.umap")
srt <- FindNeighbors(srt, reduction = "banksy0.2.pca", dims = 1:20)
srt <- FindClusters(srt, resolution = 1, algorithm = 4)
srt$banksy_clusters <- srt$seurat_clusters
Idents(srt) <- "banksy_clusters"
qs_save(srt, "tumour_normal_anno_srt.qs2")
cat("BANKSY done\n")
# ── A8. DEG + enrichment (on banksy_clusters) ─────────────────────────────────
DefaultAssay(srt) <- "SpaNorm"
Idents(srt) <- "banksy_clusters"
DEG <- FindAllMarkers(srt, test.use = "MAST", only.pos = TRUE)
saveRDS(DEG, "tumour_normal_deg_spanorm.Rds")
run_gsea_panel(DEG, gene_sets, "tumour_normal_deg_enrich.Rds")
} # end if/else (cache check)
# ── A4. Clean meta.data (remove all stale score columns before scoring) ────────
{
  drop_score <- grep(
    "_UCell$|Module[0-9]*$|^ct_[0-9]+$|^GDmod_|^module_G|^gs23_|^arch_|^gd_|^tn_|ct_$|arch_$|gd_$|gs23_$|tn_$",
    colnames(srt@meta.data), value = TRUE
  )
  all_gset_names <- unique(c(names(clean_module), names(all_marker),
                              names(genesets2023), names(groupdeg), names(meta_programs)))
  drop_overlap <- intersect(all_gset_names, colnames(srt@meta.data))
  to_drop <- unique(c(drop_score, drop_overlap))
  if (length(to_drop) > 0) {
    srt@meta.data <- srt@meta.data[, setdiff(colnames(srt@meta.data), to_drop), drop = FALSE]
    cat("Dropped", length(to_drop), "stale score columns\n")
  }
}

# ── A5. AddModuleScore_UCell for all signatures ───────────────────────────────
# UCell with name="prefix_" creates prefix_1, prefix_2, ... (numbered).
# ucell_score() calls UCell then renames to prefix_<signame>_UCell.
ucell_score <- function(srt, feats, prefix) {
  feats <- lapply(feats, function(g) intersect(as.character(g), rownames(srt)))
  feats <- feats[lengths(feats) > 0]
  if (length(feats) == 0) return(list(srt = srt, cols = character(0)))
  srt      <- AddModuleScore_UCell(srt, features = feats, name = paste0(prefix, "_UCell"), slot = "data", missing_genes = "skip" )
  new_cols <- paste0(names(feats), prefix , "_UCell")# UCell names: <signame><prefix>_UCell
  list(srt = srt, cols = new_cols)
}

DefaultAssay(srt) <- "SpaNorm"

# 1. Archetype modules (clean_module) → arch_AR_UCell, arch_Inflammation_UCell, ...
res <- ucell_score(srt, clean_module, "_arch")
srt <- res$srt; arch_mod_cols <- res$cols

# 2. G1/G2/G3 group-DEG → gd_G1_UCell, gd_G2_UCell, gd_G3_UCell
res <- ucell_score(srt, groupdeg, "_gd")
srt <- res$srt; mod_score_cols <- res$cols

# 3. genesets2023 → gs23_*_UCell
res <- ucell_score(srt, genesets2023, "_gs23")
srt <- res$srt; gs23_cols <- res$cols

# 4. Tumour + normal scores → tn_tumour_score_UCell, tn_normal_score_UCell
tn_feats <- list(
  tumour_score = tumour_markers,
  normal_score  = unique(unlist(all_marker, use.names = FALSE))
)
res <- ucell_score(srt, tn_feats, "")
srt <- res$srt; tn_cols <- res$cols

# 5. Meta-programs (from 6.2.3 metas) → <meta_name>_meta_UCell
res <- ucell_score(srt, meta_programs_unlist, "_meta")
srt <- res$srt; meta_cols <- res$cols


cat("Module scores added. Meta.data columns:", ncol(srt@meta.data), "\n")

# ── A6. Module_group annotation (from 6.2.3 metas) ───────────────────────────
labs         <- names(groupdeg)
group_combos <- unlist(lapply(seq_along(labs), function(k)
  combn(labs, k, FUN = function(x) paste(x, collapse = "/"))))
group_levels <- c("Neg", group_combos)
metas_s   <- metas[metas$slide == i, ]
mg_lookup <- setNames(as.character(metas_s$Module_group), metas_s$cell)
is_tum    <- srt$compartment == "Tumour"
module_anno         <- rep("Normal", ncol(srt))
module_anno[is_tum] <- mg_lookup[colnames(srt)[is_tum]]
n_unmatched <- sum(is_tum & is.na(module_anno))
cat(sprintf("\nTumour module_anno: %d/%d matched (%.1f%%); %d -> 'Neg'\n",
            sum(is_tum) - n_unmatched, sum(is_tum),
            100 * (sum(is_tum) - n_unmatched) / sum(is_tum), n_unmatched))
module_anno[is.na(module_anno)] <- "Neg"

present <- levels(droplevels(factor(module_anno)))
mg_pal  <- setNames(group_pal[canon(present)], present)
mg_pal["Normal"] <- "lightpink"

mg_levels       <- c("Normal", group_levels)
srt$module_anno <- factor(module_anno, levels = mg_levels[mg_levels %in% module_anno])
cat("\nmodule_anno:\n"); print(table(srt$module_anno, useNA = "ifany"))

# ── A7. Fill final_annotation NAs with module_anno ────────────────────────────
fa <- as.character(srt$final_annotation)
fa[is.na(fa)] <- paste("Tumour", srt$module_anno[is.na(fa)])
srt$final_annotation <- fa
cat("\nfinal_annotation:\n"); print(table(srt$final_annotation, useNA = "ifany"))

# ── A9. Save ───────────────────────────────────────────────────────────────────
qs_save(srt, "tumour_normal_anno_srt.qs2")
cat("==================== ", i, " processing done ====================\n")



# ════════════════════════════════════════════════════════════════════════════════
# SECTION B: VISUALIZATION
# ════════════════════════════════════════════════════════════════════════════════

ct_lvls <- sort(unique(as.character(srt$cell_type)))
ct_pal  <- setNames(as.vector(polychrome())[seq_along(ct_lvls)], ct_lvls)

# Combined palette for final_annotation (normal cell types + G-group labels)
final_anno_pal <- c(ct_pal, mg_pal)
final_anno_pal <- final_anno_pal[!duplicated(names(final_anno_pal))]

grad2 <- scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                 midpoint = 0)

# ── B0. Helper: FeaturePlot with auto-scaled height ───────────────────────────
fp_save <- function(srt, feats, red, outfile, ncol = 3, dpi = 400) {
  feats <- feats[feats %in% colnames(srt@meta.data)]
  if (length(feats) == 0) return(invisible(NULL))
  nrows <- ceiling(length(feats) / ncol)
  p <- FeaturePlot(srt, reduction = red, features = feats, ncol = ncol, order = TRUE)
  ggsave(outfile, p & grad2,
         width = ncol * 4, height = nrows * 3.5,
         dpi = dpi, limitsize = FALSE)
}

# ── B0. Helper: composition bar plots ─────────────────────────────────────────
barplot_comp <- function(meta, group_var, fill_var, fill_pal, outfile) {
  if (!group_var %in% colnames(meta) || !fill_var %in% colnames(meta))
    return(invisible(NULL))
  d <- meta[!is.na(meta[[group_var]]) & !is.na(meta[[fill_var]]), , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  lvls_g <- sort(unique(as.character(d[[group_var]])))
  lvls_f <- sort(unique(as.character(d[[fill_var]])))
  comp <- as.data.frame(table(
    group = factor(as.character(d[[group_var]]), levels = lvls_g),
    fill  = factor(as.character(d[[fill_var]]),  levels = lvls_f)
  ))
  pal <- if (!is.null(fill_pal) && length(fill_pal) > 0) {
    fill_pal
  } else {
    setNames(as.vector(polychrome())[seq_along(lvls_f)], lvls_f)
  }
  base <- ggplot(comp, aes(group, Freq, fill = fill)) +
    scale_fill_manual(values = pal, na.value = "grey80", name = fill_var) +
    labs(x = group_var) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p_count <- base + geom_col() +
    labs(y = "cells", title = paste(i, "—", group_var, "×", fill_var))
  p_prop <- base + geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) + labs(y = "proportion")
  w <- max(8, length(lvls_g) * 0.5 + 3)
  ggsave(outfile, p_count / p_prop, width = w, height = 10, dpi = 150, limitsize = FALSE)
}

# ── B1. SpaNorm UMAP ───────────────────────────────────────────────────────────
cat("SpaNorm plots...\n")

p_ct   <- DimPlot(srt, reduction = "umap", group.by = "cell_type",
                  cols = ct_pal, raster = FALSE) +
            ggtitle(paste(i, "— cell type"))
p_comp <- DimPlot(srt, reduction = "umap", group.by = "compartment",
                  raster = FALSE) +
            ggtitle(paste(i, "— compartment"))
p_cl   <- DimPlot(srt, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE, label.size = 3, repel = TRUE,
                  cols = as.vector(polychrome()), raster = FALSE) +
            ggtitle(sprintf("%s — SpaNorm clusters (n=%d)", i,
                            nlevels(factor(srt$seurat_clusters)))) +
            theme(legend.position = "none")
p_mg   <- DimPlot(srt, reduction = "umap", group.by = "module_anno",
                  cols = mg_pal, raster = FALSE) +
            ggtitle(paste(i, "— module group")) +
            theme(legend.text = element_text(size = 7))
p_fa   <- DimPlot(srt, reduction = "umap", group.by = "final_annotation",
                  cols = final_anno_pal, raster = FALSE) +
            ggtitle(paste(i, "— final annotation")) +
            theme(legend.text = element_text(size = 7))
p_sub  <- DimPlot(srt, reduction = "umap", group.by = "subclone",
                  label = TRUE, label.size = 3, repel = TRUE,
                  cols = as.vector(polychrome()), raster = FALSE) +
            ggtitle(paste(i, "— subclone"))

ggsave(file.path(spanorm_dir, "1_DimPlot_combined.png"),
       (p_ct | p_comp | p_cl) / (p_mg | p_fa | p_sub),
       width = 27, height = 14, dpi = 150)



fp_save(srt, arch_mod_cols,  "umap", file.path(spanorm_dir, "7_FeaturePlot_archetype_modules.png"))
fp_save(srt, mod_score_cols, "umap", file.path(spanorm_dir, "8_FeaturePlot_G123_scores.png"), ncol = 3)
fp_save(srt, tn_cols,        "umap", file.path(spanorm_dir, "9_FeaturePlot_tumour_normal_scores.png"), ncol = 2)

if (length(gs23_cols) > 0)
  fp_save(srt, gs23_cols, "umap",
          file.path(spanorm_dir, "10_FeaturePlot_genesets2023.png"), ncol = 4)

# ── B2. Pearson PCA UMAP ───────────────────────────────────────────────────────
cat("Pearson PCA plots...\n")

ggsave(file.path(pearson_dir, "1_DimPlot_cell_type.png"),
       DimPlot(srt, reduction = "pearsonumap", group.by = "cell_type",
               cols = ct_pal, raster = FALSE) +
         ggtitle(paste(i, "— cell type")),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(pearson_dir, "2_DimPlot_compartment.png"),
       DimPlot(srt, reduction = "pearsonumap", group.by = "compartment",
               raster = FALSE) +
         ggtitle(paste(i, "— compartment")),
       width = 8, height = 7, dpi = 150)

ggsave(file.path(pearson_dir, "3_DimPlot_clusters.png"),
       DimPlot(srt, reduction = "pearsonumap", group.by = "pearson_clusters",
               label = TRUE, label.size = 3, repel = TRUE,
               cols = as.vector(polychrome()), raster = FALSE) +
         ggtitle(sprintf("%s — Pearson clusters (n=%d)", i,
                         nlevels(factor(srt$pearson_clusters)))) +
         theme(legend.position = "none"),
       width = 8, height = 7, dpi = 150)

ggsave(file.path(pearson_dir, "4_DimPlot_module_anno.png"),
       DimPlot(srt, reduction = "pearsonumap", group.by = "module_anno",
               cols = mg_pal, raster = FALSE) +
         ggtitle(paste(i, "— module group")) +
         theme(legend.text = element_text(size = 7)),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(pearson_dir, "5_DimPlot_final_annotation.png"),
       DimPlot(srt, reduction = "pearsonumap", group.by = "final_annotation",
               cols = final_anno_pal, raster = FALSE) +
         ggtitle(paste(i, "— final annotation")) +
         theme(legend.text = element_text(size = 7)),
       width = 10, height = 7, dpi = 150)

ggsave(file.path(pearson_dir, "6_FeaturePlot_tumour_markers.png"),
       FeaturePlot(srt, reduction = "pearsonumap", features = tumour_markers,
                   ncol = 3, order = TRUE),
       width = 12, height = 12, dpi = 150)

fp_save(srt, arch_mod_cols,  "pearsonumap", file.path(pearson_dir, "7_FeaturePlot_archetype_modules.png"))
fp_save(srt, mod_score_cols, "pearsonumap", file.path(pearson_dir, "8_FeaturePlot_G123_scores.png"), ncol = 3)
fp_save(srt, tn_cols,        "pearsonumap", file.path(pearson_dir, "9_FeaturePlot_tumour_normal_scores.png"), ncol = 2)

if (length(gs23_cols) > 0)
  fp_save(srt, gs23_cols, "pearsonumap",
          file.path(pearson_dir, "10_FeaturePlot_genesets2023.png"), ncol = 4)

# Gavish meta-program FeaturePlots across all three reductions (saves to pearson_dir)
srt <- score_meta_programs(srt, meta_programs, out_dir = pearson_dir,
                            reductions = c("umap", "pearsonumap", "banksy0.2.umap"))

# ── B3. BANKSY UMAP ────────────────────────────────────────────────────────────
cat("BANKSY plots...\n")

ggsave(file.path(banksy_dir, "1_DimPlot_banksy_clusters.png"),
       DimPlot(srt, reduction = "banksy0.2.umap", group.by = "banksy_clusters",
               label = TRUE, label.size = 3, repel = TRUE,
               cols = as.vector(polychrome()), raster = FALSE) +
         ggtitle(sprintf("%s — BANKSY clusters (n=%d)", i,
                         nlevels(factor(srt$banksy_clusters)))) +
         theme(legend.position = "none"),
       width = 8, height = 7, dpi = 150)

ggsave(file.path(banksy_dir, "2_DimPlot_cell_type.png"),
       DimPlot(srt, reduction = "banksy0.2.umap", group.by = "cell_type",
               cols = ct_pal, raster = FALSE) +
         ggtitle(paste(i, "— cell type")),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(banksy_dir, "3_DimPlot_module_anno.png"),
       DimPlot(srt, reduction = "banksy0.2.umap", group.by = "module_anno",
               cols = mg_pal, raster = FALSE) +
         ggtitle(paste(i, "— module group")) +
         theme(legend.text = element_text(size = 7)),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(banksy_dir, "4_DimPlot_final_annotation.png"),
       DimPlot(srt, reduction = "banksy0.2.umap", group.by = "final_annotation",
               cols = final_anno_pal, raster = FALSE) +
         ggtitle(paste(i, "— final annotation")) +
         theme(legend.text = element_text(size = 7)),
       width = 10, height = 7, dpi = 150)

ggsave(file.path(banksy_dir, "5_FeaturePlot_tumour_markers.png"),
       FeaturePlot(srt, reduction = "banksy0.2.umap", features = tumour_markers,
                   ncol = 3, order = TRUE),
       width = 12, height = 12, dpi = 150)

fp_save(srt, arch_mod_cols,  "banksy0.2.umap", file.path(banksy_dir, "6_FeaturePlot_archetype_modules.png"))
fp_save(srt, mod_score_cols, "banksy0.2.umap", file.path(banksy_dir, "7_FeaturePlot_G123_scores.png"), ncol = 3)
fp_save(srt, tn_cols,        "banksy0.2.umap", file.path(banksy_dir, "8_FeaturePlot_tumour_normal_scores.png"), ncol = 2)

if (length(gs23_cols) > 0)
  fp_save(srt, gs23_cols, "banksy0.2.umap",
          file.path(banksy_dir, "9_FeaturePlot_genesets2023.png"), ncol = 4)

# Cluster × cell-type composition
comp <- as.data.frame(table(cluster = srt$banksy_clusters, cell_type = srt$cell_type))
p_count <- ggplot(comp, aes(cluster, Freq, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = ct_pal, name = "Cell type") +
  labs(title = paste(i, "— BANKSY cluster composition"),
       x = "BANKSY cluster", y = "cells") +
  theme_bw(base_size = 11)
p_prop <- ggplot(comp, aes(cluster, Freq, fill = cell_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = ct_pal, name = "Cell type") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "BANKSY cluster", y = "proportion") +
  theme_bw(base_size = 11)
ggsave(file.path(banksy_dir, "10_composition_barplot.png"),
       p_count / p_prop, width = 12, height = 10, dpi = 150)
write.csv(comp, file.path(banksy_dir, "banksy_cluster_celltype_composition.csv"),
          row.names = FALSE)

plot_top_deg(srt, DEG %>% filter(p_val_adj < 0.05),
             out_dir = banksy_dir, prefix = "11_DEG")

# ── B4. Spatial — tissue-image space ──────────────────────────────────────────
cat("Spatial plots...\n")

ggsave(file.path(spatial_dir, "1_ImageDimPlot_cell_type.png"),
       ImageDimPlot(srt, group.by = "cell_type", cols = ct_pal,
                    border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(spatial_dir, "2_ImageDimPlot_compartment.png"),
       ImageDimPlot(srt, group.by = "compartment",
                    border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(spatial_dir, "3_ImageDimPlot_banksy_clusters.png"),
       ImageDimPlot(srt, group.by = "banksy_clusters", cols = "polychrome",
                    border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(spatial_dir, "4_ImageDimPlot_module_anno.png"),
       ImageDimPlot(srt, group.by = "module_anno", cols = mg_pal,
                    border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(spatial_dir, "5_ImageDimPlot_final_annotation.png"),
       ImageDimPlot(srt, group.by = "final_annotation", cols = final_anno_pal,
                    border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 150)

ggsave(file.path(spatial_dir, "6_ImageFeaturePlot_tumour_markers.png"),
       ImageFeaturePlot(srt, features = tumour_markers, cols = c("white", "red")),
       width = 15, height = 12, dpi = 150)

ggsave(file.path(spatial_dir, "7_ImageFeaturePlot_archetype_modules.png"),
       ImageFeaturePlot(srt, features = arch_mod_cols),
       width = 15, height = 10, dpi = 150)

ggsave(file.path(spatial_dir, "8_ImageFeaturePlot_G123_scores.png"),
       ImageFeaturePlot(srt, features = mod_score_cols),
       width = 12, height = 5, dpi = 150)

ggsave(file.path(spatial_dir, "9_ImageFeaturePlot_tumour_normal_scores.png"),
       ImageFeaturePlot(srt, features = tn_cols),
       width = 10, height = 5, dpi = 150)

if (length(gs23_cols) > 0) {
  gs23_ncol <- min(4L, length(gs23_cols))
  gs23_nrow <- ceiling(length(gs23_cols) / gs23_ncol)
  ggsave(file.path(spatial_dir, "10_ImageFeaturePlot_genesets2023.png"),
         ImageFeaturePlot(srt, features = gs23_cols),
         width = gs23_ncol * 4, height = gs23_nrow * 3,
         dpi = 150, limitsize = FALSE)
}

# ── B5. Bar plots (grouping × cell_type and module_anno) ──────────────────────
cat("Bar plots...\n")
meta <- srt@meta.data

grouping_vars <- c("seurat_clusters", "banksy_clusters", "pearson_clusters",
                   "category", "subclone")

for (gvar in grouping_vars) {
  barplot_comp(meta, gvar, "cell_type",   ct_pal,
               file.path(barplots_dir, paste0("barplot_", gvar, "_by_cell_type.png")))
  barplot_comp(meta, gvar, "module_anno", mg_pal,
               file.path(barplots_dir, paste0("barplot_", gvar, "_by_module_anno.png")))
}

cat("==================== ", i, " visualization done ====================\n")
