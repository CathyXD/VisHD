#!/usr/bin/env Rscript
# 10.1.per_sample_tumour_normal.R   (per-sample, array 1-8)
# Combine the clean tumour cells (tumour/tumour_srt.qs2) and the finally-annotated
# normal cells (normal/normal_anno_srt.qs2) of one sample into a single object,
# re-run the 4.1/4.2 workflow (SpaNorm -> Pearson PCA -> BANKSY clustering) and
# marker-expression panels, then add a per-BANKSY-cluster cell-type composition
# bar chart. Tumour cells are labelled "Tumour <subclone>"; normal cells keep
# their `final_annotation`. Output object: LUT-245-XX/tumour_normal_anno_srt.qs2;
# all images to LUT-245-XX/final_png/.
#
#   Rscript 10.1.per_sample_tumour_normal.R <sample-index 1-8>
#
# The combined cell SET + annotations come from the two requested objects, but
# counts + the single FOV come from their shared parent tumour_subclone_srt.qs2
# (both were subset from it) — merging the two processed objects directly would
# leave two overlapping per-sample FOVs and break ImageDimPlot / do.spanorm.

library(Seurat)
library(dplyr)
library(SpaNorm,     lib.loc = "~/R_Library/4.5")
library(qs2)
library(leidenbase,  lib.loc = "~/R_Library/4.5")
library(UCell,       lib.loc = "~/R_Library/4.5")
library(ggplot2)
library(readxl)
library(purrr)
library(patchwork)
library(stringr)
library(pals)

source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")   # all_marker

# ── Reference data (mirror 4.1/4.2) ─────────────────────────────────────────
Hall <- readRDS("~/VisHD/Hall.Rds")
C6   <- readRDS("~/VisHD/C6.Rds")
C5   <- readRDS("~/VisHD/C5.Rds")
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

clean_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

meta_xlsx     <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname     <- excel_sheets(meta_xlsx)
meta_programs <- set_names(sheetname, sheetname) %>%
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) %>%
  map(~ map(.x, ~ as.character(na.omit(.x))))

SVEC_marker <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                 "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")

tumour_markers <- c("AR", "FOLH1", "KLK2", "KLK3", "KLK4", "TMPRSS2",
                    "NKX3-1", "HOXB13", "TRPM8")

# ── CLI arg ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 10.1.per_sample_tumour_normal.R <sample-index 1-8>")
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[arg]
i     <- basename(path)
cat("working at", path, "\n")
setwd(path)
dir.create("final_png", showWarnings = FALSE, recursive = TRUE)

# ── 1. Load the two requested objects + build the combined annotation ─────────
tumour_srt <- qs_read(file.path(path, "tumour", "tumour_srt.qs2"))
normal_srt <- qs_read(file.path(path, "normal", "normal_anno_srt.qs2"))
tum_cells  <- colnames(tumour_srt)
nor_cells  <- colnames(normal_srt)
cat("  tumour cells:", length(tum_cells), " normal cells:", length(nor_cells), "\n")

# cell_type: tumour -> "Tumour <subclone>"; normal -> final_annotation
tum_label <- paste0("Tumour ", as.character(tumour_srt$subclone))
names(tum_label) <- tum_cells
nor_label <- as.character(normal_srt$final_annotation)
names(nor_label) <- nor_cells

# ── 2. Rebuild the combined object from the shared single-FOV parent ──────────
parent <- qs_read(file.path(path, "tumour_subclone_srt.qs2"))
keep   <- c(tum_cells, nor_cells)                 # disjoint (Tumour vs Normal)
srt    <- subset(parent, cells = keep)
rm(parent); gc()

# unname() so Seurat's [[<- assigns by position — a partially-matching named
# vector (e.g. nor_label has no tumour barcodes) yields NA names and trips
# Seurat's `all(names(value) == colnames(x))` check.
srt$cell_type        <- unname(c(tum_label, nor_label)[colnames(srt)])
srt$final_annotation <- unname(nor_label[colnames(srt)])  # NA on tumour cells
srt$compartment      <- ifelse(colnames(srt) %in% tum_cells, "Tumour", "Normal")
cat("\nCombined object:", ncol(srt), "cells\n")
print(table(srt$cell_type, useNA = "ifany"))

# ── 3. SpaNorm + Pearson PCA + clusters (mirror 4.1/4.2) ──────────────────────
srt <- do.spanorm(srt)
srt <- do.pearson_pca(srt)
srt <- FindClusters(srt, resolution = 0.8, algorithm = 4)
spatial_plot(srt, outdir = "final_png/", name = "spanorm")
qs_save(srt, "tumour_normal_anno_srt.qs2")
cat("SpaNorm Done\n")

# ── 4. BANKSY clustering (lambda 0.2) (mirror 4.1/4.2) ────────────────────────
srt <- SeuratWrappers::RunBanksy(srt, lambda = 0.2, verbose = TRUE,
                                 use_agf = TRUE, assay = "SpaNorm", slot = "data",
                                 k_geom = c(15), assay_name = "BANKSY_0.2")
srt <- RunPCA(srt, npcs = 30, features = rownames(srt), reduction.name = "banksy0.2.pca")
srt <- RunUMAP(srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name = "banksy0.2.umap")
srt <- FindNeighbors(srt, reduction = "banksy0.2.pca", dims = 1:20)
srt <- FindClusters(srt, resolution = 1, algorithm = 4)
srt$banksy_clusters <- srt$seurat_clusters          # freeze the BANKSY clusters
Idents(srt) <- "banksy_clusters"
qs_save(srt, "tumour_normal_anno_srt.qs2")
spatial_plot(srt, outdir = "final_png/", name = "Banksy_lam0.2")
cat("BANKSY DONE\n")

# ── 5. Marker-expression panels (mirror 4.1 archetype + 4.2 normal scores) ────
# Archetype module scores (tumour programs).
srt <- AddModuleScore(srt, features = clean_module, name = "Module")
colnames(srt@meta.data)[grep("Module", colnames(srt@meta.data), fixed = TRUE)] <-
  paste0(names(clean_module), "_Module")
mod_cols <- paste0(names(clean_module), "_Module")
g <- FeaturePlot(srt, mod_cols, ncol = 3, reduction = "banksy0.2.umap") &
     scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred")
ggsave("final_png/archetype_module_exp.pdf", plot = g, width = 9, height = 6)

# Normal cell-type module scores.
# Drop signatures whose genes are absent from this sample's probe panel —
# AddModuleScore aborts on any list left with no features present.
all_marker <- lapply(all_marker, function(g) intersect(g, rownames(srt)))
all_marker <- all_marker[lengths(all_marker) > 0]
srt <- AddModuleScore(srt, features = all_marker, name = "ct_")
added <- paste0("ct_", seq_along(all_marker))
colnames(srt@meta.data)[match(added, colnames(srt@meta.data))] <- names(all_marker)
ggsave("final_png/normal_score_FeaturePlot.png",
       plot = FeaturePlot(srt, names(all_marker), cols = c("white", "red"),
                          reduction = "banksy0.2.umap"),
       width = 25, height = 15, dpi = 350)
ggsave("final_png/normal_score_ImageFeaturePlot.png",
       plot = ImageFeaturePlot(srt, names(all_marker), cols = c("white", "red")),
       width = 25, height = 15, dpi = 350)

# Gavish meta-programs + lineage markers.
srt <- score_meta_programs(srt, meta_programs, out_dir = "final_png")
ggsave("final_png/SVEC_marker.png",
       plot = FeaturePlot(srt, SVEC_marker, reduction = "banksy0.2.umap"),
       width = 12, height = 12)
ggsave("final_png/tumour_marker.png",
       plot = FeaturePlot(srt, tumour_markers, reduction = "banksy0.2.umap"),
       width = 12, height = 12)
ggsave("final_png/tumour_marker_ImageFeaturePlot.png",
       plot = ImageFeaturePlot(srt, tumour_markers), width = 15, height = 12)

# ── 6. Per-BANKSY-cluster cell-type composition bar chart ─────────────────────
comp <- as.data.frame(table(cluster   = srt$banksy_clusters,
                            cell_type = srt$cell_type))
ct_lvls <- sort(unique(as.character(comp$cell_type)))
ct_pal  <- setNames(as.vector(polychrome())[seq_along(ct_lvls)], ct_lvls)

p_count <- ggplot(comp, aes(cluster, Freq, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = ct_pal, name = "Cell type") +
  labs(title = paste(i, "— BANKSY cluster composition (counts)"),
       x = "BANKSY cluster", y = "cells") +
  theme_bw(base_size = 11)
p_prop <- ggplot(comp, aes(cluster, Freq, fill = cell_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = ct_pal, name = "Cell type") +
  labs(title = "Proportion", x = "BANKSY cluster", y = "proportion") +
  theme_bw(base_size = 11)
ggsave("final_png/banksy_cluster_celltype_composition.png",
       plot = p_count / p_prop, width = 12, height = 10, dpi = 350)

# Annotation overlays for context.
ggsave("final_png/cell_type_ImageDimPlot.png",
       plot = ImageDimPlot(srt, group.by = "cell_type", cols = ct_pal,
                           border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 350)
ggsave("final_png/banksy_clusters_ImageDimPlot.png",
       plot = ImageDimPlot(srt, group.by = "banksy_clusters", cols = "polychrome",
                           border.color = "#00000000", size = 0.3),
       width = 9, height = 7, dpi = 350)

write.csv(comp, "final_png/banksy_cluster_celltype_composition.csv", row.names = FALSE)
cat("Wrote composition bar chart + csv\n")

# ── 7. Cluster DEG (mirror 4.1/4.2 workflow) ──────────────────────────────────
DefaultAssay(srt) <- "SpaNorm"
DEG <- FindAllMarkers(srt, test.use = "MAST", only.pos = TRUE)
saveRDS(DEG, "tumour_normal_deg_spanorm.Rds")
run_gsea_panel(DEG, gene_sets, "tumour_normal_deg_enrich.Rds")
plot_top_deg(srt, DEG %>% filter(p_val_adj < 0.05),
             out_dir = "final_png", prefix = "spanorm_DEG")

qs_save(srt, "tumour_normal_anno_srt.qs2")
cat("=====================", i, "tumour+normal done ====================\n")
