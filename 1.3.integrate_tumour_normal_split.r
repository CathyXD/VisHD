library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
source("~/VisHD/functions.R")

# SCTransform/future export the count matrix as a global per worker; the default
# 500 MiB cap is too small for the integrated object.
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

in_dir  <- path.expand("~/VisHD/1.1.integrate_raw_cell")
full_path <- file.path(in_dir, "integrated_pearson_srt.qs2")
setwd(in_dir)
cat("Loading full integrated srt:", full_path, "\n")
srt <- qs_read(full_path )

ref_clusters <- as.character(c(0, 1, 12, 19, 4, 3, 5, 20, 9, 17))
svec_cluster <- "6"
normal_clusters <- c(ref_clusters, svec_cluster)

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% normal_clusters, "Normal", "Tumour")

g <- DimPlot(srt, group.by = "tumour_anno", label = TRUE) + 
scale_color_manual(values = c("Normal" = "lightgrey", "Tumour" = "red")) +
  labs(title = "Tumour vs normal cells") + 
  theme(legend.position = "none")
ggsave(filename = "7.tumour_vs_normal.png", plot = g, width = 6, height = 5, dpi = 300)

# ── Per-slide CB category from each sample's category.csv ────────────────
srt$category <- NA_character_
for (sl in unique(srt$slide)) {
  CB <- read.csv(file.path("~/VisHD", sl, "category.csv"))
  CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
  is_slide <- srt$slide == sl
  cellids   <- as.integer(gsub(paste0("^", sl, "_"), "",
                               colnames(srt)[is_slide]))
  srt$category[is_slide] <- CB$category[match(cellids, CB$cellid)]
}

g_cat <- DimPlot(srt, group.by = "category", label = TRUE, label.size = 3) +
  labs(title = "Category (CB / DT)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(filename = "7b.category.png", plot = g_cat, width = 6, height = 5, dpi = 300)

saveRDS(srt@meta.data, "tumour_anno_metadata.rds")

# tumour ========================================================================
out_dir <- file.path(in_dir, "tumour")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tumour_srt <- subset(srt, tumour_anno == "Tumour")
tumour_srt <- NormalizeData(tumour_srt)
tumour_srt <- do.pearson_pca(tumour_srt, find_hvgs = TRUE)
tumour_srt <- do.pearson_pca(tumour_srt, find_hvgs = TRUE, batch_variable = "slide")
# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(tumour_srt, reduction = "pearsonumap",
              group.by = "pearson_clusters", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (no batch)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

dp_batch <- DimPlot(tumour_srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (batch corrected)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2b_DimPlot_clusters_batch.png"), dp_batch,
       width = 6, height = 5, dpi = 400)
# ── Slide layout: UMAP (no batch) + UMAP (batch corrected) ────────────────
dp_s <- DimPlot(tumour_srt, reduction = "pearsonumap", group.by = "slide", label = FALSE) +
  ggtitle("UMAP (no batch)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

dp_s_batch <- DimPlot(tumour_srt, reduction = "pearsonbatchumap", group.by = "slide",
                      label = FALSE) +
  ggtitle("UMAP (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "2c_slide_layout.png"),
       (dp_s | dp_s_batch) + plot_annotation(title = "Integrated — cells by slide"),
       width = 12, height = 5, dpi = 400)

g_cat <- DimPlot(tumour_srt, reduction = "pearsonbatchumap",
                 group.by = "category", label = TRUE, label.size = 3) +
  labs(title = "Tumour — category (CB / DT)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(filename = file.path(out_dir, "7b.category.png"), plot = g_cat,
       width = 6, height = 5, dpi = 300)

qs_save(tumour_srt,  "tumour_srt.qs2")
srt2anndata <- srt2anndata(tumour_srt, data_assay = "Spatial", save_name = "tumour", svg_path = NULL)

# normal ========================================================================
out_dir <- file.path(in_dir, "normal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
normal_srt <- subset(srt, tumour_anno == "Normal")
normal_srt <- NormalizeData(normal_srt)
normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE)
normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE, batch_variable = "slide")
# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(normal_srt, reduction = "pearsonumap",
              group.by = "pearson_clusters", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (no batch)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

dp_batch <- DimPlot(normal_srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (batch corrected)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2b_DimPlot_clusters_batch.png"), dp_batch,
       width = 6, height = 5, dpi = 400)
# ── Slide layout: UMAP (no batch) + UMAP (batch corrected) ────────────────
dp_s <- DimPlot(normal_srt, reduction = "pearsonumap", group.by = "slide", label = FALSE) +
  ggtitle("UMAP (no batch)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

dp_s_batch <- DimPlot(normal_srt, reduction = "pearsonbatchumap", group.by = "slide",
                      label = FALSE) +
  ggtitle("UMAP (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "2c_slide_layout.png"),
       (dp_s | dp_s_batch) + plot_annotation(title = "Integrated — cells by slide"),
       width = 12, height = 5, dpi = 400)

g_cat <- DimPlot(normal_srt, reduction = "pearsonbatchumap",
                 group.by = "category", label = TRUE, label.size = 3) +
  labs(title = "Normal — category (CB / DT)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(filename = file.path(out_dir, "7b.category.png"), plot = g_cat,
       width = 6, height = 5, dpi = 300)


f <- FeaturePlot(normal_srt, features = c("tumour_score", "normal_score"),
                 reduction = "pearsonumap", ncol = 2) +
  theme(legend.key.size = unit(0.3, "cm"))
f_batch <- FeaturePlot(normal_srt, features = c("tumour_score", "normal_score"),
                 reduction = "pearsonbatchumap", ncol = 2) +
  theme(legend.key.size = unit(0.3, "cm"))

ggsave(file.path(out_dir, "3c_featureplot_tumour_normal_score.png"),
       (f /f_batch) + plot_annotation(title = "Normal — tumour/normal scores"),
       width = 8, height = 10, dpi = 400)
vln <- VlnPlot(normal_srt, features = c("tumour_score", "normal_score"),
               group.by = "pearson_clusters_batch") 

ggsave(file.path(out_dir, "3d_vlnplot_tumour_normal_score.png"), plot = vln) 
qs_save(normal_srt,  "normal_srt.qs2")
srt2anndata <- srt2anndata(normal_srt, data_assay = "Spatial", save_name = "normal", svg_path = NULL)

cat(dim(normal_srt)[2], "cells\n")
# second round of filtering to remove cluster 2, which has high tumour score and low normal score
normal_srt <- qs_read("normal_srt.qs2")
normal_srt <- subset(normal_srt, pearson_clusters_batch != "2")
out_dir <- file.path(in_dir, "normal_cleaned")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# normal_srt <- SCTransform(normal_srt, verbose = FALSE, vars.to.regress = "nCount_Spatial", assay = "Spatial")
# normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE, assay = "SCT")
# normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE, batch_variable = "slide", assay = "SCT")
normal_srt <- NormalizeData(normal_srt)
normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE)
normal_srt <- do.pearson_pca(normal_srt, find_hvgs = TRUE, batch_variable = "slide")

# ── Cluster DimPlots ───────────────────────────────────────────────────────
dp <- DimPlot(normal_srt, reduction = "pearsonumap",
              group.by = "pearson_clusters", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (no batch)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2a_DimPlot_clusters.png"), dp,
       width = 6, height = 5, dpi = 400)

dp_batch <- DimPlot(normal_srt, reduction = "pearsonbatchumap",
                    group.by = "pearson_clusters_batch", label = TRUE, label.size = 3) +
  ggtitle("Integrated — clusters (batch corrected)") +
  theme(legend.key.size = unit(0.4, "cm"))
ggsave(file.path(out_dir, "2b_DimPlot_clusters_batch.png"), dp_batch,
       width = 6, height = 5, dpi = 400)
# ── Slide layout: UMAP (no batch) + UMAP (batch corrected) ────────────────
dp_s <- DimPlot(normal_srt, reduction = "pearsonumap", group.by = "slide", label = FALSE) +
  ggtitle("UMAP (no batch)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

dp_s_batch <- DimPlot(normal_srt, reduction = "pearsonbatchumap", group.by = "slide",
                      label = FALSE) +
  ggtitle("UMAP (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 9))

ggsave(file.path(out_dir, "2c_slide_layout.png"),
       (dp_s | dp_s_batch) + plot_annotation(title = "Integrated — cells by slide"),
       width = 12, height = 5, dpi = 400)

g_cat <- DimPlot(normal_srt, reduction = "pearsonbatchumap",
                 group.by = "category", label = TRUE, label.size = 3) +
  labs(title = "Normal — category (CB / DT)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(filename = file.path(out_dir, "7b.category.png"), plot = g_cat,
       width = 6, height = 5, dpi = 300)


f <- FeaturePlot(normal_srt, features = c("tumour_score", "normal_score"),
                 reduction = "pearsonumap", ncol = 2) +
  theme(legend.key.size = unit(0.3, "cm"))
f_batch <- FeaturePlot(normal_srt, features = c("tumour_score", "normal_score"),
                 reduction = "pearsonbatchumap", ncol = 2) +
  theme(legend.key.size = unit(0.3, "cm"))

ggsave(file.path(out_dir, "3c_featureplot_tumour_normal_score.png"),
       (f /f_batch) + plot_annotation(title = "Normal — tumour/normal scores"),
       width = 8, height = 10, dpi = 400)
vln <- VlnPlot(normal_srt, features = c("tumour_score", "normal_score"),
               group.by = "pearson_clusters_batch") 

ggsave(file.path(out_dir, "3d_vlnplot_tumour_normal_score.png"), plot = vln) 
qs_save(normal_srt,  "normal_cleaned_srt.qs2")
srt2anndata <- srt2anndata(normal_srt, data_assay = "Spatial", save_name = "normal_cleaned", svg_path = NULL)

cat(dim(normal_srt)[2], "cells\n")