#!/usr/bin/env Rscript
# 8.7.integrate_normal_anno.R
# Rebuild the cross-sample integrated normal-cell object (run-once). The
# previous integration (7.1.normal_cell_integration/integrated_pearson_srt2.qs2)
# was lost to a Pawsey scratch purge. Rather than re-deriving annotations from
# scratch via 7.1 -> 8.1 -> 8.2 -> 8.4, this merges the 8 already-annotated
# per-sample LUT-245-XX/normal/normal_anno_srt.qs2 objects (final_annotation
# resolved by 8.6.final_clear_normal_persample.R) and re-embeds them jointly
# with Pearson residual PCA, no-batch and batch-corrected by slide.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pals)
  library(qs2)
  library(SeuratWrappers)
})
source("~/VisHD/functions.R")   # do.pearson_pca()
options(future.globals.maxSize = 50 * 1024^3)

RES <- 1.0   # clustering resolution on the batch-corrected Pearson graph

# ── Paths ──────────────────────────────────────────────────────────────────
paths   <- system("realpath ~/VisHD/LUT-245-*/normal/normal_anno_srt.qs2", intern = TRUE)
slides  <- basename(gsub("/normal", "", dirname(paths)))
out_dir <- path.expand("~/VisHD/8.4.final_clear_normal_integration")
png_dir <- file.path(out_dir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
integrated_file <- file.path(out_dir, "normal_srt_final_anno.qs2")

cat("Found", length(paths), "slides:", paste(slides, collapse = ", "), "\n")

srt <- NULL
if (file.exists(integrated_file)) {
  cat("Loading precomputed integrated srt from:", integrated_file, "\n")
  srt <- qs_read(integrated_file)
  if (!"integrated.rpca_umap" %in% Reductions(srt)) {
    cat("Cached object is missing the integrated.rpca_umap reduction -- rebuilding.\n")
    srt <- NULL
  }
}

if (is.null(srt)) {
  # ── Load and merge counts + metadata ────────────────────────────────────
  cat("Loading and merging", length(paths), "slides...\n")
  srt_list <- lapply(seq_along(paths), function(i) {
    cat("  Loading", slides[i], "\n")
    srt_full <- qs_read(paths[i])
    counts   <- GetAssayData(srt_full, assay = "Spatial", layer = "counts")
    meta     <- srt_full@meta.data
    meta$final_annotation <- as.character(meta$final_annotation)
    meta$slide <- slides[i]
    coords   <- GetTissueCoordinates(srt_full, which = "centroids")
    meta$x_centroid <- coords[rownames(meta), "x"]
    meta$y_centroid <- coords[rownames(meta), "y"]
    srt_mini <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
    rm(srt_full); gc()
    srt_mini
  })

  srt <- merge(srt_list[[1]], y = srt_list[-1], add.cell.ids = slides)
  rm(srt_list); gc()
  cat("Merged:", ncol(srt), "cells\n")

  # merge() already splits layers per input object (counts.1..8); join them
  # back and re-split by slide so layer names match `slide` instead of the
  # positional counts.N naming, which was mislabeling one SCT model's
  # umi.assay as "RNA" instead of "Spatial" and breaking
  # PrepSCTFindMarkers() downstream.
  srt[["Spatial"]] <- JoinLayers(srt[["Spatial"]])
  srt[["Spatial"]] <- split(srt[["Spatial"]], f = srt$slide)

  # ── 3. SCTransform -> HVGs on the SCT assay (mirror 9.3) ──────────────────────
  # NB: keep the per-slide split Spatial layers intact through IntegrateLayers
  # -- RPCAIntegration needs them. Only JoinLayers afterwards.
DefaultAssay(srt) <- "Spatial"
srt <- SCTransform(srt, assay = "Spatial", method = "glmGamPoi",
                   variable.features.n = 3000, verbose = FALSE)

# SCTransform() on a split-layer assay only patches umi.assay on
# SCTModel.list[[1]] (Seurat v5 SCTransform.Seurat); the other per-slide
# models keep an internal default of "RNA", which later trips
# PrepSCTFindMarkers's "Multiple UMI assays are used for SCTransform" check.
# Force every model's umi.assay so all 8 agree.
sct_assay <- srt[["SCT"]]
sct_models <- slot(sct_assay, "SCTModel.list")
for (m in names(sct_models)) {
  slot(sct_models[[m]], "umi.assay") <- "Spatial"
}
slot(sct_assay, "SCTModel.list") <- sct_models
srt[["SCT"]] <- sct_assay

srt <- RunPCA(srt, npcs = 30, verbose = F)
srt <- IntegrateLayers(object = srt, method = RPCAIntegration, normalization.method = "SCT",new.reduction = "integrated.rpca", verbose = F)
srt <- FindNeighbors(srt, dims = 1:30, reduction = "integrated.rpca")
srt <- FindClusters(srt, resolution = 1, cluster.name = "integrated.rpca_clusters", verbose = F)
srt <- RunUMAP(srt, dims = 1:30, reduction = "integrated.rpca", verbose = F, reduction.name = "integrated.rpca_umap")

srt <- JoinLayers(srt, assay = "Spatial", layers = c("counts"),
                   new.layer.names = c("counts"))
DefaultAssay(srt) <- "Spatial"                   
srt <- SCTransform(srt, assay = "Spatial", method = "glmGamPoi",
                   variable.features.n = 3000, verbose = FALSE)
hvg_sct <- VariableFeatures(srt, assay = "SCT")
spatial_genes <- rownames(GetAssayData(srt, assay = "Spatial", layer = "counts"))
hvg <- intersect(hvg_sct, spatial_genes)
VariableFeatures(srt, assay = "Spatial") <- hvg
cat("\nSCT HVGs:", length(hvg_sct), "-> usable on Spatial:", length(hvg), "\n")

# ── 4. Batch-corrected Pearson PCA (by slide) + new UMAP + clusters (mirror 9.3)
DefaultAssay(srt) <- "Spatial"
srt <- do.pearson_pca(srt, batch_variable = "slide", assay = "Spatial",
                      find_hvgs = FALSE, reduction_prefix = "pearsonbatch",
                      clusters_col = "pearson_clusters_batch", resolution = RES)
Idents(srt) <- "pearson_clusters_batch"
n_clu <- nlevels(factor(srt$pearson_clusters_batch))
cat("New batch-corrected Pearson clustering (res", RES, "):", n_clu, "clusters\n")

# ── 5. DE between clusters with MAST (up-regulated markers) (mirror 9.3) ──────
DefaultAssay(srt) <- "SCT"
srt <- PrepSCTFindMarkers(srt, assay = "SCT")
markers <- FindAllMarkers(srt, assay = "SCT", test.use = "MAST",
                          only.pos = TRUE, logfc.threshold = 0.25,
                          min.pct = 0.1, max.cells.per.ident = 1000,
                          verbose = FALSE)
saveRDS(markers, file.path(out_dir, "FindAllMarkers_MAST.Rds"))
write.csv(markers, file.path(out_dir, "FindAllMarkers_MAST.csv"), row.names = FALSE)

sig_deg <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)
write.csv(sig_deg, file.path(out_dir, "significant_DEGs.csv"), row.names = FALSE)
cat("Significant up-DEGs (padj<0.05, log2FC>0.25):", nrow(sig_deg),
    "across", dplyr::n_distinct(sig_deg$cluster), "clusters\n")


  qs_save(srt, integrated_file)
}

# ── Sanity plots ─────────────────────────────────────────────────────────
final_anno_levels <- sort(unique(srt$final_annotation))
final_pal  <- setNames(as.vector(polychrome())[seq_along(final_anno_levels)],
                        final_anno_levels)
slide_cols <- as.vector(brewer.set1(length(unique(srt$slide))))


dp_anno_batch <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "final_annotation",
                          cols = final_pal, label = TRUE, label.size = 3) +
  ggtitle("Integrated normal cells — final_annotation (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "1b_DimPlot_final_annotation_batch.png"), dp_anno_batch,
       width = 7, height = 5, dpi = 350)


dp_clusters_batch <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "pearson_clusters_batch",
                              cols = as.vector(polychrome()), label = TRUE, label.size = 3) +
  ggtitle("Integrated normal cells — pearson_clusters_batch (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "2b_DimPlot_clusters_batch.png"), dp_clusters_batch,
       width = 7, height = 5, dpi = 350)

dp_slide_batch <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "slide", cols = slide_cols) +
  ggtitle("Integrated normal cells — by slide (batch corrected)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "3_DimPlot_slide.png"), dp_slide_batch,
       width = 7, height = 5, dpi = 350)

dp_rpca_clusters <- DimPlot(srt, reduction = "integrated.rpca_umap", group.by = "integrated.rpca_clusters",
                             cols = as.vector(polychrome()), label = TRUE, label.size = 3) +
  ggtitle("Integrated normal cells — integrated.rpca_clusters (RPCA)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "4a_DimPlot_rpca_clusters.png"), dp_rpca_clusters,
       width = 7, height = 5, dpi = 350)

dp_rpca_anno <- DimPlot(srt, reduction = "integrated.rpca_umap", group.by = "final_annotation",
                         cols = final_pal, label = TRUE, label.size = 3) +
  ggtitle("Integrated normal cells — final_annotation (RPCA)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "4b_DimPlot_rpca_final_annotation.png"), dp_rpca_anno,
       width = 7, height = 5, dpi = 350)

dp_rpca_slide <- DimPlot(srt, reduction = "integrated.rpca_umap", group.by = "slide", cols = slide_cols) +
  ggtitle("Integrated normal cells — by slide (RPCA)") +
  theme(legend.key.size = unit(0.3, "cm"))
ggsave(file.path(png_dir, "4c_DimPlot_rpca_slide.png"), dp_rpca_slide,
       width = 7, height = 5, dpi = 350)

cat("Done. Integrated object at", integrated_file, "\n")
