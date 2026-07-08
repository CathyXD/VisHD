#!/usr/bin/env Rscript
# 8.2.DEG_cluster_annotation.R
# Re-embed the annotated normal-cell object (9.2 output) with SCT-derived HVGs +
# batch-corrected Pearson PCA (by slide), cluster at res 1.0, run MAST DE between
# clusters, overlay the prior annotations on the NEW pearsonbatchumap, and
# quantify how each cluster's significant DEGs overlap the meta-program gene sets.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(purrr)
  library(stringr)
  library(pals)
  library(qs2)
})
source("~/VisHD/functions.R")   # do.pearson_pca()
stopifnot(requireNamespace("MAST", quietly = TRUE))      # FindAllMarkers(test.use = "MAST")
stopifnot(requireNamespace("glmGamPoi", quietly = TRUE)) # fast SCTransform
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB; SCTransform future workers exceed the 500 MiB default

# ── Paths ─────────────────────────────────────────────────────────────────────
in_srt    <- path.expand("~/VisHD/8.1.scimilarity_check/normal_srt_annotated.qs2")
clu_csv   <- path.expand("~/VisHD/9.normalcell_annotation/cellsearch_cluster/cell_scimilarity_cluster_celltype.csv")
meta_xlsx <- path.expand("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx")
out_dir   <- path.expand("~/VisHD/8.2.DEG_cluster_annotation")
png_dir   <- file.path(out_dir, "png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

RES <- 1.0   # clustering resolution on the new batch-corrected Pearson graph

# ── 1. Load object + join the SCimilarity cluster-cellsearch label ────────────
cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
cat("  ", ncol(srt), "cells x", nrow(srt), "genes\n")

clu <- read.csv(clu_csv, check.names = FALSE)
srt$scim_cluster_celltype <- clu$scimilarity_cluster_celltype[match(srt$cell_ID, clu$cell_ID)]
cat("Joined scim cluster-cellsearch:", sum(!is.na(srt$scim_cluster_celltype)),
    "/", ncol(srt), "cells matched\n")

# ── 2. Simplify the SCimilarity single-cell hint into broad (first-layer) types ─
simplify_broad <- function(x) {
  xl  <- tolower(as.character(x))
  out <- rep("Other", length(xl))
  out[grepl("fibroblast", xl)]                                                   <- "Fibroblasts"
  out[grepl("myofibroblast|smooth muscle|pericyte|mural", xl)]                   <- "Smooth muscle"
  out[grepl("endotheli", xl)]                                                    <- "Endothelial"
  out[grepl("epitheli|luminal|basal cell|secretory|club cell|goblet|keratinocyte|acinar|principal cell|tubule|collecting duct|loop of henle|glandular", xl)] <- "Epithelial"
  out[grepl("macrophage|monocyte|dendritic|kupffer|microglia|myeloid|neutrophil|granulocyte", xl)] <- "Myeloid"
  out[grepl("\\bb cell|plasma cell", xl)]                                        <- "B cells"
  out[grepl("\\bt cell|thymocyte|cd4|cd8|regulatory t|gamma-delta", xl)]         <- "T cells"
  out[grepl("natural killer|innate lymphoid", xl)]                               <- "NK/ILC"
  out[grepl("mast cell", xl)]                                                    <- "Mast"
  out[grepl("glial|neuron|schwann|astrocyte|oligodendro", xl)]                   <- "Glial/Neural"
  out[grepl("muscle cell|cardiac|myocyte", xl) & out == "Other"]                 <- "Muscle"
  out[grepl("unknown", xl)]                                                      <- "Unknown"
  out[is.na(x)]                                                                  <- "Unknown"
  out
}
srt$scim_broad <- simplify_broad(srt$celltype_hint)
cat("\nscim_broad (simplified single-cell hint):\n")
print(table(srt$scim_broad, useNA = "ifany"))

# ── 3. SCTransform -> HVGs on the SCT assay ───────────────────────────────────
DefaultAssay(srt) <- "Spatial"
srt <- SCTransform(srt, assay = "Spatial", method = "glmGamPoi",
                   variable.features.n = 3000, verbose = FALSE)
hvg_sct <- VariableFeatures(srt, assay = "SCT")
# Pearson PCA runs on raw Spatial counts; keep only HVGs present there.
spatial_genes <- rownames(GetAssayData(srt, assay = "Spatial", layer = "counts"))
hvg <- intersect(hvg_sct, spatial_genes)
VariableFeatures(srt, assay = "Spatial") <- hvg
cat("\nSCT HVGs:", length(hvg_sct), "-> usable on Spatial:", length(hvg), "\n")

# ── 4. Batch-corrected Pearson PCA (by slide) + new UMAP + clusters ───────────
# Overwrites pearsonbatchpca / pearsonbatchumap / pearsonbatchgraph and
# pearson_clusters_batch with the SCT-HVG embedding.
DefaultAssay(srt) <- "Spatial"
srt <- do.pearson_pca(srt, batch_variable = "slide", assay = "Spatial",
                      find_hvgs = FALSE, reduction_prefix = "pearsonbatch",
                      clusters_col = "pearson_clusters_batch", resolution = RES)
Idents(srt) <- "pearson_clusters_batch"
n_clu <- nlevels(factor(srt$pearson_clusters_batch))
cat("New batch-corrected Pearson clustering (res", RES, "):", n_clu, "clusters\n")

# ── 5. DE between clusters with MAST (up-regulated markers) ───────────────────
DefaultAssay(srt) <- "Spatial"
srt <- NormalizeData(srt, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(srt, assay = "Spatial", test.use = "MAST",
                          only.pos = TRUE, logfc.threshold = 0.25,
                          min.pct = 0.1, max.cells.per.ident = 1000,
                          verbose = FALSE)
saveRDS(markers, file.path(out_dir, "FindAllMarkers_MAST.Rds"))
write.csv(markers, file.path(out_dir, "FindAllMarkers_MAST.csv"), row.names = FALSE)

sig_deg <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)
write.csv(sig_deg, file.path(out_dir, "significant_DEGs.csv"), row.names = FALSE)
cat("Significant up-DEGs (padj<0.05, log2FC>0.25):", nrow(sig_deg),
    "across", dplyr::n_distinct(sig_deg$cluster), "clusters\n")

# ── 6. Visualise the new clustering + prior annotations on pearsonbatchumap ───
# "Unknown" is always drawn in light grey; other levels use a polychrome palette.
plot_anno <- function(group, title, label = TRUE) {
  v  <- as.character(srt@meta.data[[group]]); v[is.na(v)] <- "Unknown"
  srt@meta.data[[group]] <- v
  others <- setdiff(sort(unique(v)), "Unknown")
  pal <- setNames(as.vector(polychrome())[seq_along(others)], others)
  if ("Unknown" %in% v) pal <- c(pal, Unknown = "lightgrey")
  DimPlot(srt, reduction = "pearsonbatchumap", group.by = group, cols = pal,
          label = label, label.size = 2, repel = TRUE) +
    ggtitle(title) + theme(legend.text = element_text(size = 6))
}

p_clu <- DimPlot(srt, reduction = "pearsonbatchumap", group.by = "pearson_clusters_batch",
                 label = TRUE, label.size = 3, repel = TRUE, cols = as.vector(polychrome())) +
  ggtitle(sprintf("Pearson batch clusters (res %.1f, n=%d)", RES, n_clu)) +
  theme(legend.position = "none")
ggsave(file.path(png_dir, "1_clusters_pearsonbatchumap.png"), p_clu,
       width = 8, height = 7, dpi = 300)

p_scim_s <- plot_anno("scim_broad",            "SCimilarity single cellsearch (broad)")
p_scim_c <- plot_anno("scim_cluster_celltype", "SCimilarity cluster cellsearch")
p_tme    <- plot_anno("celltype_annotation_batch", "TME annotation pipeline (batch)")

ggsave(file.path(png_dir, "2_scimilarity_single_broad.png"), p_scim_s, width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "3_scimilarity_cluster.png"),      p_scim_c, width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "4_TME_annotation.png"),           p_tme,    width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "5_annotation_overlay_grid.png"),
       (p_clu | p_scim_s) / (p_scim_c | p_tme), width = 18, height = 14, dpi = 300)

# ── 7. Meta-program overlap with cluster DEGs ─────────────────────────────────
# First layer = sheet (cell type); programs = columns within each sheet.
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
meta_programs <- set_names(sheetname, sheetname) |>
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  map(~ map(.x, ~ as.character(na.omit(.x))))

deg_by_cluster <- split(sig_deg$gene, sig_deg$cluster)
overlap_df <- map_dfr(names(deg_by_cluster), function(cl) {
  degs <- unique(deg_by_cluster[[cl]])
  map_dfr(names(meta_programs), function(layer) {
    map_dfr(names(meta_programs[[layer]]), function(prog) {
      pg <- meta_programs[[layer]][[prog]]
      pg <- unique(pg[pg != "" & !is.na(pg)])
      data.frame(cluster = cl, first_layer = layer, program = prog,
                 n_deg = length(degs), n_program_genes = length(pg),
                 n_overlap = length(intersect(degs, pg)),
                 stringsAsFactors = FALSE)
    })
  })
})
write.csv(overlap_df, file.path(out_dir, "cluster_metaprogram_overlap.csv"), row.names = FALSE)

# Per first-layer summary (sum of overlaps across that layer's programs)
layer_summary <- overlap_df %>%
  group_by(cluster, first_layer) %>%
  summarise(n_overlap = sum(n_overlap), .groups = "drop")
write.csv(layer_summary, file.path(out_dir, "cluster_firstlayer_overlap_summary.csv"),
          row.names = FALSE)

# ── 8. Combined faceted bar plot: facet by first layer, x = program, fill = cluster
clu_lvls <- sort(unique(overlap_df$cluster))
clu_cols <- setNames(as.vector(polychrome())[seq_along(clu_lvls)], clu_lvls)

bar <- ggplot(overlap_df, aes(x = program, y = n_overlap, fill = factor(cluster))) +
  geom_col(position = position_dodge(width = 0.85)) +
  facet_wrap(~ first_layer, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = clu_cols, name = "Cluster") +
  labs(title = "Cluster DEG overlap with meta-programs (first layer = cell type)",
       x = "Meta-program", y = "# overlapping DEGs") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.key.size = unit(0.3, "cm"), strip.text = element_text(size = 9))
ggsave(file.path(png_dir, "6_metaprogram_overlap_bar.png"),
       bar, width = 22, height = 14, dpi = 200, limitsize = FALSE)

# ── 9. Save re-embedded + annotated object ────────────────────────────────────
qs_save(srt, file.path(out_dir, "normal_srt_DEcluster.qs2"))
cat("\nDone. Outputs in", out_dir, "\n")
