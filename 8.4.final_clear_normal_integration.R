#!/usr/bin/env Rscript
# 8.4.final_clear_normal_integration.R
# Integrative stage of the final normal-cell cleanup. Loads the 9.3 re-embedded
# object, freezes the curated cluster -> cell-type map (interpreted from the 9.3
# MAST DEG + meta-program overlap on `pearson_clusters_batch`) as the FINAL
# annotation, DROPS the discarded clusters ("Removed"), then mirrors the
# 8.2.DEG_cluster_annotation.R workflow on the cleaned object: SCT-derived HVGs
# + batch-corrected Pearson PCA (by slide), cluster at res 1.0, MAST DE between
# clusters, annotation overlays on the new pearsonbatchumap, and cluster-DEG /
# meta-program overlap. Exports the per-cell final annotation as CSV (the join
# key for 9.4.2) and the re-embedded object.

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
options(future.globals.maxSize = 8 * 1024^3)  # SCTransform future workers exceed the 500 MiB default

# ── Curated cluster -> final annotation map (on the 9.3 pearson_clusters_batch) ─
# Clusters 8/10/11 are discarded (stressed / stress-heatshock / tumour contam).
cluster_map <- c(
  "0"  = "Smooth muscle",
  "1"  = "CAF",
  "2"  = "Pericyte",
  "3"  = "Macrophages",
  "4"  = "B cells",
  "5"  = "Epithelial",
  "6"  = "SVEC",
  "7"  = "Smooth muscle",
  "8"  = "Removed",
  "9"  = "SVEC",
  "10" = "Removed",
  "11" = "Removed",
  "12" = "Neurons",
  "13" = "Plasma"
)

# ── Paths ─────────────────────────────────────────────────────────────────────
in_srt    <- path.expand("~/VisHD/8.2.DEG_cluster_annotation/normal_srt_DEcluster.qs2")
meta_xlsx <- path.expand("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx")
out_dir   <- path.expand("~/VisHD/8.4.final_clear_normal_integration")
png_dir   <- file.path(out_dir, "png")
anno_csv  <- file.path(out_dir, "final_annotation.csv")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

RES <- 1.0   # clustering resolution on the new batch-corrected Pearson graph

# ── 1. Load object + apply the curated final annotation ───────────────────────
cat("Loading", in_srt, "\n")
srt <- qs_read(in_srt)
cat("  ", ncol(srt), "cells x", nrow(srt), "genes\n")

src_clu  <- as.character(srt$pearson_clusters_batch)   # basis of the curated map
unmapped <- setdiff(sort(unique(src_clu)), names(cluster_map))
if (length(unmapped))
  stop("pearson_clusters_batch has unmapped levels: ", paste(unmapped, collapse = ", "),
       " — update cluster_map at the top of 8.4.final_clear_normal_integration.R")

srt$source_cluster   <- src_clu                        # 9.3 cluster (kept for the CSV)
srt$final_annotation <- unname(cluster_map[src_clu])
cat("\nFinal annotation (before clearing Removed):\n")
print(table(srt$final_annotation, useNA = "ifany"))

# ── Representative positive markers: prostate / Wolffian epithelial populations ──

tumour_markers <- c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")

mature_luminal   <- c("KLK3", "KLK2", "ACPP", "MSMB",
                      "NKX3-1", "FOLH1", "AR", "TMPRSS2")

luminal_prog     <- c("PIGR", "SCGB1A1", "SCGB3A1", "MMP7", "LTF",
                      "OLFM4", "CP", "PSCA", "LCN2",
                      "KRT13", "KRT4", "KRT19")   # last 3 = hillock-leaning

basal_cell       <- c("TP63", "KRT5", "KRT14", "KRT15", "KRT17", "NGFR")

epididymal       <- c("CRISP1", "PATE1", "SPINK2", "WFDC2", "LCN5",
                      "LCN6", "GPX5", "RNASE10", "ADGRG2", "CTSE")

seminal_vesicle  <- c("SEMG1", "SEMG2", "PATE4", "MUC6", "PGC",
                      "PAEP", "SERPINA5", "CYP4F8")

# ── Combined named list (for AddModuleScore / UCell / DotPlot) ──

epi_markers <- list(
  Tumour = tumour_markers,
  Mature_luminal   = mature_luminal,
  Luminal_prog     = luminal_prog,
  Basal_cell       = basal_cell,
  Epididymal       = epididymal,
  Seminal_vesicle  = seminal_vesicle
)

# ── 1b. Epithelial module scores + binarised positivity (before clearing) ─────
# AddModuleScore each prostate / Wolffian epithelial signature, binarise every
# score with binarise_expression() (GMM background threshold on the per-cell
# score), then overlay the per-module pos/neg calls on the 9.3 pearsonbatchumap
# in one combined figure. Negative module scores are non-positive by definition
# (the GMM only models score > 0), so they fall to "neg" cleanly.
DefaultAssay(srt) <- "SCT"
srt <- AddModuleScore(srt, features = epi_markers, name = "epimod_",
                      assay = "SCT", seed = 42)
score_cols <- setNames(paste0("epimod_", seq_along(epi_markers)), names(epi_markers))
bin_cols   <- setNames(paste0(names(epi_markers), "_pos"),        names(epi_markers))

for (m in names(epi_markers)) {
  sc  <- setNames(srt@meta.data[[score_cols[[m]]]], colnames(srt))
  cat("\nBinarising module score:", m, "\n")
  bin <- binarise_expression(sc, verbose = TRUE)
  srt@meta.data[[bin_cols[[m]]]] <- factor(ifelse(bin == 1L, "pos", "neg"),
                                           levels = c("neg", "pos"))
}

pos_pal   <- c(neg = "lightgrey", pos = "#e41a1c")
bin_plots <- lapply(names(epi_markers), function(m) {
  DimPlot(srt, reduction = "pearsonbatchumap", group.by = bin_cols[[m]],
          cols = pos_pal, pt.size = 0.5, order = "pos") +
    ggtitle(sprintf("%s+ (%d cells)", m, sum(srt@meta.data[[bin_cols[[m]]]] == "pos"))) +
    theme(legend.position = "none", plot.title = element_text(size = 9))
})
ggsave(file.path(png_dir, "0_epi_module_binarised.png"),
       wrap_plots(bin_plots, ncol = 3), width = 15, height = 10, dpi = 400)
cat("\nWrote epithelial module binarisation overlay -> 0_epi_module_binarised.png\n")

# ── 2. Drop the discarded clusters + Tumour-module-positive cells ─────────────
# Beyond the curated "Removed" clusters, also clear any cell called positive for
# the Tumour module score (likely tumour contamination) so it never enters the
# final normal set / CSV / downstream per-sample propagation.
n_tum_pos <- sum(srt$Tumour_pos == "pos", na.rm = TRUE)
keep <- colnames(srt)[srt$final_annotation != "Removed" & srt$Tumour_pos != "pos"]
srt  <- subset(srt, cells = keep)
cat("\nDropped", n_tum_pos, "Tumour-module-positive cells\n")
cat("Kept", ncol(srt), "cells after clearing Removed + Tumour+\n")
print(table(srt$final_annotation, useNA = "ifany"))

# ── 3. SCTransform -> HVGs on the SCT assay (mirror 9.3) ──────────────────────
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

# ── 6. Visualise new clustering + annotations on pearsonbatchumap (mirror 9.3) ─
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

p_final  <- plot_anno("final_annotation",         "Final normal annotation")
p_scim_s <- plot_anno("scim_broad",               "SCimilarity single cellsearch (broad)")
p_scim_c <- plot_anno("scim_cluster_celltype",    "SCimilarity cluster cellsearch")
p_tme    <- plot_anno("celltype_annotation_batch","TME annotation pipeline (batch)")

ggsave(file.path(png_dir, "2_final_annotation.png"),         p_final,  width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "3_scimilarity_single_broad.png"), p_scim_s, width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "4_scimilarity_cluster.png"),      p_scim_c, width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "5_TME_annotation.png"),           p_tme,    width = 9, height = 7, dpi = 300)
ggsave(file.path(png_dir, "6_annotation_overlay_grid.png"),
       (p_clu | p_final) / (p_scim_s | p_scim_c | p_tme), width = 21, height = 14, dpi = 300)

# Per-slide composition of the final annotation.
fpal <- (function(v) {
  others <- sort(unique(v))
  setNames(as.vector(polychrome())[seq_along(others)], others)
})(srt$final_annotation)
comp <- as.data.frame(table(slide = srt$slide, final_annotation = srt$final_annotation))
bar_comp <- ggplot(comp, aes(x = slide, y = Freq, fill = final_annotation)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = fpal, name = "Annotation") +
  labs(title = "Final annotation composition per slide", y = "Proportion") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(png_dir, "7_composition_bar.png"), bar_comp, width = 9, height = 6, dpi = 300)

# ── 7. Meta-program overlap with cluster DEGs (mirror 9.3) ────────────────────
sheetname <- setdiff(excel_sheets(meta_xlsx), "Malignant")
# purrr:: qualified: binarise_expression() attaches mclust, whose map() would
# otherwise shadow purrr::map() on the search path here.
meta_programs <- set_names(sheetname, sheetname) |>
  purrr::map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) |>
  purrr::map(~ purrr::map(.x, ~ as.character(na.omit(.x))))

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

layer_summary <- overlap_df %>%
  group_by(cluster, first_layer) %>%
  summarise(n_overlap = sum(n_overlap), .groups = "drop")
write.csv(layer_summary, file.path(out_dir, "cluster_firstlayer_overlap_summary.csv"),
          row.names = FALSE)

# ── 8. Faceted bar plot: facet by first layer, x = program, fill = cluster ────
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
ggsave(file.path(png_dir, "8_metaprogram_overlap_bar.png"),
       bar, width = 22, height = 14, dpi = 200, limitsize = FALSE)

# ── 9. Export per-cell final annotation (join key for 9.4.2) + object ─────────
# barcode strips the "{slide}_" prefix that merge(add.cell.ids) added in 8.1.0.
anno_df <- data.frame(
  cell_ID          = srt$cell_ID,
  slide            = srt$slide,
  barcode          = str_remove(srt$cell_ID, paste0("^", srt$slide, "_")),
  source_cluster   = srt$source_cluster,
  final_annotation = srt$final_annotation,
  stringsAsFactors = FALSE
)
write.csv(anno_df, anno_csv, row.names = FALSE)
cat("Wrote", nrow(anno_df), "rows ->", anno_csv, "\n")

qs_save(srt, file.path(out_dir, "normal_srt_final_anno.qs2"))
cat("\nDone (integrative stage). Outputs in", out_dir, "\n")
