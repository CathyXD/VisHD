library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")
library(infercnv)
source("~/VisHD/functions.R")
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

# ── Sample selection (Slurm array index 1-8) ─────────────────────────────────
arg    <- as.integer(commandArgs(trailingOnly = TRUE)[1])
paths  <- system("realpath ~/VisHD/LUT-245-*/raw_srt.qs", intern = TRUE)
slides <- basename(dirname(paths))
sample <- slides[arg]
cat("Processing sample", arg, ":", sample, "\n")

out_dir <- path.expand(file.path("~/VisHD", sample, "second_round_check"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1-6. Annotate + SpaNorm + Pearson PCA (resume from tumour_anno_srt2.qs2) ──
# Resume point: if the annotated/PCA'd object already exists, load it and skip
# the integrated-id lookup, raw load, SpaNorm, Pearson PCA and UMAP plots.
anno_path <- path.expand(file.path("~/VisHD", sample, "tumour_anno_srt2.qs2"))

if (file.exists(anno_path)) {
  cat("Resuming: loading", anno_path, "(skip annotation + SpaNorm + Pearson PCA)\n")
  srt <- qs_read(anno_path)
} else {
  # ── 1. Tumour/Normal cell ids for this sample from the integrated objects ──
  # Integrated colnames are "<slide>_<barcode>"; strip the slide prefix to
  # recover the plain integer barcode used in the per-sample raw object.
  int_dir <- path.expand("~/VisHD/1.1.integrate_raw_cell")

  tumour_int <- qs_read(file.path(int_dir, "tumour_srt.qs2"))
  tumour_ids <- as.integer(gsub(paste0("^", sample, "_"), "",
                                colnames(tumour_int)[tumour_int$slide == sample]))
  rm(tumour_int); gc()

  normal_int <- qs_read(file.path(int_dir, "normal_cleaned_srt.qs2"))
  normal_ids <- as.integer(gsub(paste0("^", sample, "_"), "",
                                colnames(normal_int)[normal_int$slide == sample]))
  rm(normal_int); gc()

  cat("Integrated cells for", sample, "- Tumour:", length(tumour_ids),
      "Normal:", length(normal_ids), "\n")

  # ── 2. Load raw object, assign tumour_anno, drop cells in neither set ──────
  srt <- qread(path.expand(file.path("~/VisHD", sample, "raw_srt.qs")))  # qs::qread (.qs)
  srt <- UpdateSeuratObject(srt)   # refresh FOV class (older objects lack the coords_x_orientation slot -> subset() fails)
  ids <- as.integer(colnames(srt))
  srt$tumour_anno <- ifelse(ids %in% tumour_ids, "Tumour",
                     ifelse(ids %in% normal_ids, "Normal", NA_character_))

  keep <- colnames(srt)[!is.na(srt$tumour_anno)]
  cat("Cells kept:", length(keep), "/ dropped (in neither set):",
      ncol(srt) - length(keep), "\n")
  srt <- subset(srt, cells = keep)
  cat("tumour_anno:\n"); print(table(srt$tumour_anno))

  # ── 3. category from this sample's category.csv (matched by integer cellid) ─
  CB <- read.csv(path.expand(file.path("~/VisHD", sample, "category.csv")))
  CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
  srt$category <- CB$category[match(as.integer(colnames(srt)), CB$cellid)]
  cat("category:\n"); print(table(srt$category, useNA = "ifany"))

  # ── 4. SpaNorm, then Pearson PCA on the SpaNorm SVGs ───────────────────────
  srt <- do.spanorm(srt, outdir = out_dir)                # builds SpaNorm assay + SVGs.Rds
  VariableFeatures(srt, assay = "Spatial") <- VariableFeatures(srt, assay = "SpaNorm")
  srt <- do.pearson_pca(srt, find_hvgs = FALSE)           # pearson PCA on the SVGs (no batch)

  # ── 5. Visualisation: pearsonumap by category and tumour_anno ──────────────
  p1 <- DimPlot(srt, reduction = "pearsonumap", group.by = "category",
                label = TRUE, label.size = 3) +
    ggtitle(paste(sample, "- category (CB / DT)")) +
    theme(legend.key.size = unit(0.3, "cm"))
  p2 <- DimPlot(srt, reduction = "pearsonumap", group.by = "tumour_anno",
                label = TRUE) +
    scale_color_manual(values = c(Normal = "lightgrey", Tumour = "red")) +
    ggtitle(paste(sample, "- tumour vs normal"))
  ggsave(file.path(out_dir, "1_pearsonumap_category.png"),    p1, width = 7, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "2_pearsonumap_tumour_anno.png"), p2, width = 6, height = 5, dpi = 300)

  # ── 6. Save annotated object ─────────────────────────────────────────────
  qs_save(srt, anno_path)
}

# ── Spatial map of tumour_anno (ImageDimPlot); always (re)generated ──────────
p3 <- ImageDimPlot(srt, group.by = "tumour_anno", size = 1.2) +
  scale_fill_manual(values = c(Normal = "lightgrey", Tumour = "red")) +
  ggtitle(paste(sample, "- tumour vs normal (spatial)"))
ggsave(file.path(out_dir, "3_spatial_tumour_anno.png"), p3, width = 7, height = 6, dpi = 300)

# ── Pearson clusters: UMAP (DimPlot) + spatial (ImageDimPlot) ────────────────
p4 <- DimPlot(srt, reduction = "pearsonumap", group.by = "pearson_clusters",
              label = TRUE, label.size = 3) +
  ggtitle(paste(sample, "- pearson clusters (UMAP)"))
ggsave(file.path(out_dir, "4_pearsonumap_clusters.png"), p4, width = 7, height = 5, dpi = 300)

p5 <- ImageDimPlot(srt, group.by = "pearson_clusters", size = 1.2) +
  ggtitle(paste(sample, "- pearson clusters (spatial)"))
ggsave(file.path(out_dir, "5_spatial_pearson_clusters.png"), p5, width = 7, height = 6, dpi = 300)

# ── infercnv_group: reference vs test by per-cluster Normal proportion ───────
# Annotation is pearson_clusters. A cluster is a CNV *reference* when >90% of
# its cells are Normal (tumour_anno); every other cluster is a *test* group.
normal_prop  <- tapply(srt$tumour_anno == "Normal", srt$pearson_clusters,
                       mean, na.rm = TRUE)
cat("Normal proportion per pearson cluster:\n"); print(round(normal_prop, 3))
ref_clusters <- names(normal_prop)[normal_prop > 0.9]
if (length(ref_clusters) == 0) {
  # No cluster clears 0.9 — fall back to the single most-Normal cluster.
  ref_clusters <- names(normal_prop)[which.max(normal_prop)]
  cat(sprintf("No cluster >90%% Normal — fallback to most-Normal cluster %s (%.1f%% Normal)\n",
              ref_clusters, 100 * max(normal_prop, na.rm = TRUE)))
}
cat("Reference clusters:", paste(ref_clusters, collapse = ", "), "\n")

srt$infercnv_group <- ifelse(srt$pearson_clusters %in% ref_clusters,
                             "Reference", "Test")
cat("infercnv_group:\n"); print(table(srt$infercnv_group))

p6 <- DimPlot(srt, reduction = "pearsonumap", group.by = "infercnv_group",
              label = TRUE) +
  scale_color_manual(values = c(Reference = "steelblue", Test = "firebrick")) +
  ggtitle(paste(sample, "- infercnv_group (reference vs test)"))
ggsave(file.path(out_dir, "6_pearsonumap_infercnv_group.png"), p6, width = 7, height = 5, dpi = 300)

p7 <- ImageDimPlot(srt, group.by = "infercnv_group", size = 1.2) +
  scale_fill_manual(values = c(Reference = "steelblue", Test = "firebrick")) +
  ggtitle(paste(sample, "- infercnv_group (spatial)"))
ggsave(file.path(out_dir, "7_spatial_infercnv_group.png"), p7, width = 7, height = 6, dpi = 300)

# Persist the infercnv_group column back into the annotated object (full cells).
qs_save(srt, anno_path)

# ── 7. InferCNV (cells mode); resume from run.final.infercnv_obj if present ──
infercnv_out <- file.path(out_dir, "infercnv")
dir.create(infercnv_out, showWarnings = FALSE, recursive = TRUE)
final_rds <- file.path(infercnv_out, "run.final.infercnv_obj")

if (file.exists(final_rds)) {
  cat("Resuming: loading", final_rds, "(skip CNV inference, plot only)\n")
  infercnvobject <- readRDS(final_rds)
} else {
  # Stratified subsample to 60 k cells by pearson cluster so each reference
  # cluster keeps adequate counts (rare clusters get a 50-cell floor).
  set.seed(42)
  target_n <- 60000
  if (ncol(srt) > target_n) {
    cells_by_group <- split(colnames(srt), srt$pearson_clusters)
    group_sizes    <- lengths(cells_by_group)
    alloc <- round(group_sizes / sum(group_sizes) * target_n)
    alloc <- pmin(alloc, group_sizes)
    small <- alloc == 0 & group_sizes > 0
    alloc[small] <- pmin(50, group_sizes[small])
    keep_cells <- unlist(Map(base::sample, cells_by_group, alloc))  # `sample` var (line 15) shadows base::sample()
    srt <- subset(srt, cells = keep_cells)
    cat("Subsampled from", sum(group_sizes), "to", ncol(srt),
        "cells, stratified by pearson_clusters\n")
  }

  # Annotation = pearson_clusters; references = clusters with >90% Normal.
  annotations <- data.frame(
    cluster   = as.character(srt$pearson_clusters),
    row.names = colnames(srt)
  )
  cat("Annotation groups (count): ",
      paste(names(table(annotations$cluster)), table(annotations$cluster),
            sep = "=", collapse = ", "), "\n")

  ref_present <- intersect(ref_clusters, unique(annotations$cluster))
  cat("Reference groups:", paste(ref_present, collapse = ", "), "\n")
  if (length(ref_present) == 0) stop("No reference clusters present after subsampling")

  counts <- as.matrix(GetAssayData(srt, assay = "Spatial", layer = "counts"))

  infercnvobject <- CreateInfercnvObject(
    raw_counts_matrix       = counts,
    annotations_file        = annotations,
    delim                   = "\t",
    min_max_counts_per_cell = c(50, Inf),
    chr_exclude             = c("chrY", "chrM"),
    gene_order_file         = readRDS("~/VisHD/gene_ord2.Rds"),
    ref_group_names         = ref_present
  )

  infercnvobject <- infercnv::run(
    infercnvobject,
    cutoff             = 0.01,
    out_dir            = infercnv_out,
    analysis_mode      = "cells",
    cluster_by_groups  = TRUE,
    denoise            = FALSE,
    HMM                = FALSE,
    save_rds           = TRUE,
    plot_steps         = FALSE,
    write_phylo        = TRUE,
    write_expr_matrix  = FALSE,
    num_threads        = 8,
    resume_mode        = TRUE
  )
}

plot_per_group(
  infercnvobject,
  on_references     = TRUE,
  on_observations   = TRUE,
  base_filename     = "infercnv_per_group",
  output_format     = "png",
  write_expr_matrix = FALSE,
  png_res           = 300,
  dynamic_resize    = 0,
  useRaster         = TRUE,
  out_dir           = infercnv_out
)
