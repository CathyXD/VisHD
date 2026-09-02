#!/usr/bin/env Rscript
# 10.3.SPIAT_tumour_border.R   (per-sample, array 1-8)
# SPIAT tissue-structure workflow on the per-sample tumour+normal object
# (LUT-245-XX/tumour_normal_anno_srt.qs2): detects the tumour margin
# (Internal margin / External margin / Inside / Outside) using `cell_type`
# as the reference-cell feature column (Tumour subclones 1/2 pooled into a
# single "Tumour" reference level), then records each cell's structure
# position alongside its `final_annotation` identity.
#
# Reference: https://trigosteam.github.io/SPIAT/articles/tissue-structure.html
#
# Outputs under ~/VisHD/10.3.tumour_border/:
#   per_sample_spe/    <slide>_structure_spe.qs2   (SPIAT-formatted SPE w/ Structure col)
#   per_sample_plots/  <slide>_structure.png       (plot_cell_categories by Structure)
#   per_sample_tables/ <slide>_cell_position.csv   (per-cell: cell, x, y, slide, cell_type,
#                                                    final_annotation, category, category_bin,
#                                                    Structure)
#
#   Rscript 10.3.SPIAT_tumour_border.R <sample-index 1-8>
#
# NOTE: SPIAT must be installed into ~/R_Library/4.5,
#       e.g.  BiocManager::install("SPIAT", lib = "~/R_Library/4.5")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs2)
  library(SpatialExperiment)
  library(S4Vectors)
  library(SPIAT)
  library(alphahull, lib.loc = "~/R_Library/4.5")
})

# ── CLI arg ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 10.3.SPIAT_tumour_border.R <sample-index 1-8>")
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[arg]
i     <- basename(path)
cat("Working at", path, "\n")

out_dir  <- "~/VisHD/10.3.tumour_border"
spe_dir  <- file.path(out_dir, "per_sample_spe")
plot_dir <- file.path(out_dir, "per_sample_plots")
tbl_dir  <- file.path(out_dir, "per_sample_tables")
for (d in c(spe_dir, plot_dir, tbl_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

n_margin_layers <- 5

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load object, derive feature columns
# ══════════════════════════════════════════════════════════════════════════════
srt <- qs_read(file.path(path, "tumour_normal_anno_srt.qs2"))
if (length(SeuratObject::Layers(srt, assay = "Spatial")) > 1)
  srt <- JoinLayers(srt, assay = "Spatial")

meta <- srt@meta.data
meta$category_bin <- ifelse(meta$category == "DT", "DT",
                      ifelse(grepl("^CB", meta$category), "CB", NA))
# Pool Tumour subclones ("Tumour 1"/"Tumour 2") into a single reference level.
meta$cell_type <- ifelse(grepl("^Tumour", meta$cell_type), "Tumour", meta$cell_type)

coords <- GetTissueCoordinates(srt, which = "centroids")
idx    <- match(colnames(srt), coords$cell)
stopifnot(!anyNA(idx))
xy       <- as.matrix(coords[idx, c("x", "y")])
countmat <- GetAssayData(srt, assay = "Spatial", layer = "counts")

cat("  ", i, ":", ncol(srt), "cells,", sum(meta$cell_type == "Tumour"), "Tumour cells\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Build SPIAT-formatted SpatialExperiment
# ══════════════════════════════════════════════════════════════════════════════
formatted_image <- format_image_to_spe(
  format           = "general",
  intensity_matrix = as.matrix(countmat),
  phenotypes       = meta$cell_type,
  coord_x          = xy[, 1],
  coord_y          = xy[, 2]
)
colData(formatted_image)$cell_type        <- meta$cell_type
colData(formatted_image)$final_annotation <- meta$final_annotation
colData(formatted_image)$category         <- meta$category
colData(formatted_image)$category_bin     <- meta$category_bin
colData(formatted_image)$slide            <- i

# ══════════════════════════════════════════════════════════════════════════════
# 3. Detect tumour border / margin structure
# ══════════════════════════════════════════════════════════════════════════════
non_tumour_types <- setdiff(unique(meta$cell_type), "Tumour")

result <- tryCatch({
  b <- identify_bordering_cells(formatted_image, reference_cell = "Tumour",
                                 feature_colname = "cell_type")
  d <- calculate_distance_to_margin(b)
  define_structure(d, cell_types_of_interest = non_tumour_types,
                   feature_colname = "cell_type", n_margin_layers = n_margin_layers)
}, error = function(e) {
  message(i, ": SPIAT structure detection failed — ", conditionMessage(e))
  NULL
})

if (is.null(result)) {
  cat(i, ": no output written (structure detection failed)\n")
  quit(save = "no", status = 0)
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Save outputs
# ══════════════════════════════════════════════════════════════════════════════
qs_save(result, file.path(spe_dir, paste0(i, "_structure_spe.qs2")))

p <- tryCatch(
  plot_cell_categories(result, feature_colname = "Structure") +
    ggtitle(paste(i, "— tumour border structure")),
  error = function(e) NULL
)
if (!is.null(p))
  ggsave(file.path(plot_dir, paste0(i, "_structure.png")),
         p, width = 8, height = 7, dpi = 150, limitsize = FALSE)

cd     <- as.data.frame(colData(result))
xy_out <- as.data.frame(spatialCoords(result))
colnames(xy_out) <- c("x", "y")
cell_position <- cbind(
  cell = colnames(result),
  xy_out,
  cd[, c("slide", "cell_type", "final_annotation", "category", "category_bin", "Structure")]
)
write.csv(cell_position, file.path(tbl_dir, paste0(i, "_cell_position.csv")), row.names = FALSE)

cat(i, ": done —", nrow(cell_position), "cells;",
    "Structure levels:", paste(unique(cell_position$Structure), collapse = ", "), "\n")
