#!/usr/bin/env Rscript
# 10.1.Statial.R   (run-once, all 8 samples)
# Statial spatial analysis on the per-sample tumour+normal objects
# (LUT-245-XX/tumour_normal_anno_srt.qs2), converted to a SingleCellExperiment
# list keyed on the `final_annotation` cell identity.
#
#   * Kontextual  — context-aware CO-LOCALIZATION between tumour groups and
#                   normal cell types. Tumour groups are the G1/G2/G3(+combos)
#                   / Neg labels; the "normal" parent population is built
#                   automatically from every non-tumour annotation.
#   * SpatioMark  — continuous expression (SpaNorm logcounts) of the tumour
#                   group-DEG genes (groupdeg.rds: G1/G2/G3) *in tumour cells*
#                   as a function of proximity to each normal cell type
#                   (calcStateChanges over the distance to each normal type).
#
# Survival analysis is intentionally SKIPPED.
#
# Reference: https://bioconductor.org/packages/release/bioc/vignettes/Statial/
#            inst/doc/Statial.html
#
# Outputs under ~/VisHD/10.1.Statial/:
#   kontextual/  Kontextual result table + tumour<->normal co-localization heatmaps
#   spatiomark/  calcStateChanges result table + state-change heatmaps / plots
#
#   Rscript 10.1.Statial.R
#
# NOTE: Statial (+ its Bioconductor deps) must be installed into ~/R_Library/4.5,
#       e.g.  BiocManager::install("Statial", lib = "~/R_Library/4.5")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(qs2)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(Statial, lib.loc = "~/R_Library/4.5")
})

setwd("~/VisHD")

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(nCores) || nCores < 1) nCores <- 1
cat("Using", nCores, "core(s)\n")

# Radii (in FOV pixel units; ~0.29 µm/px, so ≈ 15 / 30 / 60 µm). Tune as needed.
radii   <- c(50, 100, 200)
maxDist <- 200          # SpatioMark: max distance considered by getDistances()

# The Normal compartment is subdivided: Immune = these cell types, Stromal = the rest.
immune_cell_types <- c("Macrophages", "B cells", "Plasma")

out_dir  <- "~/VisHD/10.1.Statial"
kon_dir  <- file.path(out_dir, "kontextual")
spm_dir  <- file.path(out_dir, "spatiomark")
for (d in c(kon_dir, spm_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Reference: tumour group-DEG gene sets (G1/G2/G3) ──────────────────────────
groupdeg <- readRDS(paste0("~/VisHD/6.2archetype_downstream_tumour/archetype_module/",
                           "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))
gene2group <- stack(lapply(groupdeg, as.character))          # values, ind(=G1/G2/G3)
gene2group <- setNames(as.character(gene2group$ind), gene2group$values)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Build a SingleCellExperiment list (one per sample) → combined multi-image SCE
# ══════════════════════════════════════════════════════════════════════════════
paths   <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
samples <- basename(paths)

build_sce <- function(path, sample) {
  srt <- qs_read(file.path(path, "tumour_normal_anno_srt.qs2"))
  if (length(SeuratObject::Layers(srt, assay = "Spatial")) > 1)
    srt <- JoinLayers(srt, assay = "Spatial")
  cell_id <- paste(sample, colnames(srt), sep = "_")
  logmat   <- GetAssayData(srt, assay = "SpaNorm", layer = "data")     # continuous expr
  colnames(logmat) <- cell_id
  countmat <- GetAssayData(srt, assay = "Spatial", layer = "counts")
  colnames(countmat) <- cell_id
  # Coords live in the FOV; rows are ordered 1:n with the barcode in coords$cell.
  coords <- GetTissueCoordinates(srt, which = "centroids")
  idx    <- match(colnames(srt), coords$cell)
  stopifnot(!anyNA(idx))

  cellty <- as.character(srt$final_annotation)
  comp   <- as.character(srt$compartment)
  # Subdivide the Normal compartment into Immune vs Stromal.
  is_normal <- comp == "Normal"
  comp[is_normal] <- ifelse(cellty[is_normal] %in% immune_cell_types,
                            "Immune", "Stromal")

  cd <- S4Vectors::DataFrame(
    cellType    = cellty,
    compartment = comp,
    imageID     = sample,
    x           = coords$x[idx],
    y           = coords$y[idx],
    row.names   = cell_id
  )
  sce <- SingleCellExperiment(
    assays  = list(logcounts = logmat, counts = countmat),
    colData = cd
  )
  cat("  ", sample, ":", ncol(sce), "cells,", nrow(sce), "genes\n")
  sce
}

cat("Building SCE list...\n")
sce_list <- Map(build_sce, paths, samples)

# Combine on common genes into one multi-image SCE (Statial iterates over imageID).
common   <- Reduce(intersect, lapply(sce_list, rownames))
sce_list <- lapply(sce_list, function(s) s[common, ])
sce      <- do.call(SingleCellExperiment::cbind, sce_list)
qs_save(sce, file.path(out_dir, "combined_sce.qs2"))
rm(sce_list); gc()
cat("Combined SCE:", ncol(sce), "cells across", length(samples), "images,",
    nrow(sce), "genes\n")

# ── Auto-detect tumour groups vs the normal parent population ──────────────────
comp_map      <- unique(as.data.frame(colData(sce)[, c("cellType", "compartment")]))
tumour_types  <- sort(comp_map$cellType[comp_map$compartment == "Tumour"])
immune_types  <- sort(comp_map$cellType[comp_map$compartment == "Immune"])
stromal_types <- sort(comp_map$cellType[comp_map$compartment == "Stromal"])
non_tumour_types <- c(immune_types, stromal_types)          # all Normal (Immune + Stromal)
all_types     <- sort(unique(as.character(sce$cellType)))
cat("\nTumour groups (", length(tumour_types), "): ",
    paste(tumour_types, collapse = ", "), "\n", sep = "")
cat("Immune types  (", length(immune_types), "): ",
    paste(immune_types, collapse = ", "), "\n", sep = "")
cat("Stromal types (", length(stromal_types), "): ",
    paste(stromal_types, collapse = ", "), "\n\n", sep = "")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Kontextual — tumour <-> normal co-localization
# ══════════════════════════════════════════════════════════════════════════════
cat("Running Kontextual...\n")

# parentCombinations builds (from, to, parent) triples: `from` drawn from each
# parent population, `to` drawn from `all`. Three parents -> tumour, immune, stromal.
parentDf <- parentCombinations(
  all     = all_types,
  tumour  = tumour_types,
  immune  = immune_types,
  stromal = stromal_types
)

# Focus on tumour <-> (immune|stromal) only; drop within-compartment pairs
# (tumour-tumour, immune/stromal-immune/stromal).
cross <- (parentDf$from %in% tumour_types     & parentDf$to %in% non_tumour_types) |
         (parentDf$from %in% non_tumour_types & parentDf$to %in% tumour_types)
parentDf <- parentDf[cross, , drop = FALSE]
cat("  ", nrow(parentDf), "tumour<->normal from/to pairs\n")

kontext <- do.call(rbind, lapply(radii, function(rr) {
  cat("  r =", rr, "\n")
  k <- Kontextual(
    cells         = sce,
    parentDf      = parentDf,
    r             = rr,
    cellType      = "cellType",
    imageID       = "imageID",
    spatialCoords = c("x", "y"),
    cores         = nCores
  )
  k$r <- rr
  k
}))
saveRDS(kontext, file.path(kon_dir, "kontextual_results.rds"))
write.csv(kontext, file.path(kon_dir, "kontextual_results.csv"), row.names = FALSE)
cat("Kontextual done. columns:", paste(colnames(kontext), collapse = ", "), "\n")

# ── Kontextual heatmaps (wrapped: never lose the saved table on a plotting error)
tryCatch({
  kt <- kontext
  # Resolve from/to columns (Statial encodes the pair in `test` as "from__to").
  if (!all(c("from", "to") %in% colnames(kt)) && "test" %in% colnames(kt))
    kt <- tidyr::separate(kt, test, into = c("from", "to"), sep = "__", remove = FALSE)
  val_col <- if ("kontextual" %in% colnames(kt)) "kontextual" else
             tail(names(Filter(is.numeric, kt)), 1)
  img_col <- if ("imageID" %in% colnames(kt)) "imageID" else "image"
  kt$value <- kt[[val_col]]
  # Orient so rows = tumour group, cols = normal type (immune types then stromal).
  kt$tumour <- ifelse(kt$from %in% tumour_types, kt$from, kt$to)
  kt$normal <- ifelse(kt$from %in% non_tumour_types, kt$from, kt$to)
  kt$normal <- factor(kt$normal, levels = c(immune_types, stromal_types))

  lim <- max(abs(kt$value), na.rm = TRUE)
  for (rr in radii) {
    d <- kt[kt$r == rr, , drop = FALSE]
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(normal, tumour, fill = value)) +
      geom_tile(colour = "grey90") +
      facet_wrap(as.formula(paste0("~", img_col))) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                           midpoint = 0, limits = c(-lim, lim), name = "Kontextual") +
      labs(title = paste0("Tumour<->Normal co-localization (r=", rr, ")"),
           x = "Normal cell type", y = "Tumour group") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(kon_dir, paste0("kontextual_heatmap_r", rr, "_by_image.png")),
           p, width = 16, height = 12, dpi = 150, limitsize = FALSE)
  }
  # Cross-image mean summary (one panel per radius).
  s <- kt %>% group_by(tumour, normal, r) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  p <- ggplot(s, aes(normal, tumour, fill = value)) +
    geom_tile(colour = "grey90") +
    facet_wrap(~ r, labeller = label_both) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                         midpoint = 0, name = "mean\nKontextual") +
    labs(title = "Tumour<->Normal co-localization (mean across samples)",
         x = "Normal cell type", y = "Tumour group") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(kon_dir, "kontextual_heatmap_mean.png"),
         p, width = 14, height = 6, dpi = 150, limitsize = FALSE)
}, error = function(e) cat("Kontextual plotting skipped:", conditionMessage(e), "\n"))

# ══════════════════════════════════════════════════════════════════════════════
# 3. SpatioMark — continuous state (gene) changes in TUMOUR cells vs proximity
#                 to each NORMAL cell type, for the group-DEG genes.
# ══════════════════════════════════════════════════════════════════════════════
cat("\nRunning SpatioMark...\n")

features <- intersect(unique(unlist(groupdeg, use.names = FALSE)), rownames(sce))
cat("  group-DEG genes present:", length(features), "/",
    length(unique(unlist(groupdeg))), "\n")

# Distance from every cell to the nearest cell of each type -> reducedDim.
sce <- getDistances(sce, maxDist = maxDist, nCores = nCores)

stateChanges <- calcStateChanges(
  cells    = sce,
  type     = "distances",
  from     = tumour_types,     # expression measured in tumour cells
  to       = non_tumour_types, # proximity to each immune/stromal cell type
  marker   = features,
  assay    = "logcounts",
  cellType = "cellType",
  imageID  = "imageID",
  minCells = 20,
  nCores   = nCores
)

# Annotate each marker with its G1/G2/G3 group of origin.
mk_col <- intersect(c("marker", "feature"), colnames(stateChanges))[1]
if (!is.na(mk_col)) stateChanges$group <- gene2group[stateChanges[[mk_col]]]

saveRDS(stateChanges, file.path(spm_dir, "state_changes.rds"))
write.csv(stateChanges, file.path(spm_dir, "state_changes.csv"), row.names = FALSE)
cat("calcStateChanges done. columns:", paste(colnames(stateChanges), collapse = ", "), "\n")

# ── SpatioMark heatmaps: coefficient (marker × normal type) per tumour group ──
tryCatch({
  sc <- as.data.frame(stateChanges)
  from_col <- intersect(c("primaryCellType", "from", "imageCellType"), colnames(sc))[1]
  to_col   <- intersect(c("otherCellType", "to"),                       colnames(sc))[1]
  coef_col <- intersect(c("coef", "tval", "estimate"),                  colnames(sc))[1]
  fdr_col  <- intersect(c("fdr", "adjPval", "pval"),                    colnames(sc))[1]
  stopifnot(!is.na(from_col), !is.na(to_col), !is.na(coef_col), !is.na(mk_col))

  sc$.from <- sc[[from_col]]; sc$.to <- sc[[to_col]]
  sc$.coef <- sc[[coef_col]]; sc$.mk <- sc[[mk_col]]
  # Mean coefficient across samples/images per (tumour group, normal type, gene).
  agg <- sc %>% group_by(.from, .to, .mk, group) %>%
    summarise(coef = mean(.coef, na.rm = TRUE), .groups = "drop")

  lim <- max(abs(agg$coef), na.rm = TRUE)
  for (ft in sort(unique(agg$.from))) {
    d <- agg[agg$.from == ft, , drop = FALSE]
    if (nrow(d) == 0) next
    d$.mk <- factor(d$.mk, levels = names(sort(gene2group[unique(d$.mk)])))
    p <- ggplot(d, aes(.to, .mk, fill = coef)) +
      geom_tile(colour = "grey92") +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                           midpoint = 0, limits = c(-lim, lim), name = "state\nchange") +
      labs(title = paste0("SpatioMark: expression of group-DEG genes in '", ft,
                          "' vs proximity to normal cells"),
           x = "Normal cell type (proximity)", y = "Group-DEG gene") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    nh <- max(6, length(unique(d$.mk)) * 0.14)
    ggsave(file.path(spm_dir, paste0("statechange_heatmap_", make.names(ft), ".png")),
           p, width = 10, height = nh, dpi = 150, limitsize = FALSE)
  }

  # plotStateChanges for the strongest tumour<->normal gene relationships.
  if (!is.na(fdr_col)) {
    top <- sc[order(sc[[fdr_col]]), , drop = FALSE]
    top <- head(top[!is.na(top$.coef), ], 12)
    img_col <- intersect(c("imageID", "image"), colnames(sc))[1]
    for (k in seq_len(nrow(top))) {
      r <- top[k, ]
      pf <- tryCatch(
        plotStateChanges(cells = sce, type = "distances",
                         image  = as.character(r[[img_col]]),
                         from   = as.character(r$.from),
                         to     = as.character(r$.to),
                         marker = as.character(r$.mk),
                         size = 1, shape = 19, interactive = FALSE,
                         plotModelFit = TRUE, method = "lm"),
        error = function(e) NULL)
      if (is.null(pf)) next
      g <- if (inherits(pf, "list") && !is.null(pf$image)) pf$image else pf
      ggsave(file.path(spm_dir, sprintf("top%02d_%s_%s_%s_%s.png", k,
                       make.names(as.character(r[[img_col]])),
                       make.names(as.character(r$.from)),
                       make.names(as.character(r$.to)),
                       make.names(as.character(r$.mk)))),
             g, width = 7, height = 5, dpi = 150)
    }
  }
}, error = function(e) cat("SpatioMark plotting skipped:", conditionMessage(e), "\n"))

cat("\n==================== 10.1.Statial done ====================\n")
