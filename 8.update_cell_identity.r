#!/usr/bin/env Rscript
# 8.update_cell_identity.r
# Manual cell-identity correction, applied AFTER 8.5.2.general_layer_analysis.R
# has produced its per-group recluster objects. Lets you relabel a cluster's
# `final_annotation` (e.g. after reviewing a general_layer group's DE/DotPlot
# output and deciding a subcluster is mislabeled) by editing `final_annotation`
# directly on the two objects everything downstream is keyed off of, instead
# of re-running the heavy 8.4 integration script that originally produced them:
#   - 8.4.final_clear_normal_integration/normal_srt_final_anno.qs2
#       (read directly, by 8.5.2, so it must be patched in-place)
#   - 8.4.final_clear_normal_integration/final_annotation.csv
#       (re-read + rejoined by barcode in every downstream per-sample/
#       integration script — see 8.6.final_clear_normal_persample.R:54-58 for
#       the canonical join pattern)
#
# `general_layer` is a pure lookup from `final_annotation` (see
# 8.5.additional_analysis.R's `general_layer_map`), so it's recomputed here
# rather than corrected independently. If a correction's new_identity isn't
# already in that lookup, you must supply `new_general_layer` explicitly.
#
# HOW TO USE: edit the `corrections` list below (one entry per cluster you
# want to relabel — mirrors the manual-mapping style of 4.1refine_subclones.R),
# then:
#   Rscript 8.update_cell_identity.r
#
# This script does NOT re-run any downstream analysis. It patches the two
# source-of-truth objects above and prints which scripts consume
# final_annotation/general_layer, so you know what to re-run to regenerate
# their visualizations/tables against the corrected identities.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(qs2)
})

# ══════════════════════════════════════════════════════════════════════════════
# EDIT THIS: one entry per cluster to relabel.
#   group        - general_layer group name (matches
#                  8.5.normal_cell_subtypes/general_layer/<group>_recluster_srt.qs2)
#   cluster_col  - "rpca_clusters" or "pearson_clusters_batch"
#   cluster_ids  - cluster value(s) (as they appear in that column) to relabel
#   new_identity - the corrected final_annotation label
#   new_general_layer - only needed if new_identity isn't one of the existing
#                  final_annotation categories already mapped in general_layer_map
# ══════════════════════════════════════════════════════════════════════════════
corrections <- list(
  # list(group = "Myeloid", cluster_col = "rpca_clusters", cluster_ids = c("3"),
  #      new_identity = "Dendritic cells", new_general_layer = "Myeloid")
)

if (!length(corrections))
  stop("`corrections` is empty — edit the list at the top of this script before running.")

# Same fixed lookup as 8.5.additional_analysis.R (kept in sync manually).
general_layer_map <- c(
  "B/T cells"     = "Lymphoid",
  "Plasma"        = "Lymphoid",
  "Macrophages"   = "Myeloid",
  "Epithelial"    = "Epithelial",
  "SVEC"          = "Epithelial",
  "Smooth muscle" = "Stromal",
  "CAF"           = "Stromal",
  "Endo/Pericyte" = "Stromal",
  "Glial cells"   = "Neural"
)

general_layer_dir <- path.expand("~/VisHD/8.5.normal_cell_subtypes/general_layer")
final_anno_dir    <- path.expand("~/VisHD/8.4.final_clear_normal_integration")
srt_path          <- file.path(final_anno_dir, "normal_srt_final_anno.qs2")
csv_path          <- file.path(final_anno_dir, "final_annotation.csv")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Resolve each correction to a set of cell_IDs + its new final_annotation
# ══════════════════════════════════════════════════════════════════════════════
resolved <- lapply(corrections, function(cor) {
  group_path <- file.path(general_layer_dir, paste0(cor$group, "_recluster_srt.qs2"))
  if (!file.exists(group_path))
    stop("No cached recluster object for group '", cor$group, "' at ", group_path,
         " — run 8.5.2.general_layer_analysis.R for that group first.")
  group_srt <- qs_read(group_path)
  if (!cor$cluster_col %in% colnames(group_srt@meta.data))
    stop("Column '", cor$cluster_col, "' not found on ", group_path)

  hit <- as.character(group_srt@meta.data[[cor$cluster_col]]) %in% as.character(cor$cluster_ids)
  if (!any(hit))
    stop("No cells matched cluster_ids ", paste(cor$cluster_ids, collapse = ","),
         " in column '", cor$cluster_col, "' of group '", cor$group, "'")

  new_gl <- cor$new_general_layer
  if (is.null(new_gl)) {
    new_gl <- unname(general_layer_map[cor$new_identity])
    if (is.na(new_gl))
      stop("new_identity '", cor$new_identity, "' is not in general_layer_map — ",
           "add `new_general_layer = \"...\"` to this correction entry.")
  }

  cat(sprintf("[%s / %s == %s] %d cells -> final_annotation = '%s' (general_layer = '%s')\n",
              cor$group, cor$cluster_col, paste(cor$cluster_ids, collapse = ","),
              sum(hit), cor$new_identity, new_gl))

  list(cell_ID = colnames(group_srt)[hit], new_identity = cor$new_identity, new_general_layer = new_gl)
})

# ══════════════════════════════════════════════════════════════════════════════
# 2. Patch normal_srt_final_anno.qs2 (read directly by 8.5.2 — must be updated
#    in-place; the CSV alone would not be picked up by that script)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nLoading", srt_path, "\n")
srt <- qs_read(srt_path)

backup_path <- sub("\\.qs2$", sprintf("_backup_%s.qs2", format(Sys.time(), "%Y%m%d_%H%M%S")), srt_path)
cat("Backing up current object to", backup_path, "\n")
qs_save(srt, backup_path)

n_updated <- 0
for (r in resolved) {
  idx <- match(r$cell_ID, colnames(srt))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) {
    warning("None of the resolved cell_IDs for identity '", r$new_identity, "' were found in ", srt_path)
    next
  }
  srt$final_annotation[idx] <- r$new_identity
  srt$general_layer[idx]    <- r$new_general_layer
  n_updated <- n_updated + length(idx)
}
cat("Updated final_annotation/general_layer for", n_updated, "cells\n")

qs_save(srt, srt_path)
cat("Saved", srt_path, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Patch final_annotation.csv (source of truth every downstream per-sample /
#    integration script re-joins by barcode — see 8.6.final_clear_normal_persample.R)
# ══════════════════════════════════════════════════════════════════════════════
if (file.exists(csv_path)) {
  anno <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  backup_csv <- sub("\\.csv$", sprintf("_backup_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")), csv_path)
  cat("\nBacking up current CSV to", backup_csv, "\n")
  write.csv(anno, backup_csv, row.names = FALSE)

  for (r in resolved) {
    idx <- match(r$cell_ID, anno$cell_ID)
    idx <- idx[!is.na(idx)]
    if (length(idx)) anno$final_annotation[idx] <- r$new_identity
  }
  write.csv(anno, csv_path, row.names = FALSE)
  cat("Saved", csv_path, "\n")
} else {
  cat("\nNo final_annotation.csv found at", csv_path, "— skipping CSV update\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Point out downstream steps that consume final_annotation/general_layer and
#    so need to be re-run to reflect this correction in their outputs
# ══════════════════════════════════════════════════════════════════════════════
cat("\n== Scripts referencing final_annotation / general_layer (re-run to refresh) ==\n")
script_files <- list.files(path.expand("~/VisHD"), pattern = "\\.(R|r|Rmd|ipynb)$", full.names = TRUE)
script_files <- script_files[basename(script_files) != "8.update_cell_identity.r"]
hits <- Filter(function(f) {
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  any(grepl("final_annotation|general_layer", txt))
}, script_files)
if (length(hits)) {
  cat(paste0("  - ", sort(basename(hits))), sep = "\n")
} else {
  cat("  (none found)\n")
}
cat("\nNote: 8.5.2's own general_layer/<group>_recluster_srt.qs2 caches were NOT\n",
    "regenerated by this script (they still reflect pre-correction clustering) —\n",
    "only their source final_annotation/general_layer labels were patched.\n",
    "Re-run 8.5.2 for an affected group if you need its DE/plots to reflect the\n",
    "new cluster composition, not just the corrected label.\n", sep = "")
