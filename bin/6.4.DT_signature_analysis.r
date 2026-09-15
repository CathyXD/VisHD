# 6.4.DT_signature_analysis.r  (per-sample, array 1-8)
# DT-only counterpart of 6.4.signature_analysis.R. Same scoring / binarisation /
# comparison workflow, but sourced from the 6.3.DT_archetype_module.r groupdeg
# (DT_G1..DT_G5 signatures derived from archetypes fit on category == "DT" cells).
# Score the DT groupdeg signatures on one sample's tumour cells, binarise each
# per-cell module score, build the combined Module_group label, run DEG +
# pathway enrichment by Module_group, and save per-cell metadata + spatial plot
# objects to disk. Skips samples already processed (meta.Rds present).
#
# Cross-sample boxplots/positivity/spatial-grid/pathway-recurrence outputs are
# produced by the companion run-once script 6.4.DT_signature_analysis_aggregate.R,
# which reads this script's per-sample outputs back in — run it once all 8 array
# tasks have completed.
#
#   Rscript 6.4.DT_signature_analysis.r <sample-index 1-8>

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(patchwork)
  library(ggpubr)
  library(RColorBrewer)
})

# ── Config ──────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

# ── CLI arg ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1]))))
  stop("Usage: Rscript 6.4.DT_signature_analysis.r <sample-index 1-8>")
s <- samples[as.numeric(args[1])]

base_dir <- "~/VisHD"
source(file.path(base_dir, "functions.R"))        # binarise_expression(), run_gsea_panel()

groupdeg <- readRDS(file.path(base_dir,
  "6.3.DT_archetype_module",
  "group_DEG_enrichment/cross_sample_summary/groupdeg.rds"))
module_names <- paste0("module_", names(groupdeg))
banksy_col   <- "BANKSY_0.2_snn_res.1"            # only BANKSY clustering in tumour_srt.qs2

# gene sets for Module_group DEG pathway enrichment (run_gsea_panel(), functions.R)
Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

outdir <- file.path(base_dir, "6.4.DT_signature_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sample_dir <- file.path(outdir, s)
meta_f     <- file.path(sample_dir, "meta.Rds")

bin_dir_s <- file.path(sample_dir, "binarisation")
deg_dir_s <- file.path(sample_dir, "module_group_deg")
png_dir_s <- file.path(sample_dir, "png")
dir.create(bin_dir_s, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir_s, showWarnings = FALSE, recursive = TRUE)
dir.create(png_dir_s, showWarnings = FALSE, recursive = TRUE)

pos_pal <- c(neg = "grey85", pos = "indianred")

# Fixed Module_group levels/palette (Neg + every signature combination) so the
# colours are consistent across samples in the aggregate spatial grid. Palette is
# generated generically (not hardcoded per group name) since the DT groupdeg has
# 5 groups.
group_combos  <- unlist(lapply(seq_along(module_names), function(k)
  combn(sub("^module_", "", module_names), k, FUN = function(x) paste(x, collapse = "/"))))
group_levels  <- c("Neg", group_combos)
single_combos <- group_combos[!grepl("/", group_combos)]        # one signature positive
multi_combos  <- group_combos[grepl("/", group_combos)]         # >= 2 signatures positive
group_pal     <- c(Neg = "lightgrey",
                    setNames(brewer.pal(max(3, length(single_combos)), "Set1")[seq_along(single_combos)],
                             single_combos),
                    setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(multi_combos)),
                             multi_combos))[group_levels]

# ── Per-sample scoring + binarisation ─────────────────────────────────────────
srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
if (!file.exists(srt_f)) stop("Missing tumour_srt.qs2 for ", s)

message("Loading ", s, " ...")
srt <- qs_read(srt_f)

score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
DefaultAssay(srt) <- score_assay
if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

if (file.exists(meta_f)) {
  message(s, " already processed (", meta_f, " found) — reusing saved metadata")
  saved_meta <- readRDS(meta_f)
  srt  <- AddMetaData(srt, saved_meta)
  keep <- module_names %in% colnames(saved_meta)
  if (!any(keep)) stop("no scored module columns found in ", meta_f)
} else {
  # restrict each signature to genes present in this object (need >= 2 for scoring)
  feats <- lapply(groupdeg, function(g) intersect(g, rownames(srt)))
  keep  <- lengths(feats) >= 2
  if (!any(keep)) stop("no scorable signatures for ", s)

  srt <- AddModuleScore(srt, features = feats[keep], name = "GDmod_", assay = score_assay)
  new_cols <- paste0("GDmod_", seq_len(sum(keep)))                  # appended in feats[keep] order
  names(srt@meta.data)[match(new_cols, names(srt@meta.data))] <- module_names[keep]

  # binarise each per-cell module score (GMM background threshold) → pos/neg factor
  for (m in module_names[keep]) {
    sc  <- setNames(srt@meta.data[[m]], colnames(srt))
    message("  binarising ", s, " ", m)
    bin <- binarise_expression(sc, verbose = TRUE,
             plot_out = file.path(bin_dir_s, sprintf("%s_binarisation.png", m)))
    srt@meta.data[[paste0(m, "_pos")]] <-
      factor(ifelse(bin == 1L, "pos", "neg"), levels = c("neg", "pos"))
  }

  # combine the per-signature pos calls into one label, e.g. "DT_G1", "DT_G1/DT_G2", "Neg"
  kept_mods <- module_names[keep]
  labs      <- sub("^module_", "", kept_mods)
  pos_mat   <- do.call(cbind, lapply(kept_mods,
                  function(m) srt@meta.data[[paste0(m, "_pos")]] == "pos"))
  srt@meta.data$Module_group <- factor(apply(pos_mat, 1, function(r) {
    hit <- labs[which(r)]                           # drops NA + neg, keeps module order
    if (!length(hit)) "Neg" else paste(hit, collapse = "/")
  }), levels = group_levels)
}

# spatial map of the combined Module_group label
sp_group_s <- ImageDimPlot(srt, group.by = "Module_group", cols = group_pal,
                   border.color = "#00000000", size = 0.4) +
  ggtitle(s) +
  theme(plot.title = element_text(size = 8), legend.text = element_text(size = 6))
ggsave(file.path(png_dir_s, "sp_group.png"), sp_group_s,
       width = 6, height = 5, dpi = 150)

# spatial map of positive cells, one ImageDimPlot per signature
sp_pos_s <- list()
for (m in paste0(module_names[keep], "_pos")) {
  npos <- sum(srt@meta.data[[m]] == "pos", na.rm = TRUE)
  mod  <- sub("_pos$", "", m)
  p    <- ImageDimPlot(srt, group.by = m, cols = pos_pal,
                 border.color = "#00000000", size = 0.4) +
    ggtitle(sprintf("%s  %s+ (%d)", s, mod, npos)) +
    theme(plot.title = element_text(size = 8), legend.position = "none")
  sp_pos_s[[mod]] <- p
  ggsave(file.path(png_dir_s, sprintf("sp_pos_%s.png", mod)), p,
         width = 6, height = 5, dpi = 150)
}

# ── DEG by Module_group + pathway enrichment on those DEGs ─────────────────
srt$Module_group <- droplevels(srt$Module_group)
deg_mg_f <- file.path(deg_dir_s, "deg_modulegroup.Rds")
if (file.exists(deg_mg_f)) {
  message("  ", s, " Module_group DEG already detected (", deg_mg_f, ") — skipping")
} else if (dplyr::n_distinct(srt$Module_group) >= 2) {
  Idents(srt) <- "Module_group"
  DEG_mg <- tryCatch(
    FindAllMarkers(srt, assay = score_assay, test.use = "MAST",
                   only.pos = TRUE, verbose = FALSE),
    error = function(e) { message("  FindAllMarkers failed for ", s, ": ",
                                  conditionMessage(e)); NULL })
  if (!is.null(DEG_mg) && nrow(DEG_mg) > 0) {
    saveRDS(DEG_mg, deg_mg_f)
    write.csv(DEG_mg, file.path(deg_dir_s, "deg_modulegroup.csv"),
              row.names = FALSE)
    run_gsea_panel(DEG_mg, gene_sets,
                   file.path(deg_dir_s, "enrich_modulegroup.Rds"))
  }
}

# ── DEG per-signature (pos vs neg, one direction: genes up in pos only) ─────
for (m in module_names[keep]) {
  lab <- sub("^module_", "", m)
  col <- paste0(m, "_pos")
  deg_pos_f <- file.path(deg_dir_s, sprintf("deg_%s_pos.Rds", lab))
  if (file.exists(deg_pos_f)) {
    message("  ", s, " ", lab, " pos/neg DEG already detected (", deg_pos_f, ") — skipping")
  } else if (dplyr::n_distinct(srt@meta.data[[col]]) >= 2) {
    Idents(srt) <- col
    DEG_pos <- tryCatch(
      FindMarkers(srt, ident.1 = "pos", ident.2 = "neg", assay = score_assay,
                  test.use = "MAST", only.pos = TRUE, verbose = FALSE),
      error = function(e) { message("  FindMarkers failed for ", s, " ", lab, ": ",
                                    conditionMessage(e)); NULL })
    if (!is.null(DEG_pos) && nrow(DEG_pos) > 0) {
      DEG_pos <- tibble::rownames_to_column(DEG_pos, "gene")
      DEG_pos$cluster <- lab
      saveRDS(DEG_pos, deg_pos_f)
      write.csv(DEG_pos, file.path(deg_dir_s, sprintf("deg_%s_pos.csv", lab)),
                row.names = FALSE)
      run_gsea_panel(DEG_pos, gene_sets,
                     file.path(deg_dir_s, sprintf("enrich_%s_pos.Rds", lab)))
    }
  }
}

# ── Save per-sample outputs for the aggregate step ────────────────────────────
saveRDS(srt@meta.data, meta_f)
qs_save(sp_group_s, file.path(sample_dir, "sp_group.qs2"))
qs_save(sp_pos_s, file.path(sample_dir, "sp_pos.qs2"))

message("Done processing ", s)
