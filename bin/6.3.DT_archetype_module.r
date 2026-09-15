
# 6.3.DT_archetype_module.r
# DT-only counterpart of 6.3.archetype_module.r. Same derive-groupgene /
# correlation-filter workflow, but sourced from the 6.2.DT_archetype_downstream.R
# outputs (archetypes fit on category == "DT" cells only):
#   - archetype_expression.csv from 6.1.DT_archetype/{sample}/
#   - archetype_group_DE_MAST.Rds from 6.2.archetype_downstream_tumour_DT/
#   - per-sample DT_srt.qs2 (already SpaNorm-normalised on the DT subset, with
#     archetype metadata) from 6.2.archetype_downstream_tumour_DT/DT_srt/
# Output: VisHD/6.3.DT_archetype_module/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(jsonlite)
  library(qs2)
})
library(dendextend, lib.loc = "~/R_Library/4.5")
source("functions.R")
# ── Config ────────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "6.3.DT_archetype_module")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── DEG markers ───────────────────────────────────────────────────────────────
markerlist <- readRDS(file.path(base_dir,
  "6.2.archetype_downstream_tumour_DT", "archetype_group_DE_MAST.Rds"))
markers <- bind_rows(markerlist, .id = "sample")
markers$cluster    <- sub("^Archetype_", "A", markers$cluster)
markers$arch_group <- paste(markers$sample, markers$cluster, sep = "_")

# ── Per-archetype-group gene sets (significant, up, expression-enriched) ───────
arch_expr_list <- list()
for (s in samples) {
  expr_f <- file.path(base_dir, "6.1.DT_archetype", s, "archetype_expression.csv")
  if (!file.exists(expr_f)) {
    message("Skipping ", s, " — missing archetype_expression.csv")
    next
  }
  expr <- read.csv(expr_f, row.names = 1, check.names = FALSE)
  rownames(expr) <- paste0(s, "_A", rownames(expr))   # label rows {sample}_A{idx}
  arch_expr_list[[s]] <- expr
}
common_genes  <- Reduce(intersect, lapply(arch_expr_list, colnames))
common_genes  <- common_genes[!grepl("^MT-", common_genes)]
arch_expr_all <- do.call(rbind, lapply(unname(arch_expr_list), function(df) df[, common_genes]))

arch_group_gene <- sapply(unique(markers$arch_group), function(arch) {
  markergene <- markers %>%
    filter(arch_group == arch, p_val_adj < 0.05, avg_log2FC > 0) %>%
    pull(gene)
  enrichgene <- colnames(arch_expr_all)[arch_expr_all[arch, ] > 0.05]
  intersect(enrichgene, markergene)
})
arch_group_gene <- arch_group_gene[lengths(arch_group_gene) > 10]

# Pairwise overlap (containment) between archetype-group gene sets
mat <- sapply(arch_group_gene, function(x)
  sapply(arch_group_gene, function(y)
    length(intersect(x, y)) / min(length(x), length(y))))

# ── derive_groupgene: cluster archetype groups, take consensus genes per cluster
derive_groupgene <- function(arch_group_gene, mat, k, cutoff,
                             show_name = FALSE, split_file = NULL) {
  pdf(NULL)                                   # extract dendrogram off-screen
  ht <- draw(Heatmap(mat))
  dev.off()
  tree  <- cutree(column_dend(ht), k = k)
  group <- split(names(tree), tree)

  if (!is.null(split_file)) {
    png(split_file, width = 1400, height = 1200, res = 150)
    draw(Heatmap(mat, column_split = k, row_split = k,
                 show_row_names = show_name, show_column_names = show_name,
                 name = "overlap"))
    dev.off()
  }

  groupgene <- lapply(group, function(x)
    names(which(table(unlist(arch_group_gene[x])) >= length(x) * cutoff)))
  keep <- lengths(groupgene) > 0
  list(groupgene = groupgene[keep], group = group[keep])
}

# ── make_exclusive: a gene shared by multiple groups is kept only in the group
# where its mean avg_log2FC across that group's archetype members is highest ────
make_exclusive <- function(groupgene, group, markers) {
  dup <- names(which(table(unlist(groupgene)) > 1))
  for (g in dup) {
    cand  <- names(groupgene)[vapply(groupgene, function(v) g %in% v, logical(1))]
    score <- vapply(cand, function(grp)
      mean(markers$avg_log2FC[markers$gene == g &
                              markers$arch_group %in% group[[grp]]], na.rm = TRUE),
      numeric(1))
    for (grp in setdiff(cand, names(which.max(score))))
      groupgene[[grp]] <- setdiff(groupgene[[grp]], g)
  }
  groupgene[lengths(groupgene) > 0]
}

# ── Run ────────────────────────────────────────────────────────────────────────
png(file.path(outdir, "overlap_density.png"), width = 900, height = 600, res = 150)
plot(density(mat), main = "Archetype-group gene-set overlap", xlab = "containment")
abline(v = 0.5, col = "red", lty = 2)
dev.off()

png(file.path(outdir, "overlap_heatmap.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(mat, name = "overlap"))
dev.off()

res       <- derive_groupgene(arch_group_gene, mat, k = 4, cutoff = 0.4,
                              split_file = file.path(outdir, "group_split_heatmap.png"))
groupgene <- res$groupgene

shared_before <- names(which(table(unlist(groupgene)) > 1))
message(length(shared_before), " genes shared across groups before cleaning")

groupgene <- make_exclusive(groupgene, res$group, markers)
stopifnot(length(names(which(table(unlist(groupgene)) > 1))) == 0)
# groupgene <- groupgene[lengths(groupgene) >= 10]   # drop groups with < 10 genes
message("groups after cleaning: ", paste(names(groupgene), collapse = ", "),
        " (", paste(lengths(groupgene), collapse = ", "), " genes)")


# ── Second round: co-expression correlation filter ─────────────────────────────
# Correlate the (normalised) per-cell expression of every group gene within each
# sample and average Pearson r across samples. Visualise the gene-gene correlation
# matrix (annotated by group origin) BEFORE and AFTER two operations: (1) drop
# genes weakly co-expressed with their own group, (2) merge groups that are highly
# correlated across each other.
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

cor_cutoff   <- 0.02                             # min mean within-group correlation to keep
all_gg_genes <- unique(unlist(groupgene, use.names = FALSE))

cor_list <- list()
for (s in samples) {
  srt_f <- file.path(base_dir, "6.2.archetype_downstream_tumour_DT", "DT_srt",
                     sprintf("%s_DT_srt.qs2", s))
  if (!file.exists(srt_f)) next
  message("Correlating group genes for ", s, " ...")
  srt <- qs_read(srt_f)
  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  g <- intersect(all_gg_genes, rownames(srt))
  expr <- t(as.matrix(GetAssayData(srt, assay = score_assay, layer = "data")[g, , drop = FALSE]))
  cc <- suppressWarnings(cor(expr, method = "spearman"))   # NA cols where a gene has 0 variance
  cm <- matrix(NA_real_, length(all_gg_genes), length(all_gg_genes),
               dimnames = list(all_gg_genes, all_gg_genes))
  cm[rownames(cc), colnames(cc)] <- cc
  cor_list[[s]] <- cm
  rm(srt, expr, cc); gc()
}

# average across samples, ignoring NA (gene absent / zero-variance in a sample)
cor_sum <- Reduce(`+`, lapply(cor_list, function(m) { m[is.na(m)] <- 0; m }))
cor_n   <- Reduce(`+`, lapply(cor_list, function(m) (!is.na(m)) * 1L))
cor_avg <- cor_sum / cor_n
cor_avg[cor_n == 0] <- NA

# ── Heatmap helper: gene-gene correlation, ordered & annotated by group origin ──
plot_cor_heatmap <- function(gg, cmat, file, title) {
  gg  <- lapply(gg, function(x) intersect(x, rownames(cmat)))
  gg  <- gg[lengths(gg) > 0]
  ord <- unlist(gg, use.names = FALSE)
  gv  <- rep(names(gg), lengths(gg))
  ng  <- length(gg)
  cols <- setNames(
    if (ng <= 8) brewer.pal(max(3, ng), "Set2")[seq_len(ng)]
    else colorRampPalette(brewer.pal(8, "Set2"))(ng), names(gg))
  m <- cmat[ord, ord]
  m[is.na(m)] <- 0                       # NA (absent / zero-variance gene) → 0 so clustering works
  png(file, width = 1600, height = 1500, res = 150)
  ht <- draw(Heatmap(
    m, name = "Spearman r", column_title = title,
    col = colorRamp2(c(-1, 0, 1), c("steelblue", "white", "indianred")),
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation  = HeatmapAnnotation(group = gv, col = list(group = cols)),
    left_annotation = rowAnnotation(group = gv, col = list(group = cols))))
  dev.off()
  invisible(ht)
}

# BEFORE: original groups / genes (cross-sample average)
ht_bf <- plot_cor_heatmap(groupgene, cor_avg,
                 file.path(outdir, "groupgene_correlation_before.png"),
                 "Group-gene correlation (before filtering)")

# per-sample correlation heatmaps (original groups / genes)
persample_dir <- file.path(outdir, "correlation_per_sample")
dir.create(persample_dir, showWarnings = FALSE, recursive = TRUE)
for (s in names(cor_list))
  plot_cor_heatmap(groupgene, cor_list[[s]],
                   file.path(persample_dir, sprintf("groupgene_correlation_%s.png", s)),
                   paste0("Group-gene correlation — ", s))

# (3) cut the filter+merge heatmap's row dendrogram into 3 final correlation groups
gene_cl <- cutree(as.hclust(row_dend(ht_bf)), k = 4)
groupgene <- setNames(split(names(gene_cl), gene_cl),
                      paste0("G", seq_len(length(unique(gene_cl)))))
message("final 4 groups (", paste(lengths(groupgene), collapse = ", "), " genes)")

# (1) gene-level filter: drop genes with weak mean within-group correlation
groupgene <- setNames(lapply(names(groupgene), function(grp) {
  gs <- intersect(groupgene[[grp]], rownames(cor_avg))
  if (length(gs) < 2) return(groupgene[[grp]])
  sub <- cor_avg[gs, gs]; diag(sub) <- NA
  gs[rowMeans(sub, na.rm = TRUE) >= cor_cutoff]
}), names(groupgene))

message("after gene filter: ", paste(names(groupgene), collapse = ", "),
        " (", paste(lengths(groupgene), collapse = ", "), " genes)")

# (2) merge groups whose genes are highly correlated across groups
merge_cutoff <- 0.1                                 # min mean between-group r to merge
if (length(groupgene) > 1) {
  grp_names <- names(groupgene)
  G <- outer(grp_names, grp_names, Vectorize(function(a, b) {
    ga <- intersect(groupgene[[a]], rownames(cor_avg))
    gb <- intersect(groupgene[[b]], rownames(cor_avg))
    mean(cor_avg[ga, gb], na.rm = TRUE)
  }))
  dimnames(G) <- list(grp_names, grp_names)
  cl <- cutree(hclust(as.dist(1 - G), method = "average"), h = 1 - merge_cutoff)
  groupgene <- lapply(split(grp_names, cl), function(grps)
    unique(unlist(groupgene[grps], use.names = FALSE)))
  names(groupgene) <- paste0("G", seq_along(groupgene))
  message("after merge: ", paste(names(groupgene), collapse = ", "),
          " (", paste(lengths(groupgene), collapse = ", "), " genes)")
}

# AFTER: filtered genes / merged groups
ht_fm <- plot_cor_heatmap(groupgene, cor_avg,
                 file.path(outdir, "groupgene_correlation_after.png"),
                 "Group-gene correlation (after filter + merge)")

# (3) cut the filter+merge heatmap's row dendrogram into 5 final correlation groups
gene_cl <- cutree(as.hclust(row_dend(ht_fm)), k = 5)
groupgene <- setNames(split(names(gene_cl), gene_cl),
                      paste0("G", seq_len(length(unique(gene_cl)))))
message("final 5 groups (", paste(lengths(groupgene), collapse = ", "), " genes)")

plot_cor_heatmap(groupgene, cor_avg,
                 file.path(outdir, "groupgene_correlation_5groups.png"),
                 "Group-gene correlation (5 groups)")

# NOTE: manual curation step from 6.3.archetype_module.r (drop/merge specific
# G-clusters) was tuned to that script's dendrogram and is NOT reapplied here —
# inspect groupgene_correlation_6groups.png and adjust group membership by hand
# before trusting downstream results if curation is warranted.

plot_cor_heatmap(groupgene, cor_avg,
                 file.path(outdir, "groupgene_correlation_final.png"),
                 "Group-gene correlation (final)")

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(groupgene, file.path(outdir, "groupgene.Rds"))
gg_df <- data.frame(
  group = rep(names(groupgene), lengths(groupgene)),
  gene  = unlist(groupgene, use.names = FALSE)
)
write.csv(gg_df, file.path(outdir, "groupgene.csv"), row.names = FALSE)

# ── Visualization: groupgene module scores across all 8 samples ────────────────
# AddModuleScore each (now exclusive) group's genes per sample, then render
# FeaturePlot (UMAP) and the spatial map. DT_srt.qs2 are FOV-based Visium HD
# objects → SpatialFeaturePlot has no compatible image, so the spatial panel
# uses ImageFeaturePlot (the FOV equivalent used elsewhere in 6.2.1).
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
})

vizdir <- file.path(outdir, "module_score_spatial")
dir.create(vizdir, showWarnings = FALSE, recursive = TRUE)

# DEG + enrichment prerequisites (per-sample DEG runs inside the scoring loop below)
suppressPackageStartupMessages({ library(clusterProfiler) })
stopifnot(requireNamespace("MAST", quietly = TRUE))
source(file.path(base_dir, "functions.R"))        # pathwayenrich_plot()
Hall <- readRDS(file.path(base_dir, "Hall.Rds"))
C5   <- readRDS(file.path(base_dir, "C5.Rds"))
C6   <- readRDS(file.path(base_dir, "C6.Rds"))
deg_dir <- file.path(outdir, "group_DEG_enrichment")
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

module_names <- paste0("module_", names(groupgene))
feat_plots <- setNames(vector("list", length(groupgene)), module_names)  # UMAP FeaturePlot
sp_plots   <- setNames(vector("list", length(groupgene)), module_names)  # spatial
dim_plots  <- setNames(vector("list", length(groupgene)), module_names)  # UMAP top-30% highlight
score_list <- list()                                                     # per-cell module scores

for (s in samples) {
  srt_f <- file.path(base_dir, "6.2.archetype_downstream_tumour_DT", "DT_srt",
                     sprintf("%s_DT_srt.qs2", s))
  if (!file.exists(srt_f)) { message("Skipping ", s, " — missing DT_srt.qs2"); next }
  message("Scoring groupgene modules for ", s, " ...")
  srt <- qs_read(srt_f)

  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  # restrict each module to genes present in this object (need >= 2 for scoring)
  feats <- lapply(groupgene, function(g) intersect(g, rownames(srt)))
  keep  <- lengths(feats) >= 2
  if (!any(keep)) { rm(srt); gc(); next }

  srt <- AddModuleScore(srt, features = feats[keep], name = "GGmod_", assay = score_assay)
  new_cols <- paste0("GGmod_", seq_len(sum(keep)))                 # appended in feats[keep] order
  names(srt@meta.data)[match(new_cols, names(srt@meta.data))] <- module_names[keep]

  # collect the per-cell module scores for export
  sc           <- srt@meta.data[, module_names[keep], drop = FALSE]
  sc$barcode   <- rownames(sc)
  sc$slide     <- s
  score_list[[s]] <- sc

  s_dir <- file.path(deg_dir, "per_sample", s)            # per-sample DEG + enrichment
  dir.create(s_dir, showWarnings = FALSE, recursive = TRUE)

  # DT_srt.qs2 was normalised/clustered by do.spanorm_dt() (RunUMAP default
  # reduction name "umap"), not the cross-sample "pearsonumap" used elsewhere.
  red <- grep("umap", Reductions(srt), ignore.case = TRUE, value = TRUE)[1]
  for (m in module_names[keep]) {
    mid <- median(srt@meta.data[[m]], na.rm = TRUE)          # centre scale on module-score median
    if (!is.na(red))
      feat_plots[[m]][[s]] <- FeaturePlot(srt, features = m, reduction = red, order = TRUE) +
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
        ggtitle(paste0(s, "  ", m)) + theme(plot.title = element_text(size = 8))
    sp_plots[[m]][[s]] <- ImageFeaturePlot(srt, features = m, size = 0.4) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
      ggtitle(paste0(s, "  ", m)) + theme(plot.title = element_text(size = 8))

    # highlight the top 30% of cells by this module's score (>= 70th percentile)
    x   <- srt@meta.data[[m]]
    thr <- quantile(x, 0.9, na.rm = TRUE)
    hi  <- colnames(srt)[!is.na(x) & x >= thr]
    if (!is.na(red))
      dim_plots[[m]][[s]] <- DimPlot(srt, reduction = red, cells.highlight = hi,
                                     cols.highlight = "indianred", cols = "grey85",
                                     sizes.highlight = 0.4, order = TRUE) +
        NoLegend() +
        ggtitle(paste0(s, "  ", m)) + theme(plot.title = element_text(size = 8))

    # DEG: top-30% module cells vs rest (MAST) + GSEA, visualised when significant
    grp  <- sub("^module_", "", m)
    thrd <- quantile(x, 0.7, na.rm = TRUE)
    srt$group_id <- ifelse(!is.na(x) & x >= thrd, grp, "rest")
    Idents(srt)  <- "group_id"
    deg <- tryCatch(FindMarkers(srt, ident.1 = grp, ident.2 = "rest", assay = score_assay,
                                test.use = "MAST", only.pos = FALSE),
                    error = function(e) NULL)
    if (!is.null(deg)) {
      deg$gene <- rownames(deg)
      write.csv(deg, file.path(s_dir, sprintf("DEG_MAST_%s_vs_rest.csv", grp)),
                row.names = FALSE)
      sig <- deg %>% filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
        arrange(desc(avg_log2FC))
      if (nrow(sig) > 0) {
        gene_list <- setNames(sig$avg_log2FC, sig$gene)
        enr <- lapply(list(Hallmark = Hall, C6 = C6, C5 = C5), function(gs)
          tryCatch(clusterProfiler::GSEA(gene_list, TERM2GENE = gs, verbose = FALSE),
                   error = function(e) NULL, warning = function(w) NULL))
        for (nm in names(enr)) {
          res_e <- enr[[nm]]
          if (is.null(res_e) || nrow(res_e@result) == 0) next
          # persist the enrichment table (NES + p.adjust per pathway) for cross-sample heatmaps
          enr_cols <- intersect(c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust"),
                                colnames(res_e@result))
          write.csv(res_e@result[, enr_cols, drop = FALSE],
                    file.path(s_dir, sprintf("enrich_%s_%s.csv", grp, nm)), row.names = FALSE)
          sig_n <- sum(res_e@result$p.adjust < 0.05, na.rm = TRUE)
          if (sig_n == 0) next                        # only visualise significant enrichment
          ggsave(file.path(s_dir, sprintf("GSEA_%s_%s.pdf", grp, nm)),
                 pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res_e),
                 width = 6, height = 10)
        }
      }
    }
  }
  rm(srt); gc()
}

# ── Aggregate into one big grid: module in row, sample in column (modules × 8) ──
build_grid <- function(plotmat) {           # plotmat[[module]][[sample]]
  cells <- list()
  for (m in module_names) for (s in samples) {
    p <- plotmat[[m]][[s]]
    cells[[length(cells) + 1]] <- if (is.null(p)) patchwork::plot_spacer() else p
  }
  wrap_plots(cells, nrow = length(module_names), ncol = length(samples), byrow = TRUE)
}

grid_w <- length(samples) * 4
grid_h <- length(module_names) * 4

if (any(lengths(feat_plots) > 0))
  ggsave(file.path(vizdir, "featureplot_all.png"), build_grid(feat_plots),
         width = grid_w, height = grid_h, dpi = 150, limitsize = FALSE)
if (any(lengths(sp_plots) > 0))
  ggsave(file.path(vizdir, "spatial_all.png"), build_grid(sp_plots),
         width = grid_w, height = grid_h, dpi = 150, limitsize = FALSE)
if (any(lengths(dim_plots) > 0))
  ggsave(file.path(vizdir, "dimplot_top30_all.png"), build_grid(dim_plots),
         width = grid_w, height = grid_h, dpi = 150, limitsize = FALSE)

# Save per-cell module scores (one row per cell, NA where a module was skipped)
if (length(score_list) > 0) {
  scores_all <- dplyr::bind_rows(score_list)
  write.csv(scores_all, file.path(vizdir, "groupgene_module_scores.csv"), row.names = FALSE)
  saveRDS(score_list, file.path(vizdir, "groupgene_module_scores.Rds"))
}

message("\nDone. Module-score figures + scores saved to ", vizdir,
        " (per-sample DEG + enrichment under ", file.path(deg_dir, "per_sample"), ")")

# ── Cross-sample summary heatmaps: per-GGmod avg_log2FC ─────────────────────────
# From the per-sample DEG (module-gene avg_log2FC) collected in the scoring loop
# above, draw one heatmap per GGmod group: genes in rows, samples in columns.
# Colour scale steelblue–white–indianred (low/0/high).
ps_dir   <- file.path(deg_dir, "per_sample")
summ_dir <- file.path(deg_dir, "cross_sample_summary")
dir.create(summ_dir, showWarnings = FALSE, recursive = TRUE)

hm_col <- colorRamp2(c(-2, 0, 2), c("steelblue", "white", "indianred"))

# assemble a (gene|pathway) x sample matrix from a per-sample list of named vectors
build_mat <- function(vlist) {
  vlist <- vlist[lengths(vlist) > 0]
  if (length(vlist) == 0) return(NULL)
  rn <- sort(unique(unlist(lapply(vlist, names))))
  m  <- matrix(NA_real_, length(rn), length(vlist), dimnames = list(rn, names(vlist)))
  for (s in names(vlist)) m[names(vlist[[s]]), s] <- vlist[[s]]
  m
}

# rows ordered by mean value (clustering disabled — many NAs across samples)
save_hm <- function(m, file, title, legend, k = 1) {
  if (is.null(m) || nrow(m) < 2) return(invisible())
  grp <- sapply(strsplit(colnames(m), split = "_"), "[", 2)
  m <- m[order(rowMeans(m, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  show_rn <- nrow(m) <= 120                                  # row labels only when legible
  w <- min(30000, max(700, ncol(m) * 110 + 300))             # clamp below cairo's ~32767 px limit
  h <- min(30000, max(500, nrow(m) * 16 + 220))
  png(file, width = w, height = h, res = 150)
  ht <- draw(Heatmap(m, name = legend, col = hm_col, column_title = title, na_col = "grey90",
                show_row_names = show_rn, column_split = grp, border = TRUE, row_split = k,
                 row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8)))
  dev.off()
  return(ht)
}

mat_list <- list()
for (grp in names(groupgene)) {
  # avg_log2FC of this module's defining genes (top-30% module cells vs rest)
  l2f <- setNames(lapply(samples, function(s) {
    f <- file.path(ps_dir, s, sprintf("DEG_MAST_%s_vs_rest.csv", grp))
    if (!file.exists(f)) return(numeric(0))
    d <- read.csv(f)
    if (nrow(d) == 0) return(numeric(0))
    d <- d[d$p_val_adj <0.05, ]
    setNames(d$avg_log2FC, d$gene)
  }), samples)
  mat <- build_mat(l2f)
  colnames(mat) <- paste(colnames(mat), grp, sep = "_")
  mat_list[[grp]] <- mat
}

# cbind all per-group matrices on the union of genes; absent gene/sample → 0
mat_list  <- mat_list[!vapply(mat_list, is.null, logical(1))]
all_genes <- sort(unique(unlist(lapply(mat_list, rownames))))
mat <- do.call(cbind, lapply(mat_list, function(m) {
  full <- matrix(0, length(all_genes), ncol(m), dimnames = list(all_genes, colnames(m)))
  full[rownames(m), ] <- ifelse(is.na(m), 0, m)   # present gene → its avg_log2FC; NA (absent in sample) → 0
  full
}))

mat <- mat[apply(mat, 1, function(x) sum(x != 0))>10, ]

ht <- save_hm(mat, file.path(summ_dir, "avglogFC_all_groups.png"),
        "Module-gene avg_log2FC (top-30% cells vs rest, sig p_adj<0.05)", "avg_log2FC", k = 15)

message("\nDone. Cross-sample avg_log2FC heatmap saved to ", summ_dir,
        " (inspect the row dendrogram to hand-pick a groupdeg selection, as done",
        " via extract_group_genes() in 6.3.archetype_module.r, if a curated",
        " groupdeg module-score visualisation is wanted downstream)")


hc <- row_dend(ht)
extract_group_genes <- function(dend, idx_list){
  lapply(idx_list, function(x){
    unlist(sapply(x, function(y) labels(dend[[y]])))
  })
}

groupdeg <- extract_group_genes(hc, list(c(12, 13), c(1, 4, 5), c(6, 14, 15), c(7, 9), c(10, 11)))
names(groupdeg) <- c( "DT_G5", "DT_G4", "DT_G3", "DT_G2", "DT_G1")



saveRDS(groupdeg, file.path(summ_dir, "groupdeg.rds"))
ht2 <- save_hm(mat[unlist(groupdeg), ], file.path(summ_dir, "avglogFC_all_groups_selected.png"),
        "Selected DEGs for group genes", "avg_log2FC", k = 5)

rowdend <- row_dend(ht2)
lapply(rowdend, function(x) labels(x))

message("\nDone. Cross-sample avg_log2FC heatmap saved to ", summ_dir)

# ── Visualization: groupdeg module scores across all 8 samples ─────────────────
# AddModuleScore each (now exclusive) groupdeg group's genes per sample, then
# render FeaturePlot (UMAP) and the spatial map. These tumour_srt.qs2 are FOV-based
# Visium HD objects → SpatialFeaturePlot has no compatible image, so the spatial
# panel uses ImageFeaturePlot (the FOV equivalent used elsewhere in 6.2.1).
deg_module_names <- paste0("module_", names(groupdeg))
deg_feat_plots <- setNames(vector("list", length(groupdeg)), deg_module_names)  # UMAP FeaturePlot
deg_sp_plots   <- setNames(vector("list", length(groupdeg)), deg_module_names)  # spatial

for (s in samples) {
  srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(srt_f)) { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
  message("Scoring groupdeg modules for ", s, " ...")
  srt <- qs_read(srt_f)

  score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
  DefaultAssay(srt) <- score_assay
  if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

  # restrict each module to genes present in this object (need >= 2 for scoring)
  feats <- lapply(groupdeg, function(g) intersect(g, rownames(srt)))
  keep  <- lengths(feats) >= 2
  if (!any(keep)) { rm(srt); gc(); next }

  srt <- AddModuleScore(srt, features = feats[keep], name = "DEGmod_", assay = score_assay)
  new_cols <- paste0("DEGmod_", seq_len(sum(keep)))               # appended in feats[keep] order
  names(srt@meta.data)[match(new_cols, names(srt@meta.data))] <- deg_module_names[keep]

  # persist the per-cell groupdeg module scores alongside the per-sample DEG outputs
  s_dir <- file.path(ps_dir, s)
  dir.create(s_dir, showWarnings = FALSE, recursive = TRUE)
  sc         <- srt@meta.data[, deg_module_names[keep], drop = FALSE]
  sc$barcode <- rownames(sc)
  sc$slide   <- s
  write.csv(sc, file.path(s_dir, "groupdeg_module_scores.csv"), row.names = FALSE)

  red <- grep("pearsonumap", Reductions(srt), ignore.case = TRUE, value = TRUE)[1]
  for (m in deg_module_names[keep]) {
    mid <- median(srt@meta.data[[m]], na.rm = TRUE)          # centre scale on module-score median
    if (!is.na(red))
      deg_feat_plots[[m]][[s]] <- FeaturePlot(srt, features = m, reduction = red, order = TRUE) +
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
        ggtitle(paste0(s, "  ", m)) + theme(plot.title = element_text(size = 8))
    deg_sp_plots[[m]][[s]] <- dark_feature_plot(srt, features = m, size = 0.4) +
      ggtitle(paste0(s, "  ", m))
  }
  rm(srt); gc()
}

# Aggregate into one big grid: module in row, sample in column (modules × 8)
build_grid_deg <- function(plotmat) {           # plotmat[[module]][[sample]]
  cells <- list()
  for (m in deg_module_names) for (s in samples) {
    p <- plotmat[[m]][[s]]
    cells[[length(cells) + 1]] <- if (is.null(p)) patchwork::plot_spacer() else p
  }
  wrap_plots(cells, nrow = length(deg_module_names), ncol = length(samples), byrow = TRUE)
}

grid_w <- length(samples) * 4
grid_h <- length(deg_module_names) * 4
if (any(lengths(deg_feat_plots) > 0))
  ggsave(file.path(summ_dir, "groupdeg_featureplot_all.png"), build_grid_deg(deg_feat_plots),
         width = grid_w, height = grid_h, dpi = 150, limitsize = FALSE)
if (any(lengths(deg_sp_plots) > 0))
  ggsave(file.path(summ_dir, "groupdeg_spatial_all.png"), build_grid_deg(deg_sp_plots),
         width = grid_w, height = grid_h, dpi = 400, limitsize = FALSE, bg = "black")

message("\nDone. groupdeg module-score figures saved to ", summ_dir)
