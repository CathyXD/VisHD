
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(jsonlite)
  library(qs2)
  library(grid)
  library(dendextend, lib.loc = "~/R_Library/4.5")
  library(ggplotify, lib.loc = "~/R_Library/4.5")
  library(ggalluvial, lib.loc = "~/R_Library/4.5")
  library(patchwork)
})

# ── Config ────────────────────────────────────────────────────────────────────
samples <- c("LUT-245-07", "LUT-245-09", "LUT-245-10", "LUT-245-11",
             "LUT-245-15", "LUT-245-16", "LUT-245-17", "LUT-245-20")

base_dir <- "~/VisHD"
outdir   <- file.path(base_dir, "6.3.archetype_module_Jaccard")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── DEG markers ───────────────────────────────────────────────────────────────
markers <- read.csv(file.path(base_dir,
  "6.2.archetype_downstream_tumour/archetype_group_DE_MAST_all_samples.csv"))
markers$arch_group <- markers$archetype   # already "{slide}_A{n}", matching arch_expr_list row labels

# ── Per-archetype-group gene sets (significant, up, expression-enriched) ───────
arch_expr_list <- list()
for (s in samples) {
  expr_f <- file.path(base_dir, "6.1.archetype", s, "archetype_expression.csv")
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
arch_expr_enrich <- apply(arch_expr_all, 1, function(x) colnames(arch_expr_all)[x > 0])   # genes enriched in this archetype group
arch_expr_enrich <- arch_expr_enrich[lengths(arch_expr_enrich) > 10]

arch_group_gene <- sapply(unique(markers$arch_group), function(arch) {
  markergene <- markers %>%
    filter(arch_group == arch, p_val_adj < 0.05, avg_log2FC > 0, pct.1 >0.3) %>%
    pull(gene)
  markergene
  # enrichgene <- colnames(arch_expr_all)[arch_expr_all[arch, ] > 0.05]
  # intersect(enrichgene, markergene)
})
arch_group_DEG <- arch_group_gene[lengths(arch_group_gene) > 10]

for(g in names(arch_group_DEG)) arch_group_DEG[[g]] <- intersect(arch_group_DEG[[g]], arch_expr_enrich[[g]])


# Pairwise overlap (Jaccard) between archetype-group gene sets
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
DEGmat <- outer(seq_along(arch_group_DEG), seq_along(arch_group_DEG),
             Vectorize(function(i, j) jaccard(arch_group_DEG[[i]], arch_group_DEG[[j]])))
dimnames(DEGmat) <- list(names(arch_group_DEG), names(arch_group_DEG))

library(igraph); library(Matrix)
find_group <- function(J, thr = 0.10) {
  A <- J; A[A < thr] <- 0; diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
  set.seed(42)
  cl <- cluster_leiden(g, objective_function = "modularity",
                       resolution_parameter = 1.0, n_iterations = 10)
  return(membership(cl))
}

# ── derive_groupgene: cluster archetype groups, take consensus genes per cluster
derive_groupgene <- function(arch_group_gene, mat, cutoff,
                             show_name = T, split_file = NULL) {
  group <- find_group(mat, thr = 0.1)

  if (!is.null(split_file)) {
    png(split_file, width = 1400, height = 1200, res = 150)
    draw(Heatmap(mat, col = colorRamp2(c(0, 0.5, 1), c("white", "#de6c01", "#800026")), 
    column_split = group, row_split = group,
                 show_row_names = show_name, show_column_names = show_name,
                 name = "jaccard"))
    dev.off()
  }
  group <- split(names(group), group)
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

# ── DT DEG markers (merge DEGs + enriched genes, DT archetype track) ───────────
# Same merge pattern as arch_group_DEG above (markers filtered by p_val_adj/
# avg_log2FC/pct.1, intersected with expression-enriched genes), sourced from
# the DT-only archetype outputs (6.2.archetype_downstream_tumour_DT /
# 6.1.DT_archetype), as in 6.3.DT_archetype_module.r.
markerlist_DT <- readRDS(file.path(base_dir,
  "6.2.archetype_downstream_tumour_DT", "archetype_group_DE_MAST.Rds"))
markers_DT <- bind_rows(markerlist_DT, .id = "sample")
markers_DT$cluster    <- sub("^Archetype_", "A", markers_DT$cluster)
markers_DT$arch_group <- paste(markers_DT$sample, markers_DT$cluster, sep = "_")

arch_expr_list_DT <- list()
for (s in samples) {
  expr_f <- file.path(base_dir, "6.1.DT_archetype", s, "archetype_expression.csv")
  if (!file.exists(expr_f)) {
    message("Skipping ", s, " — missing archetype_expression.csv (DT)")
    next
  }
  expr <- read.csv(expr_f, row.names = 1, check.names = FALSE)
  rownames(expr) <- paste0(s, "_A", rownames(expr))   # label rows {sample}_A{idx}
  arch_expr_list_DT[[s]] <- expr
}
common_genes_DT  <- Reduce(intersect, lapply(arch_expr_list_DT, colnames))
common_genes_DT  <- common_genes_DT[!grepl("^MT-", common_genes_DT)]
arch_expr_all_DT <- do.call(rbind, lapply(unname(arch_expr_list_DT), function(df) df[, common_genes_DT]))
arch_expr_enrich_DT <- apply(arch_expr_all_DT, 1, function(x) colnames(arch_expr_all_DT)[x > 0])
arch_expr_enrich_DT <- arch_expr_enrich_DT[lengths(arch_expr_enrich_DT) > 10]

arch_group_gene_DT <- sapply(unique(markers_DT$arch_group), function(arch) {
  markergene <- markers_DT %>%
    filter(arch_group == arch, p_val_adj < 0.05, avg_log2FC > 0, pct.1 > 0.3) %>%
    pull(gene)
  markergene
})
arch_group_DEG_DT <- arch_group_gene_DT[lengths(arch_group_gene_DT) > 10]

for (g in names(arch_group_DEG_DT))
  arch_group_DEG_DT[[g]] <- intersect(arch_group_DEG_DT[[g]], arch_expr_enrich_DT[[g]])

# Pairwise overlap (Jaccard) between DT archetype-group gene sets
DEGmat_DT <- outer(seq_along(arch_group_DEG_DT), seq_along(arch_group_DEG_DT),
             Vectorize(function(i, j) jaccard(arch_group_DEG_DT[[i]], arch_group_DEG_DT[[j]])))
dimnames(DEGmat_DT) <- list(names(arch_group_DEG_DT), names(arch_group_DEG_DT))

# ── CB DEG markers (merge DEGs + enriched genes, CB archetype track) ───────────
# Same merge pattern as arch_group_DEG above, sourced from the CB-only archetype
# outputs (6.2.archetype_downstream_tumour_CB / 6.1.CB_archetype), mirroring the
# DT block above (6.3.archetype_module_Jaccard.r).
markerlist_CB <- readRDS(file.path(base_dir,
  "6.2.archetype_downstream_tumour_CB", "archetype_group_DE_MAST.Rds"))
markers_CB <- bind_rows(markerlist_CB, .id = "sample")
markers_CB$cluster    <- sub("^Archetype_", "A", markers_CB$cluster)
markers_CB$arch_group <- paste(markers_CB$sample, markers_CB$cluster, sep = "_")

arch_expr_list_CB <- list()
for (s in samples) {
  expr_f <- file.path(base_dir, "6.1.CB_archetype", s, "archetype_expression.csv")
  if (!file.exists(expr_f)) {
    message("Skipping ", s, " — missing archetype_expression.csv (CB)")
    next
  }
  expr <- read.csv(expr_f, row.names = 1, check.names = FALSE)
  rownames(expr) <- paste0(s, "_A", rownames(expr))   # label rows {sample}_A{idx}
  arch_expr_list_CB[[s]] <- expr
}
common_genes_CB  <- Reduce(intersect, lapply(arch_expr_list_CB, colnames))
common_genes_CB  <- common_genes_CB[!grepl("^MT-", common_genes_CB)]
arch_expr_all_CB <- do.call(rbind, lapply(unname(arch_expr_list_CB), function(df) df[, common_genes_CB]))
arch_expr_enrich_CB <- apply(arch_expr_all_CB, 1, function(x) colnames(arch_expr_all_CB)[x > 0])
arch_expr_enrich_CB <- arch_expr_enrich_CB[lengths(arch_expr_enrich_CB) > 10]

arch_group_gene_CB <- sapply(unique(markers_CB$arch_group), function(arch) {
  markergene <- markers_CB %>%
    filter(arch_group == arch, p_val_adj < 0.05, avg_log2FC > 0, pct.1 > 0.3) %>%
    pull(gene)
  markergene
})
arch_group_DEG_CB <- arch_group_gene_CB[lengths(arch_group_gene_CB) > 10]

for (g in names(arch_group_DEG_CB))
  arch_group_DEG_CB[[g]] <- intersect(arch_group_DEG_CB[[g]], arch_expr_enrich_CB[[g]])

# Pairwise overlap (Jaccard) between CB archetype-group gene sets
DEGmat_CB <- outer(seq_along(arch_group_DEG_CB), seq_along(arch_group_DEG_CB),
             Vectorize(function(i, j) jaccard(arch_group_DEG_CB[[i]], arch_group_DEG_CB[[j]])))
dimnames(DEGmat_CB) <- list(names(arch_group_DEG_CB), names(arch_group_DEG_CB))

# ── Run ────────────────────────────────────────────────────────────────────────
analyse_modules <- function(mods,
                            min_size    = 3,
                            max_degree  = 2,
                            drop_blocks = NULL,
                            plot        = FALSE,
                            verbose     = TRUE) {

  stopifnot(is.list(mods), length(mods) >= 2)
  if (is.null(names(mods))) names(mods) <- paste0("M", seq_along(mods))
  mods <- lapply(mods, function(x) unique(as.character(x[!is.na(x)])))

  all_genes <- sort(unique(unlist(mods, use.names = FALSE)))

  # binary membership matrix: genes x modules
  memb <- sapply(mods, function(m) as.integer(all_genes %in% m))
  rownames(memb) <- all_genes

  # degree = times of overlap; origin = which modules
  deg    <- rowSums(memb)
  origin <- apply(memb, 1, function(x) paste(names(mods)[x == 1], collapse = "&"))

  gene_tbl <- data.frame(gene = all_genes, degree = deg, origin = origin,
                         row.names = NULL, stringsAsFactors = FALSE)

  # disjoint blocks, largest first
  blocks <- split(gene_tbl$gene, gene_tbl$origin)
  blocks <- blocks[order(-lengths(blocks))]
  block_deg <- lengths(strsplit(names(blocks), "&", fixed = TRUE))

  # per-module composition
  composition <- lapply(mods, function(m)
    sort(table(gene_tbl$origin[gene_tbl$gene %in% m]), decreasing = TRUE))

  # filter to sub-modules for scoring
  keep <- block_deg <= max_degree &
          lengths(blocks) >= min_size &
          !(names(blocks) %in% drop_blocks)
  sub_mods <- blocks[keep]

  if (verbose) {
    cat("--- block sizes ---\n");  print(lengths(blocks))
    cat("\n--- degree ---\n");     print(table(gene_tbl$degree))
    cat("\n--- retained ---\n");   print(lengths(sub_mods))
    cat("\nkept", sum(lengths(sub_mods)), "of", length(all_genes), "genes\n")
  }

  if (plot) {
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      cm <- ComplexHeatmap::make_comb_mat(mods)
      print(ComplexHeatmap::UpSet(cm,
              comb_order = order(-ComplexHeatmap::comb_size(cm))))
    } else message("ComplexHeatmap not installed; skipping plot.")
  }

  invisible(list(gene_tbl    = gene_tbl,
                 memb        = memb,
                 blocks      = blocks,
                 composition = composition,
                 sub_mods    = sub_mods,
                 dropped     = blocks[!keep]))
}

colfun <- colorRamp2(c(0, 0.5, 1), c("white", "gold", "red2"))
png(file.path(outdir, "overlap_heatmap_DEG.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(DEGmat, name = "jaccard", col = colfun))
dev.off()

res       <- derive_groupgene(arch_group_DEG, DEGmat, cutoff = 0.3, split_file = file.path(outdir, "group_split_heatmap_DEG.png"))
DEGgroup <- res$group
DEGgroup <- DEGgroup[lengths(DEGgroup)>1]
   # 3 groups of archetype groups
groupgene <- res$groupgene[names(DEGgroup)]
names(groupgene) <- paste("WT", names(groupgene), sep = "-")

# ── DEGmat_DT: heatmap + derive_groupgene for DT archetype-group gene sets ──────
png(file.path(outdir, "overlap_heatmap_DEG_DT.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(DEGmat_DT, name = "jaccard", col = colfun))
dev.off()

res_DT      <- derive_groupgene(arch_group_DEG_DT, DEGmat_DT, cutoff = 0.3, split_file = file.path(outdir, "group_split_heatmap_DEG_DT.png"))
DEGgroup_DT <- res_DT$group
DEGgroup_DT <- DEGgroup_DT[lengths(DEGgroup_DT)>1]
   # 3 groups of archetype groups
groupgene_DT <- res_DT$groupgene[names(DEGgroup_DT)]
names(groupgene_DT) <- paste("DT", names(groupgene_DT), sep = "-")

# ── DEGmat_CB: heatmap + derive_groupgene for CB archetype-group gene sets ──────
png(file.path(outdir, "overlap_heatmap_DEG_CB.png"), width = 1400, height = 1200, res = 150)
draw(Heatmap(DEGmat_CB, name = "jaccard", col = colfun))
dev.off()

res_CB      <- derive_groupgene(arch_group_DEG_CB, DEGmat_CB, cutoff = 0.3, split_file = file.path(outdir, "group_split_heatmap_DEG_CB.png"))
DEGgroup_CB <- res_CB$group
DEGgroup_CB <- DEGgroup_CB[lengths(DEGgroup_CB)>1]
   # groups of archetype groups
groupgene_CB <- res_CB$groupgene[names(DEGgroup_CB)]
names(groupgene_CB) <- paste("CB", names(groupgene_CB), sep = "-")

total_groupgene <- c(groupgene, groupgene_DT, groupgene_CB)
# res <- analyse_modules(total_groupgene, min_size = 3, max_degree = 2, drop_blocks = NULL,
#                 plot = FALSE, verbose = TRUE)

## Module gene lists -- CB (pre-treatment) and DT (post-treatment)
## WT modules excluded. CB-4 / CB-7 keep their original NMF factor names here;
## cb_vs_dt_module_comparison.R renames them to CB-1 / CB-2 in section 1.

groupgene <- list(

  ## DT-1  (n = 18)
  `DT-1` = c(
    "ACP3", "ACSL3", "AZGP1", "CKB", "FASN", "KLK2", "KLK3", "KLK4", "LENG8",
    "MLPH", "MSMB", "NDRG1", "NPDC1", "PDE4D", "PLPP1", "SAT1", "SLC39A6",
    "TMPRSS2"
  ),

  ## DT-2  (n = 25)
  `DT-2` = c(
    "A2M", "ACTB", "ACTG2", "CALD1", "COL1A1", "COL1A2", "COL3A1", "COL6A2",
    "COL6A3", "CSRP1", "DES", "EEF2", "FBLN1", "IGFBP7", "MGP", "MSMB",
    "MYH11", "MYL9", "PSAP", "SPARC", "TAGLN", "TIMP1", "TMSB4X", "TPM2",
    "VIM"
  ),

  ## DT-3  (n = 48)
  `DT-3` = c(
    "ACP3", "APLP2", "ATP2C1", "B2M", "BCAM", "CD74", "CKB", "CPE", "DHRS7",
    "F3", "FOLH1", "FTH1", "GOLM1", "GP2", "HSD17B4", "HSP90AA1", "KLK3",
    "KLK4", "MSMB", "MUC12", "NDRG1", "NEDD4L", "NPY", "NR4A1", "NUPR1",
    "PDE4D", "PDIA3", "PDK4", "PLA2G2A", "PLA2G7", "PPM1B", "PPP3CA",
    "PRUNE2", "PSAP", "RDH11", "REXO2", "RHOU", "SCUBE2", "SELENOW",
    "SERINC5", "SLC44A4", "STEAP2", "TAPBP", "TCIM", "TMSB10", "VEGFA",
    "XBP1", "ZBTB16"
  ),

  ## DT-4  (n = 33)
  `DT-4` = c(
    "ACTB", "ATF3", "CDKN1A", "CEBPD", "COL3A1", "DNAJA4", "DNAJB1", "DUSP1",
    "EGR1", "ERRFI1", "FASN", "FOS", "FOSB", "GDF15", "HSP90AA1", "HSPA1B",
    "HSPA8", "IER2", "JUN", "JUNB", "JUND", "KLF6", "MCL1", "MGP", "NR4A1",
    "PDE4D", "RHOB", "SGK1", "TAGLN", "VEGFA", "VMP1", "ZBTB16", "ZFP36"
  ),

  ## DT-5  (n = 67)
  `DT-5` = c(
    "ABCC4", "ACP3", "ADGRF5", "ALDH1A3", "AMACR", "APP", "ASS1", "ATP5MC2",
    "BCAM", "CACNA1D", "CCNG1", "CCNI", "CD24", "CD44", "CD9", "CIRBP",
    "EEF1B2", "EEF1G", "EEF2", "EIF1", "EIF3A", "ELOVL5", "EPCAM", "FASN",
    "FOXA1", "FTL", "GDF15", "GOLM1", "GPX4", "HDLBP", "HEBP2", "HSP90AB1",
    "HSPA8", "KLK12", "KLK2", "KLK3", "KRT18", "MBOAT2", "NCAPD3", "NDUFA4",
    "NEDD4L", "NKX3-1", "NME2", "NME4", "NPDC1", "NPY", "PPDPF", "PPP3CA",
    "PRAC1", "RACK1", "REXO2", "RGS10", "SCARB2", "SLC25A6", "SLC38A11",
    "SPOCK1", "TM9SF3", "TMPRSS2", "TOMM20", "TOMM7", "TRGC2", "TTC3",
    "TUT7", "UBA52", "XBP1", "YBX3", "YWHAE"
  ),

  ## CB-4  (n = 58)
  `CB-4` = c(
    "ABCC4", "ACP3", "ALDH1A3", "AMACR", "ATP2C1", "AZGP1", "B2M", "BCAM",
    "BRI3", "CD9", "CKB", "DDX5", "DHCR24", "DHRS7", "EEF1G", "EEF2",
    "ELAPOR1", "ELK4", "ERG", "F3", "FASN", "FXYD3", "GOLM1", "HERPUD1",
    "KLK2", "KLK4", "LCP1", "LRIG1", "MSMB", "NKX3-1", "NPDC1", "NPY",
    "NUPR1", "OAZ1", "OGT", "OR51E2", "RACK1", "RDH11", "REXO2", "SLC30A4",
    "SLC39A6", "SLC44A4", "SLC45A3", "SMS", "SPINT2", "SPON2", "ST6GAL1",
    "STEAP2", "TACSTD2", "TAGLN", "TMPRSS2", "TMSB4X", "TPD52", "TSPAN1",
    "TUT7", "UBA52", "XBP1", "YBX3"
  ),

  ## CB-7  (n = 19)
  `CB-7` = c(
    "ACTB", "AMACR", "ANPEP", "AZGP1", "CCNI", "EEF1G", "EEF2", "KLK2",
    "KLK3", "KLK4", "LCP1", "NAAA", "NKX3-1", "NPDC1", "SELENOP", "SLC15A2",
    "TFF3", "TMSB4X", "TRPM8"
  )
)

## Consistent palette across every panel.
## Purple = pre-treatment, coral = post-treatment, grey = lost / absent.
COL_CB   <- "#7F77DD"
COL_DT   <- "#D85A30"
COL_LOST <- "#888780"

stopifnot(exists("groupgene"), is.list(groupgene))


## -----------------------------------------------------------------------------
## 1. Rename CB modules to sequential order, drop WT, fix column order
## -----------------------------------------------------------------------------

## CB modules came out of the NMF as factors 4 and 7. Renaming them CB-1/CB-2
## for readability, but the mapping is written to disk so the original factor
## index is always recoverable.
cb_map <- c("CB-4" = "CB-1", "CB-7" = "CB-2")

idx <- match(names(cb_map), names(groupgene))
if (any(!is.na(idx))) {
  names(groupgene)[idx[!is.na(idx)]] <- cb_map[!is.na(idx)]
}

cb_key <- data.frame(new_label  = unname(cb_map),
                     original   = names(cb_map),
                     nmf_factor = c(4L, 7L))
write.csv(cb_key, file.path(outdir, "cb_module_rename_key.csv"), row.names = FALSE)

CB   <- c("CB-1", "CB-2")
DT   <- paste0("DT-", 1:5)
mods <- groupgene[c(CB, DT)]          # WT modules dropped here

## Descriptive labels. Index numbering is arbitrary WITHIN each arm and carries
## no cross-arm correspondence -- CB-1's closest relative is DT-5, not DT-1.
## Use these for figures; keep the CB-n / DT-n codes for tables and code.
mod_label <- c(
  "CB-1" = "ERG-luminal",
  "CB-2" = "Benign-like luminal",
  "DT-1" = "Secretory (KLK-high)",
  "DT-2" = "Stroma / smooth muscle",
  "DT-3" = "PSMA-high",
  "DT-4" = "Stress / AP-1",
  "DT-5" = "Metabolic / translational")

arm <- factor(ifelse(names(mods) %in% CB, "Pre (CB)", "Post (DT)"),
              levels = c("Pre (CB)", "Post (DT)"))

message("Module sizes:")
print(lengths(mods))


## -----------------------------------------------------------------------------
## 2. Overlap statistics
## -----------------------------------------------------------------------------

gene_pool <- sort(unique(unlist(mods)))
n_pool    <- length(gene_pool)

## Membership frequency: how many modules contain each gene. Genes in >=4
## modules are the shared prostate-luminal / translation backbone. They inflate
## every epithelial Jaccard and are down-weighted in the weighted metric below.
gene_freq <- sapply(gene_pool, function(g) sum(sapply(mods, function(m) g %in% m)))
backbone  <- names(gene_freq)[gene_freq >= 4]
message("Backbone genes (>=4 modules): ", paste(backbone, collapse = ", "))

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## Frequency-weighted overlap (IDF-style). Down-weights genes that appear in
## many modules, so similarity reflects identity rather than shared backbone.
idf <- log(length(mods) / gene_freq)
weighted_overlap <- function(a, b) {
  sh <- intersect(a, b); un <- union(a, b)
  if (!length(un)) return(0)
  sum(idf[sh]) / sum(idf[un])
}

## Hypergeometric enrichment: is the overlap larger than expected given the
## detected-gene background? Set bg to the number of genes tested in the NMF.
bg <- n_pool   # replace with nrow(obj) for a proper background
hyper_p <- function(a, b) {
  phyper(length(intersect(a, b)) - 1, length(a), bg - length(a), length(b),
         lower.tail = FALSE)
}

grid_pairs <- expand.grid(CB = CB, DT = DT, stringsAsFactors = FALSE)
overlap_tbl <- grid_pairs |>
  rowwise() |>
  mutate(
    n_shared   = length(intersect(mods[[CB]], mods[[DT]])),
    jaccard    = jaccard(mods[[CB]], mods[[DT]]),
    weighted_J = weighted_overlap(mods[[CB]], mods[[DT]]),
    p_hyper    = hyper_p(mods[[CB]], mods[[DT]]),
    shared_genes = paste(sort(intersect(mods[[CB]], mods[[DT]])), collapse = ";")
  ) |>
  ungroup() |>
  mutate(fdr = p.adjust(p_hyper, method = "BH"))

write.csv(overlap_tbl, file.path(outdir, "CBvsDT_overlap_statistics.csv"),
          row.names = FALSE)
print(overlap_tbl |> select(CB, DT, n_shared, jaccard, weighted_J, fdr))


## -----------------------------------------------------------------------------
## 3. Gene flow CB -> DT
## -----------------------------------------------------------------------------

dt_all <- unique(unlist(mods[DT]))

flow <- lapply(CB, function(cb) {
  g <- mods[[cb]]
  tibble(from = cb,
         to   = c(DT, "No DT module"),
         n    = c(sapply(DT, function(d) sum(g %in% mods[[d]])),
                  sum(!g %in% dt_all)))
}) |> bind_rows() |> filter(n > 0)

## NOTE: genes that enter more than one DT module are counted more than once,
## so band totals exceed module size. That multi-counting IS the fission signal,
## but it must be stated in the figure legend.
multi <- sapply(CB, function(cb)
  sum(sapply(mods[[cb]], function(g) sum(sapply(mods[DT], function(m) g %in% m))) > 1))
message("Genes entering >1 DT module -- ",
        paste(sprintf("%s: %d", CB, multi), collapse = "; "))

write.csv(flow, file.path(outdir, "CBtoDT_gene_flow.csv"), row.names = FALSE)


## -----------------------------------------------------------------------------
## 4. Panel A : alluvial ribbon
## -----------------------------------------------------------------------------

flow_plot <- flow |>
  mutate(from_lab = paste0(from, "\n", mod_label[from]),
         to_lab   = ifelse(to == "No DT module", "No DT module",
                           paste0(to, "\n", mod_label[to])))

lvl_from <- unique(flow_plot$from_lab)
lvl_to   <- c(paste0(DT, "\n", mod_label[DT]), "No DT module")

p_ribbon <- ggplot(flow_plot,
                   aes(axis1 = factor(from_lab, lvl_from),
                       axis2 = factor(to_lab,   lvl_to),
                       y = n)) +
  geom_alluvium(aes(fill = from), alpha = .45, width = .18, knot.pos = .4) +
  geom_stratum(width = .18, fill = "grey96", colour = "grey40", linewidth = .3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.4, lineheight = .9) +
  scale_x_discrete(limits = c("Pre (CB)", "Post (DT)"), expand = c(.15, .05)) +
  scale_fill_manual(values = setNames(c(COL_CB, "#1D9E75"), CB)) +
  labs(y = "Shared genes", x = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())


## -----------------------------------------------------------------------------
## 5. Panel B : marker membership dot matrix
## -----------------------------------------------------------------------------
## Ordered marker panel. Block order and within-block gene order are both
## deliberate -- do not sort. Tier assignment drives the side annotation.

markers <- list(
  ## ---- TIER I : tumour identity & therapeutic targets ----------------------
  ## PSMA first; it is itself AR-regulated, hence the block directly below.
  "Radioligand / ADC target"  = c("FOLH1", "TACSTD2"),
  ## Pioneer/lineage TFs, then AR target genes, then AR-regulated transporters.
  "AR-regulated / luminal"    = c("FOXA1", "NKX3-1",
                                  "KLK3", "KLK4", "TMPRSS2", "SLC45A3", "AZGP1",
                                  "STEAP2", "ABCC4", "AMACR"),
  "ERG-fusion clone"          = c("ERG", "TPD52", "OR51E2"),
  "Epithelial identity"       = c("EPCAM", "KRT18", "CD24", "CD9"),
  "Basal / stem-associated"   = c("CD44", "ALDH1A3"),
  ## Expected to render empty -- NE markers first, then basal.
  "Lineage plasticity"        = c("ASCL1", "INSM1", "CHGA", "SYP", "DLL3",
                                  "POU2F3", "SOX2", "KRT5", "TP63"),
  "Benign / normal-adjacent"  = c("TRPM8", "ANPEP", "SLC15A2", "SELENOP",
                                  "NAAA", "MSMB"),
  
  ## ---- TIER II : treatment response ----------------------------------------
  ## Ordered by proximity to the radiation insult.
  "p53 / DNA damage"          = c("CDKN1A", "CCNG1", "GDF15"),
  ## AP-1 family, then other IEG TFs, then negative-feedback effectors.
  "AP-1 / immediate-early"    = c("FOS", "JUN", "JUNB", "JUND", "ATF3",
                                  "EGR1", "IER2", "ZFP36", "DUSP1"),
  "Heat shock / proteotoxic"  = c("HSPA1B", "HSP90AA1", "DNAJB1", "DNAJA4"),
  "ER stress / UPR"           = c("XBP1", "NUPR1"),
  ## PDK4 placed last so it abuts the metabolic tier below.
  "Glucocorticoid response"   = c("ZBTB16", "SGK1", "KLF6", "CEBPD", "PDK4"),
  
  ## ---- TIER III : adaptive metabolism --------------------------------------
  "Hypoxia / vascular"        = c("VEGFA", "NDRG1"),
  ## Synthesis -> elongation -> membrane remodelling -> peroxidation defence.
  "Lipid / redox metabolism"  = c("FASN", "ELOVL5", "MBOAT2", "GPX4", "ASS1"),
  
  ## ---- TIER IV : microenvironment ------------------------------------------
  ## MHC-I loading complex first, then MHC-II invariant chain.
  "Antigen presentation"      = c("B2M", "TAPBP", "PDIA3", "CD74"),
  "Fibroblast / ECM"          = c("COL1A1", "COL3A1", "FBLN1", "SPARC",
                                  "TIMP1", "VIM"),
  "Smooth muscle"             = c("ACTG2", "MYH11", "DES"))

tier <- c(
  rep("I. Identity & target",  7),
  rep("II. Treatment response", 5),
  rep("III. Adaptive metabolism", 2),
  rep("IV. Microenvironment",   3))
names(tier) <- names(markers)

## Ordered factors so ComplexHeatmap respects the sequence above.
gvec       <- unlist(markers, use.names = FALSE)
block_f    <- factor(rep(names(markers), lengths(markers)), levels = names(markers))
tier_f     <- factor(tier[as.character(block_f)], levels = unique(tier))
tier_cols <- c("I. Identity & target"     = "#7F77DD",
               "II. Treatment response"   = "#D85A30",
               "III. Adaptive metabolism" = "#1D9E75",
               "IV. Microenvironment"     = "#888780")

mat <- sapply(mods, function(m) as.integer(gvec %in% m))
rownames(mat) <- gvec
mat[, arm == "Post (DT)"] <- mat[, arm == "Post (DT)"] * 2L

## --- top annotation: module size, split by whether the gene appears in the panel

n_total <- lengths(mods)                                  # genes per module
n_shown <- sapply(mods, function(m) sum(m %in% gvec))     # of those, in the panel
size_mat <- cbind(shown = n_shown, hidden = n_total - n_shown)

stopifnot(identical(rownames(size_mat), colnames(mat)))   # order must match

top_anno <- HeatmapAnnotation(
  ## arm strip sits directly above the columns so pre/post reads without
  ## consulting the column split title
  arm  = arm,
  size = anno_barplot(
    size_mat,
    gp          = gpar(fill = c("grey35", "grey85"), col = NA),
    bar_width   = 0.75,
    height      = unit(14, "mm"),
    add_numbers = TRUE,                                   # total n above each bar
    numbers_gp  = gpar(fontsize = 6),
    axis_param  = list(gp = gpar(fontsize = 6))),
  col = list(arm = c("Pre (CB)" = COL_CB, "Post (DT)" = COL_DT)),
  annotation_name_gp   = gpar(fontsize = 7),
  annotation_name_side = "left",
  show_legend          = c(arm = FALSE),
  gap = unit(1.5, "mm"))

ht_markers <- Heatmap(mat,
        col              = c("0" = "grey96", "1" = COL_CB, "2" = COL_DT),
        row_split        = block_f,          # ordered factor -> keeps your sequence
        row_order        = seq_along(gvec),  # and keeps within-block gene order
        column_split     = arm,
        top_annotation   = top_anno,
        cluster_rows     = FALSE, cluster_columns = FALSE,
        row_title_rot    = 0, row_title_gp = gpar(fontsize = 9),
        row_names_gp     = gpar(fontsize = 6.5, fontface = "italic"),
        column_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 9),
        rect_gp          = gpar(col = "white", lwd = 0.8),
        left_annotation  = rowAnnotation(
          tier = tier_f,
          col  = list(tier = tier_cols),
          show_annotation_name = FALSE,
          annotation_legend_param = list(tier = list(title = NULL,
                                                     labels_gp = gpar(fontsize = 9)))),
        show_heatmap_legend = FALSE)


## -----------------------------------------------------------------------------
## 6. Combined figure
## -----------------------------------------------------------------------------

p_dots <- as.ggplot(grid.grabExpr(draw(ht_markers)))

fig <- (p_ribbon | p_dots) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

## 180 mm = Nature Communications double-column width.
ggsave(file.path(outdir, "Fig_CB_vs_DT_modules.pdf"), fig,
       width = 180, height = 150, units = "mm", device = cairo_pdf)
ggsave(file.path(outdir, "Fig_CB_vs_DT_modules.png"), fig,
       width = 180, height = 150, units = "mm", dpi = 400)


dt_all  <- unique(unlist(mods[grep("DT", names(mods))]))
orphans <- setdiff(union(mods[["CB-1"]], mods[["CB-2"]]), dt_all)

merged <- c(mods[DT], list("CB-only" = sort(orphans)))

## disjointness is a property of the construction -- assert it so a future
## edit to the module lists cannot silently break it
stopifnot(all(sapply(mods[DT], function(m) length(intersect(m, orphans))) == 0))

merge_key <- data.frame(
  module      = names(merged),
  n_genes     = lengths(merged),
  arm_origin  = c(rep("post (DT)", 5), "pre (CB) only"),
  source      = c(DT, "CB-1 (21) + CB-2 (6) + both (1)"))
write.csv(merge_key, file.path(outdir, "merged_module_key.csv"), row.names = FALSE)

print(lengths(merged))


groupgene <- mods

# # ── Second round: co-expression correlation filter ─────────────────────────────
# # Correlate the (normalised) per-cell expression of every group gene within each
# # sample and average Pearson r across samples. Visualise the gene-gene correlation
# # matrix (annotated by group origin) BEFORE and AFTER two operations: (1) drop
# # genes weakly co-expressed with their own group, (2) merge groups that are highly
# # correlated across each other.
# suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

# cor_cutoff   <- 0.02                             # min mean within-group correlation to keep
# all_gg_genes <- unique(unlist(groupgene, use.names = FALSE))

# cor_list <- list()
# for (s in samples) {
#   srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
#   if (!file.exists(srt_f)) next
#   message("Correlating group genes for ", s, " ...")
#   srt <- qs_read(srt_f)
#   score_assay <- if ("SpaNorm" %in% Assays(srt)) "SpaNorm" else "Spatial"
#   DefaultAssay(srt) <- score_assay
#   if (score_assay == "Spatial") srt <- NormalizeData(srt, verbose = FALSE)

#   g <- intersect(all_gg_genes, rownames(srt))
#   expr <- t(as.matrix(GetAssayData(srt, assay = score_assay, layer = "data")[g, , drop = FALSE]))
#   cc <- suppressWarnings(cor(expr, method = "spearman"))   # NA cols where a gene has 0 variance
#   cm <- matrix(NA_real_, length(all_gg_genes), length(all_gg_genes),
#                dimnames = list(all_gg_genes, all_gg_genes))
#   cm[rownames(cc), colnames(cc)] <- cc
#   cor_list[[s]] <- cm
#   rm(srt, expr, cc); gc()
# }

# # average across samples, ignoring NA (gene absent / zero-variance in a sample)
# cor_sum <- Reduce(`+`, lapply(cor_list, function(m) { m[is.na(m)] <- 0; m }))
# cor_n   <- Reduce(`+`, lapply(cor_list, function(m) (!is.na(m)) * 1L))
# cor_avg <- cor_sum / cor_n
# cor_avg[cor_n == 0] <- NA

# # ── Heatmap helper: gene-gene correlation, ordered & annotated by group origin ──
# plot_cor_heatmap <- function(gg, cmat, file, title) {
#   gg  <- lapply(gg, function(x) intersect(x, rownames(cmat)))
#   gg  <- gg[lengths(gg) > 0]
#   ord <- unlist(gg, use.names = FALSE)
#   gv  <- rep(names(gg), lengths(gg))
#   ng  <- length(gg)
#   cols <- setNames(
#     if (ng <= 8) brewer.pal(max(3, ng), "Set2")[seq_len(ng)]
#     else colorRampPalette(brewer.pal(8, "Set2"))(ng), names(gg))
#   m <- cmat[ord, ord]
#   m[is.na(m)] <- 0                       # NA (absent / zero-variance gene) → 0 so clustering works
#   png(file, width = 1600, height = 1500, res = 150)
#   ht <- draw(Heatmap(
#     m, name = "Spearman r", column_title = title,
#     col = colorRamp2(c(-1, 0, 0.5), c("steelblue", "white", "indianred")),
#     show_row_names = FALSE, show_column_names = FALSE,
#     top_annotation  = HeatmapAnnotation(group = gv, col = list(group = cols)),
#     left_annotation = rowAnnotation(group = gv, col = list(group = cols))))
#   dev.off()
#   invisible(ht)
# }

# # BEFORE: original groups / genes (cross-sample average)
# ht_bf <- plot_cor_heatmap(groupgene, cor_avg,
#                  file.path(outdir, "groupgene_correlation_before.png"),
#                  "Group-gene correlation (before filtering)")

# # per-sample correlation heatmaps (original groups / genes)
# persample_dir <- file.path(outdir, "correlation_per_sample")
# dir.create(persample_dir, showWarnings = FALSE, recursive = TRUE)
# for (s in names(cor_list))
#   plot_cor_heatmap(groupgene, cor_list[[s]],
#                    file.path(persample_dir, sprintf("groupgene_correlation_%s.png", s)),
#                    paste0("Group-gene correlation — ", s))

# # (3) cut the filter+merge heatmap's row dendrogram into 4 final correlation groups
# gene_cl <- cutree(as.hclust(row_dend(ht_bf)), k = 6)
# groupgene <- setNames(split(names(gene_cl), gene_cl),
#                       paste0("G", seq_len(length(unique(gene_cl)))))
# message("final 6 groups (", paste(lengths(groupgene), collapse = ", "), " genes)")

# # (1) gene-level filter: drop genes with weak mean within-group correlation
# groupgene <- setNames(lapply(names(groupgene), function(grp) {
#   gs <- intersect(groupgene[[grp]], rownames(cor_avg))
#   if (length(gs) < 2) return(groupgene[[grp]])
#   sub <- cor_avg[gs, gs]; diag(sub) <- NA
#   gs[rowMeans(sub, na.rm = TRUE) >= cor_cutoff]
# }), names(groupgene))
# groupgene <- groupgene[lengths(groupgene) >= 10]
# message("after gene filter: ", paste(names(groupgene), collapse = ", "),
#         " (", paste(lengths(groupgene), collapse = ", "), " genes)")

# # (2) merge groups whose genes are highly correlated across groups
# merge_cutoff <- 0.1                                 # min mean between-group r to merge
# if (length(groupgene) > 1) {
#   grp_names <- names(groupgene)
#   G <- outer(grp_names, grp_names, Vectorize(function(a, b) {
#     ga <- intersect(groupgene[[a]], rownames(cor_avg))
#     gb <- intersect(groupgene[[b]], rownames(cor_avg))
#     mean(cor_avg[ga, gb], na.rm = TRUE)
#   }))
#   dimnames(G) <- list(grp_names, grp_names)
#   cl <- cutree(hclust(as.dist(1 - G), method = "average"), h = 1 - merge_cutoff)
#   groupgene <- lapply(split(grp_names, cl), function(grps)
#     unique(unlist(groupgene[grps], use.names = FALSE)))
#   names(groupgene) <- paste0("G", seq_along(groupgene))
#   message("after merge: ", paste(names(groupgene), collapse = ", "),
#           " (", paste(lengths(groupgene), collapse = ", "), " genes)")
# }

# # AFTER: filtered genes / merged groups
# ht_fm <- plot_cor_heatmap(groupgene, cor_avg,
#                  file.path(outdir, "groupgene_correlation_after.png"),
#                  "Group-gene correlation (after filter + merge)")

# # (3) cut the filter+merge heatmap's row dendrogram into 4 final correlation groups
# gene_cl <- cutree(as.hclust(row_dend(ht_fm)), k = 5)
# groupgene <- setNames(split(names(gene_cl), gene_cl),
#                       paste0("G", seq_len(length(unique(gene_cl)))))
# message("final 5 groups (", paste(lengths(groupgene), collapse = ", "), " genes)")

# plot_cor_heatmap(groupgene, cor_avg,
#                  file.path(outdir, "groupgene_correlation_5groups.png"),
#                  "Group-gene correlation (5 groups)")

# # manual curation: drop G1, merge G3 + G4, then renumber groups in order
# groupgene$G3 <- unique(c(groupgene$G3, groupgene$G4))
# groupgene$G1 <- NULL
# groupgene$G4 <- NULL
# groupgene <- setNames(groupgene, paste0("G", seq_along(groupgene)))
# message("after curation: ", paste(names(groupgene), collapse = ", "),
#         " (", paste(lengths(groupgene), collapse = ", "), " genes)")

# plot_cor_heatmap(groupgene, cor_avg,
#                  file.path(outdir, "groupgene_correlation_final.png"),
#                  "Group-gene correlation (final)")

# # ── Save ─────────────────────────────────────────────────────────────────────
# saveRDS(groupgene, file.path(outdir, "groupgene.Rds"))
# gg_df <- data.frame(
#   group = rep(names(groupgene), lengths(groupgene)),
#   gene  = unlist(groupgene, use.names = FALSE)
# )
# write.csv(gg_df, file.path(outdir, "groupgene.csv"), row.names = FALSE)

# ── Visualization: groupgene module scores across all 8 samples ────────────────
# AddModuleScore each (now exclusive) group's genes per sample, then render
# FeaturePlot (UMAP) and the spatial map. These tumour_srt.qs2 are FOV-based
# Visium HD objects → SpatialFeaturePlot has no compatible image, so the spatial
# panel uses ImageFeaturePlot (the FOV equivalent used elsewhere in 6.2.1).
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
  srt_f <- file.path(base_dir, s, "tumour", "tumour_srt.qs2")
  if (!file.exists(srt_f)) { message("Skipping ", s, " — missing tumour_srt.qs2"); next }
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

  red <- grep("pearsonumap", Reductions(srt), ignore.case = TRUE, value = TRUE)[1]
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

# ── Cross-sample summary heatmaps: per-GGmod avg_log2FC + enrichment NES ────────
# From the per-sample DEG (module-gene avg_log2FC) and GSEA NES collected in the
# scoring loop above, draw one heatmap per GGmod group: genes / pathways in rows,
# samples in columns. Colour scale steelblue–white–indianred (low/0/high).
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
  # (1) avg_log2FC of this module's defining genes (top-30% module cells vs rest)
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

  # # (2) GSEA NES per collection; keep pathways significant (p.adjust < 0.05) in >= 1 sample
  # for (nm in c("Hallmark", "C6", "C5")) {
  #   tabs <- setNames(lapply(samples, function(s) {
  #     f <- file.path(ps_dir, s, sprintf("enrich_%s_%s.csv", grp, nm))
  #     if (!file.exists(f)) return(NULL)
  #     d <- read.csv(f)
  #     if (nrow(d) == 0) NULL else d
  #   }), samples)
  #   nes <- lapply(tabs, function(d) if (is.null(d)) numeric(0) else setNames(d$NES, d$ID))
  #   m <- build_mat(nes)
  #   if (is.null(m)) next
  #   sig_ids <- unique(unlist(lapply(tabs, function(d)
  #     if (is.null(d)) character(0) else d$ID[!is.na(d$p.adjust) & d$p.adjust < 0.05])))
  #   m <- m[rownames(m) %in% sig_ids, , drop = FALSE]
  #   if (nrow(m) < 2) next
  #   if (nrow(m) > 50) {                       # cap to most-recurrent / strongest pathways
  #     ord <- order(rowSums(!is.na(m)), rowMeans(abs(m), na.rm = TRUE), decreasing = TRUE)
  #     m <- m[ord[seq_len(50)], , drop = FALSE]
  #   }
  #   save_hm(m, file.path(summ_dir, sprintf("enrichNES_%s_%s.png", grp, nm)),
  #           paste0(grp, " — ", nm, " GSEA NES (sig in >= 1 sample)"), "NES")
  # }
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
hc <- row_dend(ht)
extract_group_genes <- function(dend, idx_list){
  lapply(idx_list, function(x){
    unlist(sapply(x, function(y) labels(dend[[y]])))
  })
}
groupdeg <- extract_group_genes(hc, list(c(1, 3), c(6, 7, 8, 10, 11), c(2, 4, 13, 14)))
names(groupdeg) <- c("G3", "G2", "G1")
groupdeg <- setNames(lapply(names(groupdeg), function(grp){
  gen1 <- groupdeg[[grp]]
  submat <- mat_list[[grp]]
  gen2 <- rownames(submat[complete.cases(submat), , drop = F])
  intersect(gen1, gen2)
}), names(groupdeg))

saveRDS(groupdeg, file.path(summ_dir, "groupdeg.rds"))
ht2 <- save_hm(mat[unlist(groupdeg), ], file.path(summ_dir, "avglogFC_all_groups_selected.png"),
        "Selected DEGs for group genes", "avg_log2FC", k = 3)

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
    deg_sp_plots[[m]][[s]] <- ImageFeaturePlot(srt, features = m, size = 0.4) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred", midpoint = mid) +
      ggtitle(paste0(s, "  ", m)) + theme(plot.title = element_text(size = 8))
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
         width = grid_w, height = grid_h, dpi = 150, limitsize = FALSE)

message("\nDone. groupdeg module-score figures saved to ", summ_dir)
