# =============================================================================
# TME Cluster-Level Annotation Pipeline for VisHD / SpaNorm  (v4)
# =============================================================================
#
# Key changes from v3
# ───────────────────
#  Step 2 (re-clustering) removed entirely. The function now accepts existing
#  cluster labels via the `cluster_col` parameter (any meta.data column) and
#  annotates those clusters directly.
#
#  `cluster_col` (default "seurat_clusters") – name of the meta.data column
#  that holds pre-computed cluster IDs. Any column is accepted, allowing the
#  user to supply Leiden clusters, manually merged clusters, spatial domains,
#  or any other grouping.
#
# Pipeline steps
# ──────────────
#  STEP 1 │ Soft QC  – label cells; no removal
#  STEP 2 │ Annotate clusters using the supplied cluster column
#          │   2a. pseudo-bulk cluster means
#          │   2b. z-score across clusters
#          │   2c. score with PRIMARY markers; apply exclusivity penalty
#          │   2d. compute confidence; if below threshold -> try SECONDARY
#          │   2e. record qc_low fraction per cluster
#          │   2f. propagate labels + qc_low_frac to cells
# =============================================================================

library(Seurat)
library(Matrix)
library(dplyr)

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

.row_zscore <- function(mat) {
  if (inherits(mat, "sparseMatrix")) mat <- as.matrix(mat)
  mu  <- rowMeans(mat, na.rm = TRUE)
  sig <- apply(mat, 1, sd, na.rm = TRUE)
  sig[sig == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sig, "/")
}

.sigmoid <- function(x) 1 / (1 + exp(-x))

.hdr <- function(txt, width = 72) {
  bar <- paste(rep("-", width), collapse = "")
  message("\n", bar, "\n  ", txt, "\n", bar)
}

# Score clusters x cell-types given a z-scored cluster-mean matrix and a
# filtered marker list.  Returns a list: score_mat, adj_score_mat, coverage_mat.
.score_clusters <- function(z_cluster,
                            marker_lists,
                            valid_types,
                            exclusivity_weight,
                            trim) {
  
  unique_cls <- colnames(z_cluster)
  n_cls      <- length(unique_cls)
  
  score_mat    <- matrix(0, nrow = n_cls, ncol = length(valid_types),
                         dimnames = list(unique_cls, valid_types))
  coverage_mat <- matrix(0, nrow = n_cls, ncol = length(valid_types),
                         dimnames = list(unique_cls, valid_types))
  
  for (ct in valid_types) {
    mkrs  <- marker_lists[[ct]]
    sub_z <- t(z_cluster[mkrs, , drop = FALSE])   # clusters x markers
    score_mat[, ct]    <- apply(sub_z, 1, mean, trim = trim)
    coverage_mat[, ct] <- rowMeans(sub_z > 0, na.rm = TRUE)
  }
  
  # Exclusivity adjustment
  adj_score_mat <- score_mat
  for (ci in seq_len(n_cls)) {
    row_s <- score_mat[ci, ]
    for (ct in valid_types) {
      comp_max <- max(row_s[names(row_s) != ct], na.rm = TRUE)
      adj_score_mat[ci, ct] <- row_s[ct] - exclusivity_weight * max(comp_max, 0)
    }
  }
  
  list(score     = score_mat,
       adj_score = adj_score_mat,
       coverage  = coverage_mat)
}

# Given adj_score_mat + coverage_mat, return per-cluster best label, runner-up,
# and three-component confidence.  Returns a named list of vectors (not a
# data.frame) so that cluster-ID subsetting with [cl] always works reliably.
# Robust to: single cell type, NaN/Inf scores, ties, all-zero scores.
.pick_best <- function(adj_score_mat, coverage_mat) {
  n_cls       <- nrow(adj_score_mat)
  valid_types <- colnames(adj_score_mat)
  n_types     <- length(valid_types)
  cls_names   <- rownames(adj_score_mat)   # cluster IDs as character
  
  # Replace any NA/NaN in scores with a large negative so comparisons are safe
  adj_score_mat[!is.finite(adj_score_mat)] <- -1e6
  
  best_idx   <- apply(adj_score_mat, 1, which.max)
  best_label <- valid_types[best_idx]
  best_score <- adj_score_mat[cbind(seq_len(n_cls), best_idx)]
  
  # Runner-up: only meaningful when >= 2 types exist
  if (n_types >= 2) {
    runner_idx <- apply(adj_score_mat, 1, function(x) {
      x[which.max(x)] <- -1e6
      which.max(x)
    })
    runner_label <- valid_types[runner_idx]
    runner_score <- adj_score_mat[cbind(seq_len(n_cls), runner_idx)]
  } else {
    runner_label <- best_label
    runner_score <- rep(0, n_cls)
  }
  
  strength <- .sigmoid(best_score)
  coverage <- coverage_mat[cbind(seq_len(n_cls), best_idx)]
  # coverage can be NA if only one marker — replace with 0
  coverage[!is.finite(coverage)] <- 0
  
  # Margin: gap between best and runner-up, normalised to [0, 1]
  margin_raw <- best_score - runner_score
  margin_raw[!is.finite(margin_raw)] <- 0
  denom         <- pmax(abs(best_score), 1e-3)           # avoid near-zero denom
  margin_scaled <- pmin(pmax(margin_raw / denom, 0), 1)  # clamp to [0, 1]
  if (n_types == 1) margin_scaled <- rep(1, n_cls)       # no competition
  
  confidence <- (strength + coverage + margin_scaled) / 3
  confidence[!is.finite(confidence)] <- 0                # final safety net
  
  # Return as named vectors keyed by cluster ID — safe for [cl] subsetting
  list(
    best_label   = setNames(best_label,          cls_names),
    runner_label = setNames(runner_label,         cls_names),
    best_score   = setNames(round(best_score, 4), cls_names),
    confidence   = setNames(round(confidence, 4), cls_names)
  )
}


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

#' Soft-QC + two-tier cluster annotation for VisHD SpaNorm objects
#'
#' @param obj                 Seurat object with a SpaNorm assay.
#' @param tme_markers         Named list of PRIMARY marker gene vectors.
#' @param secondary_genes     Either:
#'   \itemize{
#'     \item A plain character vector – used ONLY for QC (broader panel);
#'           no fallback annotation is performed.
#'     \item A named list of gene vectors (same format as tme_markers) –
#'           used for BOTH QC and fallback annotation when primary confidence
#'           falls below conf_threshold.
#'   }
#'   Pass NULL (default) to skip secondary QC and fallback entirely.
#' @param cluster_col         Name of the meta.data column containing
#'                            pre-computed cluster IDs to annotate.
#'                            Accepts any grouping: seurat_clusters, Leiden
#'                            clusters, spatial domains, manual merges, etc.
#'                            (default "seurat_clusters").
#' @param assay               Assay name (default "SpaNorm").
#' @param data_slot           Slot: "data" or "scale.data" (default "data").
#'
#' QC parameters (labelling only -- no cells are removed)
#' @param expr_min_val        Threshold above which a gene counts as expressed
#'                            (default 0 = any non-zero value).
#' @param primary_expr_frac   Fraction of PRIMARY panel genes a cell must
#'                            express to pass primary QC (default 0.10).
#' @param secondary_expr_frac Fraction of SECONDARY panel a cell must express
#'                            to pass secondary QC (default 0.05).
#'
#' Annotation parameters
#' @param min_markers         Minimum usable markers to score a cell type
#'                            (default 2).
#' @param conf_threshold      Clusters below this confidence trigger fallback
#'                            to secondary; if still below -> "Unknown"
#'                            (default 0.15).
#' @param exclusivity_weight  Competitor-score penalty weight (default 0.30).
#' @param detection_min       Minimum per-gene detection rate to use a marker
#'                            (default 0.01).
#' @param trim                Trim fraction for mean across markers (0.10).
#'
#' @return Annotated Seurat object (ALL cells retained) with new meta.data:
#'   primary_expr_frac    fraction of PRIMARY panel expressed per cell
#'   secondary_expr_frac  fraction of SECONDARY panel expressed (NA if unused)
#'   qc_label             "qc_pass" | "qc_low_primary" |
#'                        "qc_low_secondary" | "qc_low_both"
#'   celltype_annotation  best cell type label
#'   celltype_confidence  confidence score [0, 1]
#'   celltype_score_raw   winning adjusted score
#'   celltype_runner_up   second-best label
#'   annotation_source    "primary" | "secondary" | "unknown"
#'   cluster_qc_low_frac  fraction of qc_low_* cells in the cluster

tme_cluster_annotation_pipeline <- function(
    obj,
    tme_markers,
    secondary_genes      = NULL,
    cluster_col          = "seurat_clusters",
    assay                = "SpaNorm",
    data_slot            = "data",
    # QC
    expr_min_val         = 0,
    primary_expr_frac    = 0.10,
    secondary_expr_frac  = 0.05,
    # Annotation
    min_markers          = 2,
    conf_threshold       = 0.15,
    exclusivity_weight   = 0.30,
    detection_min        = 0.01,
    trim                 = 0.10
) {
  
  # ============================================================================
  # 0. Validate & parse inputs
  # ============================================================================
  stopifnot(
    inherits(obj, "Seurat"),
    assay %in% names(obj@assays),
    is.list(tme_markers), !is.null(names(tme_markers)),
    is.character(cluster_col), length(cluster_col) == 1,
    cluster_col %in% colnames(obj@meta.data)
  )
  
  if (!cluster_col %in% colnames(obj@meta.data))
    stop("[pipeline] cluster_col '", cluster_col,
         "' not found in meta.data. Available columns: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  
  # Determine mode of secondary_genes
  has_secondary      <- !is.null(secondary_genes)
  secondary_is_named <- has_secondary &&
    is.list(secondary_genes) &&
    !is.null(names(secondary_genes))
  
  secondary_panel_vec <- if (has_secondary) unique(unlist(secondary_genes)) else character(0)
  
  DefaultAssay(obj) <- assay
  n_cells <- ncol(obj)
  message("[pipeline] Input: ", n_cells, " cells | assay: ", assay,
          " | slot: ", data_slot, " | cluster_col: '", cluster_col, "'")
  
  # ============================================================================
  # STEP 1 -- Soft QC (label only, retain all cells)
  # ============================================================================
  .hdr("STEP 1 | Soft QC -- labelling cells (no removal)")
  
  expr_mat  <- GetAssayData(obj, assay = assay, layer = data_slot)
  all_genes <- rownames(expr_mat)
  
  # -- Primary panel ----------------------------------------------------------
  primary_panel <- intersect(unique(unlist(tme_markers)), all_genes)
  n_primary     <- length(primary_panel)
  if (n_primary == 0) stop("[pipeline] No primary marker genes found in assay.")
  message("[pipeline] Primary panel  : ", length(unique(unlist(tme_markers))),
          " listed -> ", n_primary, " found in assay")
  
  pm_mat       <- expr_mat[primary_panel, , drop = FALSE]
  pm_bin       <- if (inherits(pm_mat, "sparseMatrix"))
    Matrix::colMeans(pm_mat > expr_min_val)
  else
    colMeans(pm_mat > expr_min_val, na.rm = TRUE)
  primary_frac <- pm_bin                       # named vector: cells
  pass_primary <- primary_frac >= primary_expr_frac
  
  # -- Secondary panel --------------------------------------------------------
  if (has_secondary) {
    sec_panel <- intersect(secondary_panel_vec, all_genes)
    n_sec     <- length(sec_panel)
    message("[pipeline] Secondary panel: ", length(secondary_panel_vec),
            " listed -> ", n_sec, " found in assay")
    
    if (n_sec == 0) {
      warning("[pipeline] No secondary genes found in assay; secondary QC skipped.")
      has_secondary      <- FALSE
      secondary_is_named <- FALSE
      secondary_frac     <- setNames(rep(NA_real_, n_cells), colnames(obj))
      pass_secondary     <- setNames(rep(TRUE,     n_cells), colnames(obj))
    } else {
      sm_mat         <- expr_mat[sec_panel, , drop = FALSE]
      sm_bin         <- if (inherits(sm_mat, "sparseMatrix"))
        Matrix::colMeans(sm_mat > expr_min_val)
      else
        colMeans(sm_mat > expr_min_val, na.rm = TRUE)
      secondary_frac <- sm_bin
      pass_secondary <- secondary_frac >= secondary_expr_frac
    }
  } else {
    secondary_frac <- setNames(rep(NA_real_, n_cells), colnames(obj))
    pass_secondary <- setNames(rep(TRUE,     n_cells), colnames(obj))
  }
  
  # -- Assign QC labels -------------------------------------------------------
  pp <- pass_primary[colnames(obj)]
  ps <- pass_secondary[colnames(obj)]
  
  qc_label <- dplyr::case_when(
    pp &  ps ~ "qc_pass",
    !pp &  ps ~ "qc_low_primary",
    pp & !ps ~ "qc_low_secondary",
    !pp & !ps ~ "qc_low_both",
    TRUE      ~ "qc_pass"
  )
  
  obj$primary_expr_frac   <- round(primary_frac[colnames(obj)],   4)
  obj$secondary_expr_frac <- round(secondary_frac[colnames(obj)],  4)
  obj$qc_label            <- qc_label
  
  # -- QC summary -------------------------------------------------------------
  tab_qc <- sort(table(qc_label), decreasing = TRUE)
  message("[pipeline] QC label distribution (ALL cells retained):")
  print(tab_qc)
  message("[pipeline] Primary   threshold: >= ",
          round(100 * primary_expr_frac, 1), "% of ", n_primary, " genes")
  if (has_secondary)
    message("[pipeline] Secondary threshold: >= ",
            round(100 * secondary_expr_frac, 1), "% of ", n_sec, " genes")
  
  obj@misc$qc_stats <- list(
    n_cells               = n_cells,
    primary_panel_genes   = primary_panel,
    secondary_panel_genes = if (has_secondary) sec_panel else NULL,
    primary_expr_frac     = primary_expr_frac,
    secondary_expr_frac   = if (has_secondary) secondary_expr_frac else NULL,
    qc_label_table        = tab_qc
  )
  
  # ============================================================================
  # STEP 2 -- Cluster annotation using supplied cluster_col
  # ============================================================================
  .hdr("STEP 2 | Cluster annotation  (primary markers; fallback to secondary)")
  
  expr_qc     <- GetAssayData(obj, assay = assay, layer = data_slot)
  cluster_ids <- as.character(obj@meta.data[[cluster_col]])
  unique_cls  <- as.character(sort(unique(cluster_ids)))
  n_clusters  <- length(unique_cls)
  
  message("[pipeline] Annotating ", n_clusters, " clusters from '",
          cluster_col, "'")
  print(table(cluster_ids))
  
  # -- 2a. Pseudo-bulk cluster means -----------------------------------------
  message("[pipeline] Computing pseudo-bulk means ...")
  cluster_means <- sapply(unique_cls, function(cl) {
    idx <- which(cluster_ids == cl)
    if (length(idx) == 1) as.numeric(expr_qc[, idx])
    else Matrix::rowMeans(expr_qc[, idx, drop = FALSE], na.rm = TRUE)
  })
  rownames(cluster_means) <- rownames(expr_qc)   # genes x clusters
  
  # -- 2b. Z-score cluster means across clusters -----------------------------
  message("[pipeline] Z-scoring cluster means ...")
  z_cluster <- .row_zscore(cluster_means)
  
  # -- 2c. Per-gene detection: use cluster-level mean expression -------------
  # Using global single-cell detection rate is too stringent for sparse VisHD
  # data. Instead we ask: is the gene detectably expressed in at least one
  # cluster (i.e. cluster mean > detection_min)?  This retains rare-cell-type
  # markers that are lowly detected globally but enriched in one cluster.
  # cluster_means is genes x clusters (dense), already computed above.
  gene_max_cluster_mean <- rowMaxs_safe <- apply(cluster_means, 1, max,
                                                 na.rm = TRUE)
  # A gene is "usable" if its maximum cluster mean exceeds detection_min
  gene_usable <- gene_max_cluster_mean > detection_min
  
  # Helper: filter one marker list to genes that are usable
  .filter_markers <- function(mkr_list) {
    lapply(mkr_list, function(genes) {
      found <- intersect(genes, rownames(expr_qc))
      found <- found[gene_usable[found]]
      if (length(found) < min_markers) return(NULL)
      found
    })
  }
  
  # -- 2d. Filter primary markers --------------------------------------------
  primary_filtered <- .filter_markers(tme_markers)
  valid_primary    <- names(Filter(Negate(is.null), primary_filtered))
  
  dropped_p <- setdiff(names(tme_markers), valid_primary)
  if (length(dropped_p) > 0)
    message("[pipeline] Primary types with insufficient markers: ",
            paste(dropped_p, collapse = ", "))
  message("[pipeline] Primary scoring: ", length(valid_primary), " types: ",
          paste(valid_primary, collapse = ", "))
  
  if (length(valid_primary) == 0)
    stop("[pipeline] No primary cell types have >= ", min_markers,
         " usable markers. Lower detection_min or min_markers.")
  
  # -- 2e. Filter secondary markers (named list only) ------------------------
  has_secondary_annot <- secondary_is_named
  valid_secondary     <- character(0)
  res_secondary       <- NULL
  
  if (has_secondary_annot) {
    secondary_filtered <- .filter_markers(secondary_genes)
    valid_secondary    <- names(Filter(Negate(is.null), secondary_filtered))
    dropped_s <- setdiff(names(secondary_genes), valid_secondary)
    if (length(dropped_s) > 0)
      message("[pipeline] Secondary types with insufficient markers: ",
              paste(dropped_s, collapse = ", "))
    message("[pipeline] Secondary scoring: ", length(valid_secondary),
            " types: ", paste(valid_secondary, collapse = ", "))
  }
  
  # -- 2f. Score with primary markers ----------------------------------------
  res_primary  <- .score_clusters(z_cluster, primary_filtered,
                                  valid_primary, exclusivity_weight, trim)
  best_primary <- .pick_best(res_primary$adj_score, res_primary$coverage)
  # best_primary rows = clusters
  
  # -- 2g. Score with secondary markers (pre-compute once if available) ------
  if (has_secondary_annot && length(valid_secondary) > 0) {
    res_secondary  <- .score_clusters(z_cluster, secondary_filtered,
                                      valid_secondary, exclusivity_weight, trim)
    best_secondary <- .pick_best(res_secondary$adj_score, res_secondary$coverage)
  }
  
  # -- 2h. Per-cluster: apply fallback logic ---------------------------------
  # best_primary and best_secondary are named lists of vectors, keyed by
  # cluster ID string — safe for direct [cl] subsetting.
  annotation_source <- setNames(rep("primary", n_clusters), unique_cls)
  final_label       <- best_primary$best_label[unique_cls]
  final_confidence  <- best_primary$confidence[unique_cls]
  final_score       <- best_primary$best_score[unique_cls]
  final_runner      <- best_primary$runner_label[unique_cls]
  
  # Diagnostic: show raw primary confidence per cluster
  message("[pipeline] Primary confidence per cluster:")
  print(round(final_confidence, 3))
  
  for (cl in unique_cls) {
    prim_conf <- final_confidence[[cl]]
    # Treat non-finite confidence as below threshold
    if (is.finite(prim_conf) && prim_conf >= conf_threshold) next
    
    if (has_secondary_annot && length(valid_secondary) > 0) {
      sec_conf <- best_secondary$confidence[[cl]]
      if (is.finite(sec_conf) && sec_conf >= conf_threshold) {
        final_label[[cl]]      <- best_secondary$best_label[[cl]]
        final_confidence[[cl]] <- sec_conf
        final_score[[cl]]      <- best_secondary$best_score[[cl]]
        final_runner[[cl]]     <- best_secondary$runner_label[[cl]]
        annotation_source[[cl]]<- "secondary"
      } else {
        final_label[[cl]]      <- "Unknown"
        final_confidence[[cl]] <- max(
          c(prim_conf, sec_conf)[is.finite(c(prim_conf, sec_conf))],
          0  # fallback if both non-finite
        )
        annotation_source[[cl]]<- "unknown"
      }
    } else {
      final_label[[cl]]      <- "Unknown"
      final_confidence[[cl]] <- if (is.finite(prim_conf)) prim_conf else 0
      annotation_source[[cl]]<- "unknown"
    }
  }
  
  # -- 2i. QC-low fraction per cluster ---------------------------------------
  is_qc_low       <- obj$qc_label != "qc_pass"
  qc_low_frac     <- sapply(unique_cls, function(cl) {
    idx <- which(cluster_ids == cl)
    round(mean(is_qc_low[idx], na.rm = TRUE), 4)
  })
  
  # -- 2j. Build cluster annotation table ------------------------------------
  cluster_annot <- data.frame(
    cluster             = unique_cls,
    celltype_annotation = final_label[unique_cls],
    celltype_confidence = round(final_confidence[unique_cls], 4),
    celltype_score_raw  = round(final_score[unique_cls], 4),
    celltype_runner_up  = final_runner[unique_cls],
    annotation_source   = annotation_source[unique_cls],
    cluster_qc_low_frac = qc_low_frac[unique_cls],
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )
  
  message("\n[pipeline] Cluster annotation table:")
  print(cluster_annot[, c("cluster", "celltype_annotation",
                          "celltype_confidence", "annotation_source",
                          "cluster_qc_low_frac")])
  
  n_unknown <- sum(cluster_annot$celltype_annotation == "Unknown")
  message("[pipeline] Unknown clusters: ", n_unknown, " / ", n_clusters)
  
  # -- 2k. Propagate cluster labels to every cell ----------------------------
  .mk_map <- function(col) setNames(cluster_annot[[col]], cluster_annot$cluster)
  
  obj$celltype_annotation  <- as.character(.mk_map("celltype_annotation")[cluster_ids])
  obj$celltype_confidence  <- as.character(.mk_map("celltype_confidence")[cluster_ids])
  obj$celltype_score_raw   <- as.character(.mk_map("celltype_score_raw")[cluster_ids])
  obj$celltype_runner_up   <- as.character(.mk_map("celltype_runner_up")[cluster_ids])
  obj$annotation_source    <- as.character(.mk_map("annotation_source")[cluster_ids])
  obj$cluster_qc_low_frac  <- as.character(.mk_map("cluster_qc_low_frac")[cluster_ids])
  
  # -- 2l. Store full score tables in @misc ----------------------------------
  obj@misc$cluster_annotation      <- cluster_annot
  obj@misc$cluster_means           <- as.data.frame(cluster_means)
  obj@misc$primary_score_raw       <- as.data.frame(res_primary$score)
  obj@misc$primary_score_adj       <- as.data.frame(res_primary$adj_score)
  obj@misc$primary_coverage        <- as.data.frame(res_primary$coverage)
  obj@misc$primary_markers_used    <- primary_filtered[valid_primary]
  
  if (!is.null(res_secondary)) {
    obj@misc$secondary_score_raw    <- as.data.frame(res_secondary$score)
    obj@misc$secondary_score_adj    <- as.data.frame(res_secondary$adj_score)
    obj@misc$secondary_coverage     <- as.data.frame(res_secondary$coverage)
    obj@misc$secondary_markers_used <- secondary_filtered[valid_secondary]
  }
  
  # ============================================================================
  # Final summary
  # ============================================================================
  .hdr("DONE")
  message("[pipeline] Cell-level annotation summary:")
  print(sort(table(obj$celltype_annotation), decreasing = TRUE))
  message("\n[pipeline] Annotation source (cells):")
  print(table(obj$annotation_source))
  message("\n[pipeline] Mean confidence across all cells: ",
          round(mean(obj$celltype_confidence, na.rm = TRUE), 3))
  message("[pipeline] Cells in clusters with >20% qc_low: ",
          sum(obj$cluster_qc_low_frac > 0.20, na.rm = TRUE))
  
  return(obj)
}



#' Dot plot: marker expression per annotated cell type
#'
#' @param obj           Annotated Seurat object.
#' @param tme_markers   Primary or secondary marker list.
#' @param max_per_type  Max markers to show per cell type (default 4).
plot_marker_dotplot <- function(obj, tme_markers, max_per_type = 4) {
  features <- unique(unlist(lapply(tme_markers, function(g)
    head(intersect(g, rownames(obj)), max_per_type))))
  
  DotPlot(obj, features = features,
          group.by = "celltype_annotation",
          assay    = DefaultAssay(obj)) +
    ggplot2::coord_flip() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ggtitle("Marker Expression per Annotated Cell Type")
}


# =============================================================================
# Example usage
# =============================================================================

# tme_markers <- list(
#   Epithelial  = c("EPCAM","KRT8","KRT18","CDH1","MUC1"),
#   T_Cells_Pan = c("CD3D","CD3E","CD3G","CD2"),
#   CD8_T       = c("CD8A","CD8B","GZMK","GZMB"),
#   CD4_T       = c("CD4","IL7R","LTB"),
#   Tregs       = c("FOXP3","IL2RA","IKZF2","BATF"),
#   B_Cells     = c("CD19","MS4A1","CD79A"),
#   Plasma      = c("JCHAIN","MZB1","SDC1"),
#   NK_Cells    = c("NCAM1","KLRB1","NKG7","GNLY"),
#   Myeloid_Pan = c("CD14","LYZ","CD68","CSF1R"),
#   TAMs        = c("CD163","MRC1","APOE","C1QA"),
#   DCs         = c("HLA-DRA","CD1C","CLEC9A","THBD"),
#   CAFs        = c("COL1A1","DCN","LUM","FAP","ACTA2"),
#   Endothelial = c("PECAM1","VWF","ENG","CLDN5"),
#   Pericytes   = c("RGS5","MCAM","CSPG4")
# )
#
# # ── Option A: secondary = plain vector (QC only, no fallback annotation) ──
# secondary_qc_panel <- c(
#   "ACTB","GAPDH","B2M","HSP90AB1","RPL13","RPS18","MALAT1","NEAT1",
#   "EEF1A1","TPT1","FTL","FTH1","HSPA1A","HSPA1B"
# )
#
# # ── Option B: secondary = named list (QC + fallback annotation) ────────────
# secondary_markers <- list(
#   Epithelial_broad = c("EPCAM","KRT5","KRT14","KRT17","KRT19","CLDN4","OCLN",
#                        "TJP1","GJB2","GRHL2","ESRP1","ELF3"),
#   Immune_broad     = c("PTPRC","HLA-A","HLA-B","HLA-C","B2M","LAPTM5",
#                        "TYROBP","FCER1G","LST1","AIF1"),
#   Stromal_broad    = c("VIM","FN1","THY1","S100A4","COL3A1","PDGFRA",
#                        "PDGFRB","SPARC","POSTN","IGFBP3","IGFBP4"),
#   Myeloid_broad    = c("CD14","ITGAM","ITGAX","CX3CR1","CSF1R","FCGR3A",
#                        "S100A8","S100A9","VCAN","FCN1","LILRB2"),
#   Lymphoid_broad   = c("CD3D","CD3E","CD19","MS4A1","NCAM1","IL2RA",
#                        "CD27","CD38","SELL","CCR7","IL7R")
# )
#
# # ── Run the pipeline ───────────────────────────────────────────────────────
# # Use existing seurat_clusters (default)
# visHD_annotated <- tme_cluster_annotation_pipeline(
#   obj                 = visHD_obj,
#   tme_markers         = tme_markers,
#   secondary_genes     = secondary_markers,
#   cluster_col         = "seurat_clusters",   # or e.g. "leiden_res0.5", "spatial_domain"
#   assay               = "SpaNorm",
#   data_slot           = "data",
#   expr_min_val        = 0,
#   primary_expr_frac   = 0.10,
#   secondary_expr_frac = 0.05,
#   min_markers         = 2,
#   conf_threshold      = 0.15,
#   exclusivity_weight  = 0.30,
#   detection_min       = 0.01,
#   trim                = 0.10
# )
#
# # Cluster table with qc_low_frac and annotation_source
# visHD_annotated@misc$cluster_annotation
#
# # Overview plot (6 panels)
# plot_tme_pipeline(visHD_annotated, cluster_col = "seurat_clusters")
#
# # Marker dot plot
# plot_marker_dotplot(visHD_annotated, tme_markers)
#
# # Which clusters fell back to secondary?
# subset(visHD_annotated@misc$cluster_annotation,
#        annotation_source %in% c("secondary", "unknown"))