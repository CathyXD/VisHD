spatial_plot <- function(srt, outdir, name) {
  require(Seurat)
  require(ggplot2)
  require(pals)
  require(patchwork)
  
  tryCatch(
    {
      f1 <- ImageDimPlot(srt, cols = "polychrome",border.color = "#00000000", size = 0.3)
      f2 <- ImageDimPlot(srt,  split.by = "seurat_clusters", cols = "polychrome",border.color = "#00000000", size = 0.3)
      ggsave(plot = f1 + f2 + plot_layout(ncol = 2), paste0(outdir, name, "_ImageDimPlot.png"), width = 10, height = 5, dpi = 350, create.dir =T)
      
      f1 <- ImageDimPlot(srt, group.by = "SR_Cluster", cols = "polychrome",border.color = "#00000000", size = 0.3)
      f2 <- ImageDimPlot(srt, group.by = "SR_Cluster",  split.by = "SR_Cluster", cols = "polychrome",border.color = "#00000000", size = 0.3)
      ggsave(plot = f1 + f2 + plot_layout(ncol = 2), paste0(outdir, name, "SR_Cluster_ImageDimPlot.png"), width = 10, height = 5, dpi = 350, create.dir =T)
      
      
      f3 <- ImageFeaturePlot(srt, c("AR", "FOLH1", "ASCL1", "COL1A1","cell_area", "nucleus_area","nFeature_Spatial", "nCount_Spatial"), size = 0.3)
      ggsave(plot = f3 +  plot_layout(ncol = 4), paste0(outdir, name, "_ImageFeaturePlot.png"), width = 13, height = 6, dpi = 350, create.dir =T)
      
      g1 <- DimPlot(srt, cols = "polychrome", label = T)
      g2 <- DimPlot(srt,  group.by = "SR_Cluster", cols = "polychrome", label = T)
      ggsave(plot = g1 + g2 + plot_layout(ncol = 2), paste0(outdir, name, "_DimPlot.png"), width = 10, height = 5, dpi = 350, create.dir =T)
      
      g3 <- FeaturePlot(srt, c("AR", "FOLH1", "ASCL1", "COL1A1","cell_area", "nucleus_area","nFeature_Spatial", "nCount_Spatial"))
      ggsave(plot = g3 + plot_layout(ncol = 4), paste0(outdir, name, "_FeaturePlot.png"), width = 12, height = 5, dpi = 350)
      
      p1 <- VlnPlot(srt, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, pt.size = 0, group.by = "seurat_clusters")
      p2 <- VlnPlot(srt, features = c("AR", "FOLH1", "ASCL1", "COL1A1"), ncol = 4, pt.size = 0, group.by = "seurat_clusters")        # the function you want to run first
      ggsave(plot = (p1|p2) + plot_layout(width = c(1, 2)), paste0(outdir, name, "_VlnPlot.png"), width = 15, height = 4, dpi = 350, create.dir =T)
    },
    error = function(e) {
      message("Main function failed: ", e$message)
    }
  )
}

createSPEObject <- function(countMat, 
                            spatial_locs, 
                            cell_metadata,
                            normalize = TRUE) {
  
  require(SpatialExperiment)
  require(S4Vectors)
  require(scran)
  
  rd = S4Vectors::DataFrame(symbol = rownames(countMat))
  
  data_spe = SpatialExperiment(
    assays = list(counts = countMat),
    colData = S4Vectors::DataFrame(cell_metadata),
    rowData = rd,
    spatialCoords = as.matrix(spatial_locs)
  )
  
  if (normalize == TRUE) {
    message("performing normalization...............")
    set.seed(123)
    qclus <- quickCluster(data_spe)
    data_spe <- computeSumFactors(data_spe, cluster = qclus)
    data_spe$sizeFactor <- pmax(1e-08,data_spe$sizeFactor)
    data_spe <- scater::logNormCounts(data_spe)
  }
  
  return(data_spe)
}

do.spanorm <- function(srt, outdir = "."){
  require(SpaNorm, lib.loc = "~/R_Library/4.5")
  require(leidenbase,  lib.loc = "~/R_Library/4.5")
  # Seurat v5 keeps per-sample layers after merge — collapse so GetAssayData works
  if (length(SeuratObject::Layers(srt, assay = "Spatial")) > 1) {
    srt <- SeuratObject::JoinLayers(srt, assay = "Spatial")
  }
  countMat <- GetAssayData(srt, assay = "Spatial", layer = "counts")

  # After cross-sample merge there is no image — fall back to centroids stashed in meta.data
  spatial_coord <- tryCatch(
    GetTissueCoordinates(srt)[, c("x", "y")],
    error = function(e) {
      if (all(c("x_centroid", "y_centroid") %in% colnames(srt@meta.data))) {
        setNames(srt@meta.data[, c("x_centroid", "y_centroid")], c("x", "y"))
      } else {
        stop("do.spanorm: no image and no x_centroid/y_centroid in meta.data")
      }
    }
  )
  rownames(spatial_coord) <- colnames(countMat)
  set.seed(1)
  spe <- createSPEObject(countMat, spatial_locs= spatial_coord, cell_metadata = spatial_coord, normalize = FALSE)
  spe <- SpaNorm(spe, df.tps = 4, sample.p = 0.25, gene.model = "nb")
  spe <- SpaNormSVG(spe)
  saveRDS(rowData(spe), file.path(outdir, "SVGs.Rds"))
  spe <- SpaNormPCA(spe, ncomponents = 30, svg.fdr = 0.2, nsvgs = Inf)
  pca_embeddings <- reducedDim(spe, "PCA")
  rownames(pca_embeddings) <- colnames(srt)
  topsvg <- topSVGs(spe, n = Inf, fdr = 0.2)
  
  srt[["SpaNorm"]] <- CreateAssayObject(data =  assay(spe, "logcounts"))
  DefaultAssay(srt) <- "SpaNorm"
  
  VariableFeatures(srt) <- rownames(topsvg)
  srt[["pca"]] <- CreateDimReducObject(
    embeddings = pca_embeddings,
    key = "PC_",              # Use "UMAP_" for UMAP, "tSNE_" for tSNE
    assay = DefaultAssay(srt) # Likely "Spatial" or "RNA"
  )
  srt <- srt %>% RunUMAP(dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% 
    FindClusters(resolution = 1,algorithm = 4)
}

do.pearson_pca <- function(srt,
                           batch_variable   = NULL,
                           nfeatures        = 5000,
                           assay            = "Spatial",
                           find_hvgs        = NULL,
                           reduction_prefix = NULL,
                           clusters_col     = NULL,
                           resolution       = 0.8) {
  require(scPearsonPCA, lib.loc = "~/R_Library/4.5")
  require(Seurat)
  require(data.table)

  batched <- !is.null(batch_variable)
  if (is.null(reduction_prefix)) {
    reduction_prefix <- if (batched) "pearsonbatch" else "pearson"
  }
  if (is.null(clusters_col)) {
    clusters_col <- if (batched) "pearson_clusters_batch" else "pearson_clusters"
  }

  counts_mat <- GetAssayData(srt, assay = assay, layer = "counts")
  tc         <- Matrix::colSums(counts_mat)

  # Auto-skip FindVariableFeatures if VFs already exist in `assay`; explicit
  # find_hvgs = TRUE/FALSE overrides that decision.
  existing_hvgs <- tryCatch(VariableFeatures(srt, assay = assay),
                            error = function(e) character(0))
  if (is.null(find_hvgs)) find_hvgs <- length(existing_hvgs) == 0
  if (find_hvgs) {
    srt  <- FindVariableFeatures(srt, nfeatures = nfeatures, assay = assay)
    hvgs <- VariableFeatures(srt, assay = assay)
  } else {
    if (length(existing_hvgs) == 0) {
      stop(sprintf("find_hvgs = FALSE but no VariableFeatures set on assay '%s'", assay))
    }
    hvgs <- existing_hvgs
  }

  # scPearsonPCA rejects zero-variance genes (check_x_is_dgcmatrix); drop them.
  hvg_var <- Matrix::rowMeans(counts_mat[hvgs, ]^2) - Matrix::rowMeans(counts_mat[hvgs, ])^2
  if (any(hvg_var <= 0)) {
    message(sprintf("Dropping %d zero-variance HVGs before Pearson PCA", sum(hvg_var <= 0)))
    hvgs <- hvgs[hvg_var > 0]
  }

  if (batched) {
    if (!batch_variable %in% colnames(srt@meta.data)) {
      stop(sprintf("batch_variable '%s' not found in srt@meta.data", batch_variable))
    }
    srt@meta.data$cell_ID <- rownames(srt@meta.data)
    obs <- data.table(srt@meta.data)[, c("cell_ID", batch_variable), with = FALSE]

    pcaobj <- sparse_quasipoisson_pca_seurat_batch(
      counts_mat[hvgs, ],
      totalcounts    = tc,
      grate          = gene_frequency(counts_mat, obs = obs,
                                      cellid_colname = "cell_ID",
                                      batch_variable = batch_variable)[hvgs, ],
      obs            = obs,
      batch_variable = batch_variable,
      cellid_colname = "cell_ID",
      scale.max      = 10, do.scale = TRUE, do.center = TRUE
    )
  } else {
    pcaobj <- sparse_quasipoisson_pca_seurat(
      counts_mat[hvgs, ],
      totalcounts = tc,
      grate       = gene_frequency(counts_mat)[hvgs],
      scale.max   = 10, do.scale = TRUE, do.center = TRUE
    )
  }

  umapobj <- make_umap(pcaobj)

  srt[[paste0(reduction_prefix, "pca")]]   <- pcaobj$reduction.data
  srt[[paste0(reduction_prefix, "umap")]]  <- umapobj$ump
  srt[[paste0(reduction_prefix, "graph")]] <- Seurat::as.Graph(umapobj$grph)

  srt <- FindClusters(srt, graph = paste0(reduction_prefix, "graph"),
                      resolution = resolution)
  srt@meta.data[[clusters_col]] <- srt@meta.data$seurat_clusters

  srt
}

do.banksy <- function(srt){
  require(Banksy)
  require(SummarizedExperiment)
  require(SpatialExperiment)
  require(scuttle)
  require(scater)
  require(cowplot)
  require(ggplot2)
  countMat <- GetAssayData(srt, assay = "SpaNorm", layer = "data")
  spatial_coord <- GetTissueCoordinates(srt)[, c("x", "y")]
  se <- SpatialExperiment(assay = list(data = countMat), spatialCoords = as.matrix(spatial_coord))
  lambda <- c(0.2, 0.8)
  k_geom <- c(15, 30)
  se <- Banksy::computeBanksy(se, assay_name = "data", compute_agf = TRUE, k_geom = 15)
  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)
  se <- Banksy::connectClusters(se)
}


transfer_visiumhd_to_cells <- function(
    srt_cell,
    srt_bin,
    annotation_cols = NULL,
    cells_fov       = NULL,
    bins_fov        = NULL,
    verbose         = TRUE
) {
  
  # ── 0. Package checks ───────────────────────────────────────────────────────
  for (pkg in c("Seurat", "SeuratObject", "FNN")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg))
    }
  }
  
  # ── 1. Resolve FOV names ────────────────────────────────────────────────────
  .get_fov <- function(obj, fov_arg, label) {
    available <- names(obj@images)
    if (is.null(available) || length(available) == 0) {
      stop(sprintf("No FOV/image slots found in %s.", label))
    }
    if (is.null(fov_arg)) {
      fov_arg <- available[1]
      if (verbose) message(sprintf("[%s] Using FOV: '%s'", label, fov_arg))
    } else if (!fov_arg %in% available) {
      stop(sprintf(
        "FOV '%s' not found in %s. Available: %s",
        fov_arg, label, paste(available, collapse = ", ")
      ))
    }
    fov_arg
  }
  
  cells_fov <- .get_fov(srt_cell, cells_fov, "srt_cell")
  bins_fov  <- .get_fov(srt_bin,  bins_fov,  "srt_bin")
  
  # ── 2. Extract centroids ────────────────────────────────────────────────────
  .extract_centroids <- function(obj, fov, label) {
    coords <- tryCatch(
      Seurat::GetTissueCoordinates(obj, image = fov, which = "centroids"),
      error = function(e) {
        # Fallback: some Seurat versions use different slot structures
        img <- obj@images[[fov]]
        if (inherits(img, "FOV")) {
          SeuratObject::GetTissueCoordinates(img, which = "centroids")
        } else {
          stop(sprintf(
            "Could not extract centroids from %s FOV '%s'. Error: %s",
            label, fov, conditionMessage(e)
          ))
        }
      }
    )
    
    # Normalise column names — Seurat versions differ (x/y vs imagerow/imagecol)
    col_map <- list(
      x = c("x", "imagerow", "row"),
      y = c("y", "imagecol", "col")
    )
    for (target in names(col_map)) {
      candidates <- col_map[[target]]
      hit <- intersect(candidates, tolower(names(coords)))
      if (length(hit) > 0 && !target %in% names(coords)) {
        names(coords)[tolower(names(coords)) == hit[1]] <- target
      }
    }
    
    if (!all(c("x", "y") %in% names(coords))) {
      stop(sprintf(
        "Cannot find x/y coordinate columns in %s centroids. Found: %s",
        label, paste(names(coords), collapse = ", ")
      ))
    }
    
    coords
  }
  
  if (verbose) message("Extracting centroids...")
  cell_coords <- .extract_centroids(srt_cell, cells_fov, "srt_cell")
  bin_coords  <- .extract_centroids(srt_bin,  bins_fov,  "srt_bin")
  
  # Barcode column — Seurat stores it as "cell" or uses rownames
  .get_ids <- function(coords, obj_barcodes) {
    id_col <- intersect(c("cell", "barcode", "id"), names(coords))
    if (length(id_col) > 0) {
      coords[[id_col[1]]]
    } else {
      # Fall back to rownames of the Seurat object (order matches GetTissueCoordinates)
      obj_barcodes
    }
  }
  
  cell_ids <- .get_ids(cell_coords, colnames(srt_cell))
  bin_ids  <- .get_ids(bin_coords,  colnames(srt_bin))
  
  if (verbose) {
    message(sprintf("  Cells: %d | Bins: %d", length(cell_ids), length(bin_ids)))
  }
  
  # ── 3. Resolve annotation columns ──────────────────────────────────────────
  # Seurat internal columns to always exclude
  seurat_internals <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA",
    "nCount_Spatial", "nFeature_Spatial",
    "nCount_VisiumHD", "nFeature_VisiumHD"
  )
  
  if (is.null(annotation_cols)) {
    annotation_cols <- setdiff(names(srt_bin@meta.data), seurat_internals)
    if (verbose) {
      message(sprintf(
        "No annotation_cols specified. Transferring: %s",
        paste(annotation_cols, collapse = ", ")
      ))
    }
  }
  
  missing_cols <- setdiff(annotation_cols, names(srt_bin@meta.data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in srt_bin@meta.data: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # ── 4. Nearest-neighbor search (k-d tree) ──────────────────────────────────
  if (verbose) message("Running nearest-neighbor search (FNN k-d tree)...")
  
  bin_mat  <- as.matrix(bin_coords[,  c("x", "y")])
  cell_mat <- as.matrix(cell_coords[, c("x", "y")])
  
  nn <- FNN::get.knnx(
    data  = bin_mat,   # reference: bin centroids
    query = cell_mat,  # query: cell centroids
    k     = 1
  )
  
  nn_bin_row <- nn$nn.index[, 1]   # row index into bin_coords
  nn_dist    <- nn$nn.dist[, 1]    # Euclidean distance
  
  # # ── 5. Apply max_dist filter ────────────────────────────────────────────────
  # too_far <- nn_dist > max_dist
  # if (verbose && any(too_far)) {
  #   message(sprintf(
  #     "  %d / %d cells (%.1f%%) exceeded max_dist=%.1f and will receive NA.",
  #     sum(too_far), length(too_far),
  #     100 * mean(too_far), max_dist
  #   ))
  # }
  
  # ── 6. Build transfer data frame ────────────────────────────────────────────
  matched_bin_ids        <- bin_ids[nn_bin_row]
  # matched_bin_ids[too_far] <- NA
  
  transfer_df <- data.frame(
    nn_bin_id = matched_bin_ids,
    nn_dist   = round(nn_dist, 4),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )
  
  # Pull annotation values from bin metadata
  bin_meta <- srt_bin@meta.data[, annotation_cols, drop = FALSE]
  
  for (col in annotation_cols) {
    vals <- bin_meta[[col]]
    
    # Map bin row index → annotation value
    transferred <- vals[nn_bin_row]
    
    # Preserve factor levels if applicable
    if (is.factor(vals)) {
      transferred <- factor(transferred, levels = levels(vals))
    }
    
    # # Apply max_dist mask
    # if (is.factor(transferred)) {
    #   transferred[too_far] <- NA
    # } else {
    #   transferred[too_far] <- NA
    # }
    
    transfer_df[[col]] <- transferred
  }
  
  # ── 7. Add to srt_cell metadata ─────────────────────────────────────────
  if (verbose) message("Adding transferred annotations to srt_cell@meta.data...")
  
  # Align to Seurat cell order (cell_ids may not match colnames order)
  transfer_df_aligned <- transfer_df[colnames(srt_cell), , drop = FALSE]
  
  srt_cell <- Seurat::AddMetaData(srt_cell, metadata = transfer_df_aligned)
  
  # ── 8. Summary ───────────────────────────────────────────────────────────────
  if (verbose) {
    n_matched <- sum(!is.na(transfer_df_aligned$nn_bin_id))
    message(sprintf(
      "\nTransfer complete.\n  Matched: %d / %d cells (%.1f%%)\n  Median NN distance: %.2f µm\n  Columns added: %s",
      n_matched, ncol(srt_cell),
      100 * n_matched / ncol(srt_cell),
      median(nn_dist),
      paste(c("nn_bin_id", "nn_dist", annotation_cols), collapse = ", ")
    ))
  }
  
  srt_cell
}




pathwayenrich_plot <- function(top_n, gsea_result, pvalue_show= F, plot = c("enrich", "network")[1], save.path = NULL){
  require(patchwork)
  require(ggplot2)
  require(enrichplot)
  top5_pathways <- gsea_result@result$ID[1:min(top_n, nrow(gsea_result@result))]
  if ("enrich" %in% plot){
  
  
  p <- gseaplot2(
    gsea_result,
    geneSetID = top5_pathways,
    title     = paste("Top", top_n, " Enriched Pathways"),
    base_size = 12,
    pvalue_table = pvalue_show, 
    color = as.character(pals::polychrome(top_n))
  )
  
  p[[1]] <- p[[1]] + theme(
    legend.position  = "bottom",
    legend.direction = "horizontal"
  )+ guides(color = guide_legend(ncol =1))
  
  # Reassemble with patchwork
  
  p <- patchwork::wrap_plots(p, ncol = 1)
  return(p)
  }
  
  if("network" %in% plot){
    p <- cnetplot(
      gsea_result,
      foldChange   = geneList,
      showCategory = top_n,
      node_label   = "all",# ensure all labels are shown
      size_item = 0.7,
      color_category     = "red"
    ) +
      scale_color_gradient2(
        name  = "Fold Change",
        low   = "#4575b4", mid = "white", high = "#d73027"
      ) +
      ggtitle("GSEA Concept Network") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
  }
  .return_genesetname <- function(id){
    test <- unique(sapply(strsplit(id, split = "_"), '[', 1))
    if (test[1] %in% unique(sapply(strsplit(Hall$gs_name, split = "_"), '[', 1))) return("Hall")
    if (test[1] %in% unique(sapply(strsplit(C6$gs_name, split = "_"), '[', 1))) return("C6")
    if (test[1] %in% unique(sapply(strsplit(C5$gs_name, split = "_"), '[', 1))) return("C5")    
  }
  if(!is.null(save.path)){
    name = paste0(save.path, .return_genesetname(top5_pathways), ".pdf")
    ggsave(name, width = 8, height = 8, plot = p)
  }
}

binarise_expression <- function(expr,
                                ref_cells  = NULL,
                                z_thresh   = 2,
                                min_expr   = 0,
                                plot_out   = NULL,
                                verbose    = TRUE) {
  require(mclust)
  stopifnot(is.numeric(expr), length(expr) > 0)
  
  gene_name <- if (!is.null(names(expr))) "gene" else "gene"
  
  # ── 1. log1p-transform for modelling ─────────────────────
  log_expr <- log1p(expr)
  
  # non-zero values only (zero-inflation is handled separately)
  nz_idx   <- expr > min_expr
  nz_log   <- log_expr[nz_idx]
  
  if (sum(nz_idx) < 10) {
    warning("Fewer than 10 non-zero cells – returning all zeros.")
    return(setNames(integer(length(expr)), names(expr)))
  }
  
  # ── 2. Characterise the background (low) component ───────
  
  if (!is.null(ref_cells)) {
    
    # --- Reference cluster supplied ---
    if (is.null(names(expr))) {
      stop("`expr` must be a *named* vector when `ref_cells` is provided.")
    }
    ref_in_data <- ref_cells[ref_cells %in% names(expr)]
    if (length(ref_in_data) == 0) {
      stop("None of the supplied `ref_cells` match names in `expr`.")
    }
    
    ref_log  <- log_expr[ref_in_data]
    ref_log  <- ref_log[ref_log > log1p(min_expr)]   # non-zero only
    
    if (length(ref_log) < 5) {
      warning("Fewer than 5 non-zero reference cells – falling back to mixture model.")
      ref_cells <- NULL   # trigger mixture model below
    } else {
      mu_low    <- mean(ref_log)
      sd_low    <- sd(ref_log)
      method    <- "reference_cluster"
      mc        <- NULL
    }
  }
  
  if (is.null(ref_cells)) {
    
    # --- Gaussian mixture model on non-zero log-values ---
    # mclust "V" = unequal variances; try 1 or 2 components
    mc <- tryCatch(
      mclust::Mclust(nz_log, G = 1:2, model = "V", verbose = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(mc) || mc$G == 1) {
      # Only one component found → use mean/sd of all non-zero values
      # and shift threshold based on skewness
      mu_low  <- mean(nz_log)
      sd_low  <- sd(nz_log)
      method  <- "single_component"
      message("Note: Only one mixture component detected. ",
              "Threshold set at mean + ", z_thresh, " * SD of non-zero values.")
    } else {
      # Take the lower-mean component as "background"
      low_comp <- which.min(mc$parameters$mean)
      mu_low   <- mc$parameters$mean[low_comp]
      sd_low   <- sqrt(mc$parameters$variance$sigmasq[
        min(low_comp, length(mc$parameters$variance$sigmasq))
      ])
      method   <- "mixture_model"
    }
  }
  
  # ── 3. Compute threshold (on log scale, back-transform) ──
  threshold_log <- mu_low + z_thresh * sd_low
  threshold_raw <- expm1(threshold_log)   # inverse of log1p
  
  if (verbose) {
    cat("── Binarisation summary ─────────────────────────────\n")
    cat("  Method          :", method, "\n")
    cat("  Background mean (log1p) :", round(mu_low, 3), "\n")
    cat("  Background SD   (log1p) :", round(sd_low, 3), "\n")
    cat("  z threshold     :", z_thresh, "\n")
    cat("  Threshold (log1p):", round(threshold_log, 3), "\n")
    cat("  Threshold (raw) :", round(threshold_raw, 3), "\n")
    n_pos <- sum(expr > threshold_raw)
    cat("  Positive cells  :", n_pos, "/", length(expr),
        sprintf("(%.1f%%)\n", 100 * n_pos / length(expr)))
    cat("─────────────────────────────────────────────────────\n")
  }
  
  # ── 4. Assign binary labels ───────────────────────────────
  binary <- as.integer(expr > threshold_raw)
  names(binary) <- names(expr)
  
  # ── 5. Diagnostic plot ────────────────────────────────────
  if (!is.null(plot_out)) {
    .make_diagnostic_plot(
      expr          = expr,
      log_expr      = log_expr,
      nz_log        = nz_log,
      binary        = binary,
      threshold_log = threshold_log,
      threshold_raw = threshold_raw,
      mu_low        = mu_low,
      sd_low        = sd_low,
      mc            = mc,
      ref_cells     = ref_cells,
      method        = method,
      z_thresh      = z_thresh,
      plot_out      = plot_out
    )
    if (verbose) cat("  Diagnostic plot :", plot_out, "\n")
  }
  
  return(binary)
}


# ============================================================
#  INTERNAL – diagnostic plot
# ============================================================

.make_diagnostic_plot <- function(expr, log_expr, nz_log, binary,
                                  threshold_log, threshold_raw,
                                  mu_low, sd_low, mc, ref_cells,
                                  method, z_thresh, plot_out) {
  
  df_all <- data.frame(
    log_expr = log_expr,
    binary   = factor(binary, levels = c(0, 1),
                      labels = c("Negative", "Positive"))
  )
  
  # Panel A – raw expression histogram
  p1 <- ggplot(data.frame(expr = expr), aes(x = expr)) +
    geom_histogram(bins = 80, fill = "#4a90d9", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = threshold_raw, colour = "#e74c3c",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = threshold_raw, y = Inf,
             label = paste0(" threshold\n ", round(threshold_raw, 2)),
             vjust = 1.4, hjust = -0.05, colour = "#e74c3c", size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Raw expression (all cells)",
         x = "Expression", y = "Count") +
    theme_bw(base_size = 11)
  
  # Panel B – log1p non-zero distribution with model overlay
  x_seq <- seq(min(nz_log) - 0.5, max(nz_log) + 0.5, length.out = 400)
  n_nz  <- length(nz_log)
  bw    <- diff(range(nz_log)) / 40          # approx histogram bin width
  
  dens_df <- data.frame(x = x_seq)
  
  if (!is.null(mc) && mc$G == 2) {
    prop  <- mc$parameters$pro
    means <- mc$parameters$mean
    vars  <- mc$parameters$variance$sigmasq
    if (length(vars) == 1) vars <- rep(vars, 2)
    dens_df$low  <- prop[which.min(means)] *
      dnorm(x_seq, means[which.min(means)], sqrt(vars[which.min(means)])) *
      n_nz * bw
    dens_df$high <- prop[which.max(means)] *
      dnorm(x_seq, means[which.max(means)], sqrt(vars[which.max(means)])) *
      n_nz * bw
    dens_df$total <- dens_df$low + dens_df$high
    overlay <- TRUE
  } else {
    overlay <- FALSE
  }
  
  p2 <- ggplot(data.frame(nz_log = nz_log), aes(x = nz_log)) +
    geom_histogram(bins = 60, fill = "#7fb3d3", colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = threshold_log, colour = "#e74c3c",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = threshold_log, y = Inf,
             label = paste0(" threshold\n log1p=", round(threshold_log, 2)),
             vjust = 1.4, hjust = -0.05, colour = "#e74c3c", size = 3)
  
  if (overlay) {
    p2 <- p2 +
      geom_line(data = dens_df, aes(x = x, y = low),
                colour = "#2ecc71", linewidth = 0.9, linetype = "solid") +
      geom_line(data = dens_df, aes(x = x, y = high),
                colour = "#e67e22", linewidth = 0.9, linetype = "solid") +
      geom_line(data = dens_df, aes(x = x, y = total),
                colour = "black", linewidth = 0.7, linetype = "dotted")
  }
  
  p2 <- p2 +
    labs(title = paste0("log1p expression – non-zero cells (method: ", method, ")"),
         subtitle = paste0("Background: μ=", round(mu_low, 2),
                           "  σ=", round(sd_low, 2),
                           "  threshold = μ + ", z_thresh, "σ"),
         x = "log1p(Expression)", y = "Count") +
    theme_bw(base_size = 11)
  
  # Panel C – binarised result
  tbl <- table(binary)
  df_bar <- data.frame(
    label = c("Negative (0)", "Positive (1)"),
    count = as.integer(tbl[c("0", "1")])
  )
  df_bar$pct <- 100 * df_bar$count / sum(df_bar$count)
  
  p3 <- ggplot(df_bar, aes(x = label, y = count, fill = label)) +
    geom_col(width = 0.5, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
              vjust = 1, size = 3.5) +
    scale_fill_manual(values = c("Negative (0)" = "#95a5a6",
                                 "Positive (1)" = "#e74c3c")) +
    labs(title = "Binarised cell counts", x = NULL, y = "Count") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 10))
  
  combined <- (p1 | p2) / p3 +
    plot_annotation(
      title    = "Gene expression binarisation diagnostics",
      subtitle = paste0("Total cells: ", length(expr),
                        "  |  Threshold (raw): ", round(threshold_raw, 3)),
      theme    = theme(plot.title = element_text(face = "bold", size = 13))
    )
  print(combined)
  
  ggsave(plot_out, combined, width = 12, height = 8, dpi = 150)
}

plot_cnv_heatmap <- function(mat, labels, title = "CNV Subclusters",
                             out_file = NULL, chr_annot = NULL) {
  
  # Order cells by cluster
  cell_order <- order(labels)
  mat_ord    <- mat[cell_order, ]
  labs_ord   <- labels[cell_order]
  
  # Colour scale: diverging around 0 (log2-ratio) or around 1 (ratio)
  centre_val <- median(mat)
  col_fun <- circlize::colorRamp2(
    c(centre_val - 1, centre_val, centre_val + 1),
    c("#2166ac", "white", "#b2182b")
  )
  
  # Row annotation – cluster identity
  n_cl   <- length(unique(labs_ord))
  cl_col <- setNames(
    colorRampPalette(brewer.pal(min(n_cl, 9), "Set1"))(n_cl),
    sort(unique(labs_ord))
  )
  row_ann <- ComplexHeatmap::rowAnnotation(
    Cluster = factor(labs_ord),
    col     = list(Cluster = cl_col),
    show_annotation_name = FALSE
  )
  
  # Optional column annotation – chromosome bands
  col_ann <- NULL
  if (!is.null(chr_annot)) {
    genes_in_mat <- colnames(mat_ord)
    chr_vec      <- chr_annot[genes_in_mat]
    chr_vec[is.na(chr_vec)] <- "unknown"
    chr_col <- setNames(
      colorRampPalette(brewer.pal(8, "Pastel2"))(length(unique(chr_vec))),
      unique(chr_vec)
    )
    col_ann <- ComplexHeatmap::HeatmapAnnotation(
      Chromosome = chr_vec,
      col        = list(Chromosome = chr_col),
      show_annotation_name = FALSE
    )
  }
  
  ht <- ComplexHeatmap::Heatmap(
    mat_ord,
    name              = "Rel. CN",
    col               = col_fun,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    left_annotation   = row_ann,
    top_annotation    = col_ann,
    row_split         = factor(labs_ord),
    column_title      = title,
    use_raster        = TRUE,
    raster_quality    = 2
  )
  
  if (!is.null(out_file)) {
    pdf(out_file, width = 14, height = 8)
    ComplexHeatmap::draw(ht)
    dev.off()
    message("Heatmap saved to: ", out_file)
  } else {
    ComplexHeatmap::draw(ht)
  }
  
  invisible(ht)
}


cellhighlight_imagedim <- function(srt, cond){
  srt$highlight <- ifelse(cond, 1, 0)
  p <- ImageDimPlot(srt, group.by = "highlight")
  return(p)
}


# Spatially map single/multi-gene expression on a dark background: point colour
# AND point size both scale continuously with expression (white = low, firebrick
# = high). Mimics ImageFeaturePlot()'s spatial-expression role via a custom
# geom_point renderer so point size can be driven by expression too.
dark_feature_plot <- function(srt, features, assay = "SpaNorm", layer = "data",
                               image = NULL, order_high_on_top = TRUE,
                               pt_size_range = c(0.05, 1.2),
                               low_color = "#080144", high_color = "white",
                               bg_color = "black", log1p_transform = FALSE,
                               ncol = NULL, size_legend = FALSE,
                               title_wrap_width = 10, base_title_size = 12,
                               wrapped_title_size = 8) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)

  # ── 1. Tissue coordinates, aligned to srt cell order by barcode ──────────
  coords <- tryCatch(
    Seurat::GetTissueCoordinates(srt, image = image, which = "centroids"),
    error = function(e) {
      if (all(c("x_centroid", "y_centroid") %in% colnames(srt@meta.data))) {
        data.frame(cell = colnames(srt), x = srt$x_centroid, y = srt$y_centroid,
                   stringsAsFactors = FALSE)
      } else {
        stop("dark_feature_plot: no FOV/image on `srt` and no x_centroid/",
             "y_centroid meta.data columns to fall back on.")
      }
    }
  )
  id_col <- intersect(c("cell", "barcode", "id"), names(coords))
  if (length(id_col) == 0)
    stop("dark_feature_plot: coordinates table has no cell/barcode column.")
  idx <- match(colnames(srt), coords[[id_col[1]]])
  if (anyNA(idx))
    stop(sprintf("dark_feature_plot: %d cells have no matching coordinate row.",
                 sum(is.na(idx))))
  coords <- coords[idx, ]

  # ── 2. Expression ─────────────────────────────────────────────────────────
  expr_mat <- GetAssayData(srt, assay = assay, layer = layer)
  features <- features[features %in% rownames(expr_mat)]
  if (length(features) == 0)
    stop("dark_feature_plot: none of `features` found in assay '", assay, "'.")

  dark_theme <- theme_void() +
    theme(
      plot.background   = element_rect(fill = bg_color, colour = bg_color),
      panel.background  = element_rect(fill = bg_color, colour = bg_color),
      legend.background = element_rect(fill = bg_color, colour = NA),
      legend.key        = element_rect(fill = bg_color, colour = NA),
      legend.text       = element_text(colour = "white"),
      legend.title      = element_text(colour = "white"),
      plot.title        = element_text(colour = "white", hjust = 0.5, face = "bold")
    )

  plots <- lapply(features, function(feat) {
    v <- expr_mat[feat, ]
    if (log1p_transform) v <- log1p(v)
    df <- data.frame(x = coords$x, y = coords$y, expr = v)
    if (order_high_on_top) df <- df[base::order(df$expr), ]  # high expr drawn last (on top)

    wrapped   <- nchar(feat) > title_wrap_width
    title_txt <- if (wrapped) paste(strwrap(feat, width = title_wrap_width), collapse = "\n") else feat
    title_sz  <- if (wrapped) wrapped_title_size else base_title_size

    ggplot(df, aes(x, y, colour = expr, size = expr)) +
      geom_point(shape = 16) +
      scale_colour_gradient(low = low_color, high = high_color, name = "Expression") +
      scale_size_continuous(range = pt_size_range,
                            guide = if (size_legend) "legend" else "none") +
      scale_y_reverse() +
      coord_fixed() +
      dark_theme +
      theme(plot.title = element_text(size = title_sz)) +
      ggtitle(title_txt)
  })
  names(plots) <- features

  patchwork::wrap_plots(plots, ncol = ncol) &
    theme(plot.background = element_rect(fill = bg_color, colour = bg_color))
}



srt2anndata <- function(srt,
                        count_assay = "Spatial",
                        data_assay  = "SpaNorm",
                        save_name,
                        svg_path    = "SVGs.Rds") {
  require(Seurat)
  require(dplyr)
  require(anndataR, lib.loc = "~/R_Library/4.5")
  require(Matrix)

  # 1. Matrices: align to the SpaNorm gene set, transpose to (cells × genes),
  #    and coerce to CsparseMatrix (anndataR's HDF5 writer rejects dgRMatrix
  #    that Matrix::t() can produce, tripping .validate_aligned_array()).
  counts_full  <- GetAssayData(srt, assay = count_assay, layer = "counts")
  data_full    <- GetAssayData(srt, assay = data_assay,  layer = "data")
  common_genes <- intersect(rownames(data_full), rownames(counts_full))
  counts_mat   <- as(Matrix::t(counts_full[common_genes, ]), "CsparseMatrix")
  data_mat     <- as(Matrix::t(data_full[common_genes, ]),   "CsparseMatrix")

  # 2. obs: sanitise column names — HDF5 dataset names cannot contain '/'.
  metadata <- srt[[]]
  colnames(metadata) <- gsub("/", "_", colnames(metadata), fixed = TRUE)
  colnames(metadata) <- gsub("[|[:space:]]+", "_", colnames(metadata))
  colnames(metadata) <- gsub("_+", "_", colnames(metadata))

  # 3. obsm: export every reduction as X_<sanitised-name> (scanpy convention).
  red_names <- Reductions(srt)
  obsm_list <- setNames(
    lapply(red_names, function(r) Embeddings(srt, r)),
    paste0("X_", gsub("[.-]", "_", tolower(red_names)))
  )

  # 4. var: flag SVGs (fdr < 0.05) as highly_variable.
  gene_names <- colnames(data_mat)
  var_df <- data.frame(
    highly_variable = rep(FALSE, length(gene_names)),
    row.names       = gene_names
  )
  if (!is.null(svg_path) && file.exists(svg_path)) {
    svgs_df <- as.data.frame(readRDS(svg_path))
    hvg <- svgs_df %>% filter(!is.na(svg.fdr), svg.fdr < 0.05) %>% pull(symbol)
    var_df$highly_variable <- gene_names %in% hvg
    cat(sum(var_df$highly_variable), "SVGs (fdr<0.05) marked as highly_variable\n")
    write.csv(svgs_df, sub("\\.Rds$", ".csv", svg_path), row.names = TRUE)
  } else {
    # No SVGs supplied: fall back to the variable features of the data assay.
    hvg <- tryCatch(VariableFeatures(srt, assay = data_assay),
                    error = function(e) character(0))
    var_df$highly_variable <- gene_names %in% hvg
    cat("No SVGs supplied —", sum(var_df$highly_variable),
        "variable features from assay '", data_assay,
        "' marked as highly_variable\n", sep = "")
  }

  # 5. Spatial coords → obsm[['spatial']] (Squidpy/Scanpy convention).
  # Prefer FOV centroids; after cross-sample merge those are gone, so fall
  # back to x_centroid / y_centroid stashed in @meta.data.
  spatial_coords <- tryCatch({
    coords_df <- GetTissueCoordinates(srt)
    rownames(coords_df) <- coords_df$cell
    as.matrix(coords_df[rownames(metadata), c("x", "y")])
  }, error = function(e) NULL)
  if (is.null(spatial_coords) &&
      all(c("x_centroid", "y_centroid") %in% colnames(metadata))) {
    spatial_coords <- as.matrix(metadata[, c("x_centroid", "y_centroid")])
  }
  if (!is.null(spatial_coords)) obsm_list[["spatial"]] <- spatial_coords

  # 6. Build AnnData.
  adata <- AnnData(
    X      = data_mat,
    layers = list(counts = counts_mat),
    obs    = metadata,
    var    = var_df,
    obsm   = obsm_list
  )

  # 7. Write h5ad — disable HDF5 locking (Lustre) and remove any stale file.
  Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
  out_path <- normalizePath(paste0(save_name, ".h5ad"), mustWork = FALSE)
  if (file.exists(out_path)) file.remove(out_path)
  adata$write_h5ad(out_path, mode = "w")
  cat("All Done ->", out_path, "\n")
}


filter_artefacts_knn <- function(srt,
                                 k              = 10,      # number of neighbours
                                 min_neighbours  = 3,      # min neighbours to be retained
                                 image_name     = NULL) {
  
  img    <- image_name %||% Images(srt)[1]
  coords <- GetTissueCoordinates(srt, image = img)
  xy     <- as.matrix(coords[, c("x", "y")])
  
  # ── 1. Build k-NN graph ───────────────────────────────────────────────────
  knn_result  <- dbscan::kNN(xy, k = k)
  nn_dists    <- knn_result$dist   # matrix: cells × k
  
  # ── 2. Count how many neighbours are within a data-driven distance cutoff ─
  #   Use the median of all k-NN distances × 3 as a loose "in-tissue" threshold
  dist_cutoff <- median(nn_dists, na.rm = TRUE) * 3
  neighbour_counts <- rowSums(nn_dists <= dist_cutoff)
  
  # ── 3. Add metadata and flag ──────────────────────────────────────────────
  srt$n_spatial_neighbours <- neighbour_counts
  srt$is_artefact          <- neighbour_counts < min_neighbours
  
  message(sprintf(
    "kNN filter: %d cells flagged as artefacts (neighbour count < %d, cutoff = %.1f px)",
    sum(srt$is_artefact), min_neighbours, dist_cutoff
  ))
  
  # ── 4. Return clean object ────────────────────────────────────────────────
  return(subset(srt, cells = colnames(srt)[!srt$is_artefact]))
}


# Plot pre-computed UCell meta-program scores (one PNG per sheet per reduction).
# Scores must already be in obj@meta.data, computed via:
#   res <- ucell_score(obj, meta_programs_unlist, "_meta"); meta_cols <- res$cols
# meta_programs is the original nested list (sheet → program → genes) used only
# to group columns by sheet; meta_cols is the vector of column names to plot.
score_meta_programs <- function(obj, meta_programs, meta_cols, out_dir,
                                reductions = "umap",
                                tag = "") {
  require(Seurat)
  require(ggplot2)

  for (sheet in names(meta_programs)) {
    # Columns for this sheet: ucell_score named them _meta<SheetName>.<ProgramName>_UCell
    sheet_pat  <- paste0("^", make.names(sheet))
    sheet_cols <- grep(sheet_pat, meta_cols, value = TRUE)
    sheet_cols <- sheet_cols[sheet_cols %in% colnames(obj@meta.data)]
    if (length(sheet_cols) == 0) next

    ncol_g <- min(3L, length(sheet_cols))
    nrow_g <- ceiling(length(sheet_cols) / ncol_g)
    base   <- make.names(sheet)
    plots_s <- lapply(sheet_cols, function(f) {
      v   <- obj@meta.data[[f]]
      mid <- (max(v, na.rm = TRUE) + min(v, na.rm = TRUE)) * 0.3
      tt <- gsub(paste0(sheet, "."), fix = T, "", f)
      tt <- gsub(".", " ",fix = T, tt)
      tt <- gsub("_meta_UCell", "", tt)
      FeaturePlot(obj, features = f, reduction = reductions) +
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                              midpoint = mid) + ggtitle(tt)
    })
    g <- patchwork::wrap_plots(plots_s, ncol = ncol_g)
    ggsave(file.path(out_dir, paste0("meta_", base, tag, reductions, ".png")),
           plot = g, width = ncol_g * 4, height = nrow_g * 3.5,
           dpi = 300, limitsize = FALSE)
  }
}


# Run GSEA for every cluster against each gene-set collection and save the result.
# `gene_sets` is a named list (e.g. list(Hallmark=Hall, C6=C6, C5=C5)).
run_gsea_panel <- function(deg, gene_sets, out_path) {
  if (is.null(deg) || nrow(deg) == 0) {
    message("run_gsea_panel: empty DEG, skipping")
    return(invisible(NULL))
  }
  deg <- deg %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC))
  if (nrow(deg) == 0) {
    message("run_gsea_panel: no DEG passing p_val_adj < 0.05, skipping")
    return(invisible(NULL))
  }
  gene_list <- lapply(split(deg, deg$cluster),
                      function(x) setNames(x$avg_log2FC, x$gene))

  enrich_all <- lapply(gene_sets, function(geneset) {
    out <- lapply(names(gene_list), function(cn) {
      tryCatch(clusterProfiler::GSEA(gene_list[[cn]], TERM2GENE = geneset),
               error = function(e) {
                 message(sprintf("GSEA cluster '%s' failed: %s",
                                 cn, conditionMessage(e)))
                 NULL
               })
    })
    names(out) <- names(gene_list)
    out
  })
  saveRDS(enrich_all, out_path)
  invisible(enrich_all)
}


# Save Image + UMAP FeaturePlots of the top-N DEG per cluster.
plot_top_deg <- function(srt, deg, out_dir, prefix = "spanorm_DEG",
                         n = 5, pct_diff_min = 0.2) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  if (is.null(deg) || nrow(deg) == 0) return(invisible(NULL))
  top_genes <- deg %>%
    dplyr::filter(abs(pct.1 - pct.2) > pct_diff_min) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = n) %>%
    dplyr::pull(gene)
  if (length(top_genes) == 0) return(invisible(NULL))

  f_img <- ImageFeaturePlot(srt, top_genes, size = 0.3)
  ggsave(file.path(out_dir, paste0(prefix, "_ImageFeaturePlot.png")),
         plot = f_img + plot_layout(ncol = 5), limitsize = FALSE,
         width = 15, height = 50, dpi = 350)
  f_um <- FeaturePlot(srt, top_genes)
  ggsave(file.path(out_dir, paste0(prefix, "_FeaturePlot.png")),
         plot = f_um, limitsize = FALSE, width = 15, height = 30, dpi = 350)
  invisible(top_genes)
}
