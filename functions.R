spatial_plot <- function(srt, outdir, name) {
  require(Seurat)
  require(ggplot2)
  require(pals)
  require(patchwork)
  
  tryCatch(
    {
      f1 <- ImageDimPlot(srt, cols = "polychrome",border.color = "#00000000", size = 0.3)
      f2 <- ImageDimPlot(srt,  split.by = "seurat_clusters", cols = "polychrome",border.color = "#00000000", size = 0.3)
      ggsave(plot = f1 + f2 + plot_layout(ncol = 2), paste0(outdir, name, "_ImageDimPlot.png"), width = 10, height = 5, dpi = 350)
      
      f1 <- ImageDimPlot(srt, group.by = "SR_Cluster", cols = "polychrome",border.color = "#00000000", size = 0.3)
      f2 <- ImageDimPlot(srt, group.by = "SR_Cluster",  split.by = "SR_Cluster", cols = "polychrome",border.color = "#00000000", size = 0.3)
      ggsave(plot = f1 + f2 + plot_layout(ncol = 2), paste0(outdir, name, "SR_Cluster_ImageDimPlot.png"), width = 10, height = 5, dpi = 350)
      
      
      f3 <- ImageFeaturePlot(srt, c("AR", "FOLH1", "ASCL1", "COL1A1","cell_area", "nucleus_area","nFeature_Spatial", "nCount_Spatial"), size = 0.3)
      ggsave(plot = f3 +  plot_layout(ncol = 4), paste0(outdir, name, "_ImageFeaturePlot.png"), width = 13, height = 6, dpi = 350)
      
      g1 <- DimPlot(srt, cols = "polychrome", label = T)
      g2 <- DimPlot(srt,  group.by = "SR_Cluster", cols = "polychrome", label = T)
      ggsave(plot = g1 + g2 + plot_layout(ncol = 2), paste0(outdir, name, "_DimPlot.png"), width = 10, height = 5, dpi = 350)
      
      g3 <- FeaturePlot(srt, c("AR", "FOLH1", "ASCL1", "COL1A1","cell_area", "nucleus_area","nFeature_Spatial", "nCount_Spatial"))
      ggsave(plot = g3 + plot_layout(ncol = 4), paste0(outdir, name, "_FeaturePlot.png"), width = 12, height = 5, dpi = 350)
      
      p1 <- VlnPlot(srt, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, pt.size = 0, group.by = "seurat_clusters")
      p2 <- VlnPlot(srt, features = c("AR", "FOLH1", "ASCL1", "COL1A1"), ncol = 4, pt.size = 0, group.by = "seurat_clusters")        # the function you want to run first
      ggsave(plot = (p1|p2) + plot_layout(width = c(1, 2)), paste0(outdir, name, "_VlnPlot.png"), width = 15, height = 4, dpi = 350)
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

do.spanorm <- function(srt){
  require(SpaNorm)
  
  countMat <- GetAssayData(srt, assay = "Spatial", layer = "counts")
  
  spatial_coord <- GetTissueCoordinates(srt)[, c("x", "y")]
  rownames(spatial_coord) <- colnames(countMat)
  set.seed(1)
  spe <- createSPEObject(countMat, spatial_locs= spatial_coord, cell_metadata = spatial_coord, normalize = FALSE)
  spe <- SpaNorm(spe)
  spe <- SpaNormSVG(spe)
  saveRDS(rowData(spe), "SVGs.Rds")
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

return_genesetname <- function(id){
  test <- unique(sapply(strsplit(id, split = "_"), '[', 1))
  if (test %in% unique(sapply(strsplit(Hall$gs_name, split = "_"), '[', 1))) return("Hall")
  if (test %in% unique(sapply(strsplit(C6$gs_name, split = "_"), '[', 1))) return("C6")
  if (test %in% unique(sapply(strsplit(C5$gs_name, split = "_"), '[', 1))) return("C5")    
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
  )+ guides(color = guide_legend(ncol =3))
  
  # Reassemble with patchwork
  
  patchwork::wrap_plots(p, ncol = 1)}
  
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
    
    print(p)
  }
  
  if(!is.null(save.path)){
    name = paste0(save.path, return_genesetname(top5_pathways), ".pdf")
    ggsave(name, width = 8, height = 8, plot = p)
  }
}


