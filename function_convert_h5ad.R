
srt2anndata <- function(srt,
                        count_assay = "Spatial",
                        data_assay = "SpaNorm",
                        save_name,
                        svg_path = "SVGs.Rds") {
  require(Seurat)
  require(dplyr)
  require(anndataR)
  require(qs2)
  require(Matrix)

  # 1. Extract matrices and transpose to (cells x features) for AnnData
  counts_mat <- t(srt[[count_assay]]$counts)
  data_mat   <- t(srt[[data_assay]]$data)

  # 2. Extract metadata
  metadata <- srt[[]]

  # 3. Build var with highly_variable from SVGs (fdr < 0.05)
  gene_names <- colnames(data_mat)
  var_df <- data.frame(row.names = gene_names, highly_variable = FALSE)
  if (!is.null(svg_path) && file.exists(svg_path)) {
    svgs_df <- as.data.frame(readRDS(svg_path))
    hvg <- svgs_df %>% filter(!is.na(svg.fdr), svg.fdr < 0.05) %>% pull(symbol)
    var_df$highly_variable <- gene_names %in% hvg
    cat(sum(var_df$highly_variable), "SVGs (fdr<0.05) marked as highly_variable\n")
    # Also write a plain CSV so Python can read it without R/pyreadr
    write.csv(svgs_df, sub("\\.Rds$", ".csv", svg_path), row.names = TRUE)
  } else {
    cat("SVGs.Rds not found at", svg_path, "— highly_variable set to FALSE\n")
  }

  # 4. Extract dimensionality reductions
  obsm_list <- list(
    X_pca          = Embeddings(srt, "pca"),
    X_umap         = Embeddings(srt, "umap"),
    X_banksy_pca   = Embeddings(srt, "banksy0.2.pca"),
    X_banksy_umap  = Embeddings(srt, "banksy0.2.umap")
  )

  # 5. Extract Spatial Coordinates
  coords_df <- GetTissueCoordinates(srt)
  rownames(coords_df) <- coords_df$cell
  spatial_coords <- as.matrix(coords_df[rownames(metadata), c("x", "y")])
  obsm_list[["spatial"]] <- spatial_coords

  # 6. Construct the AnnData object
  adata <- AnnData(
    X      = data_mat,
    layers = list(counts = counts_mat),
    obs    = metadata,
    var    = var_df,
    obsm   = obsm_list
  )

  # 7. Write to h5ad
  adata$write_h5ad(paste0(save_name, ".h5ad"), mode = "w")
  cat("All Done")
}
