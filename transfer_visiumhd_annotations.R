# =============================================================================
# Transfer VisiumHD 16µm Bin Annotations to Cell-Segmented Cells
# Method: Point-in-Bin via barcode-decoded grid indices
#
# KEY FINDING FROM DIAGNOSTICS:
#   - Coordinates are in PIXELS, not µm
#   - Bin step size ≈ 58.33 pixels per 16µm bin  (= 16 / 0.2742 µm/px)
#   - X axis is FLIPPED (x decreases as grid col increases)
#   - half-bin tolerance = 8µm × (58.33px/16µm) ≈ 29.17 pixels
#
# Strategy: decode the grid row/col from each bin's barcode
#   (format: s_016um_RRRRR_CCCCC-1), snap each cell centroid to the
#   nearest grid index using the inferred pixel-per-bin step, then
#   join on (grid_row, grid_col). This is exact and O(n).
#
# Input:  Two Seurat objects — VisiumHD 16µm bins + cell segmentation
# Dependencies: Seurat, SeuratObject, data.table
# =============================================================================


#' Transfer VisiumHD 16µm bin annotations to segmented cells
#'
#' Decodes grid row/col from bin barcodes (format: s_016um_RRRRR_CCCCC-1),
#' snaps each cell centroid to the nearest bin grid position using the
#' inferred pixel step size, and joins on exact grid index. This correctly
#' handles non-uniform coordinate spacing and flipped axes.
#'
#' @param seurat_cells    Seurat object from cell segmentation. Centroids read
#'                        from @images[[fov]]@boundaries$centroids@coords
#' @param seurat_bins     Seurat object from VisiumHD 16µm binning.
#' @param annotation_cols Character vector of metadata columns in seurat_bins
#'                        to transfer. NULL = all non-internal columns.
#' @param bin_size_um     Bin side length in µm. Default: 16.
#' @param cells_fov       FOV name in seurat_cells. Default: first FOV.
#' @param bins_fov        FOV name in seurat_bins. Default: first FOV.
#' @param max_offset_frac Maximum allowed offset from bin center as a fraction
#'                        of the bin step size. Default: 0.5 (exactly one bin).
#'                        Increase slightly (e.g. 0.51) to catch floating-point
#'                        edge cases at bin boundaries.
#' @param verbose         Print progress and summary. Default: TRUE.
#'
#' @return seurat_cells with annotation columns added to meta.data, plus QC
#'         columns: matched_bin_id, bin_offset_x_px, bin_offset_y_px,
#'         matched_grid_row, matched_grid_col.
#'
#' @examples
#' srt_cell <- transfer_visiumhd_point_in_bin(
#'   seurat_cells    = srt_cell,
#'   seurat_bins     = srt_bin,
#'   annotation_cols = c("seurat_clusters", "BANKSY_0.2_snn_res.1"),
#' )
#' table(is.na(srt_cell$matched_bin_id))
#' ImageDimPlot(srt_cell, group.by = "seurat_clusters")

transfer_visiumhd_point_in_bin <- function(
    seurat_cells,
    seurat_bins,
    annotation_cols   = NULL,
    bin_size_um       = 16,
    cells_fov         = NULL,
    bins_fov          = NULL,
    max_offset_frac   = 0.5,
    verbose           = TRUE
) {
  
  # ── 0. Packages ─────────────────────────────────────────────────────────────
  for (pkg in c("Seurat", "SeuratObject", "data.table")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Install required package: install.packages('%s')", pkg))
  }
  
  # ── 1. Resolve FOV names ────────────────────────────────────────────────────
  .resolve_fov <- function(obj, fov_arg, label) {
    avail <- names(obj@images)
    if (!length(avail)) stop(sprintf("No FOV slots in %s.", label))
    if (is.null(fov_arg)) {
      fov_arg <- avail[1]
      if (verbose) message(sprintf("[%s] FOV: '%s'", label, fov_arg))
    } else if (!fov_arg %in% avail) {
      stop(sprintf("FOV '%s' not in %s. Available: %s",
                   fov_arg, label, paste(avail, collapse = ", ")))
    }
    fov_arg
  }
  
  cells_fov <- .resolve_fov(seurat_cells, cells_fov, "seurat_cells")
  bins_fov  <- .resolve_fov(seurat_bins,  bins_fov,  "seurat_bins")
  
  # ── 2. Extract cell centroids ───────────────────────────────────────────────
  if (verbose) message("Extracting cell centroids...")
  
  cell_boundaries <- seurat_cells@images[[cells_fov]]@boundaries
  if (!"centroids" %in% names(cell_boundaries))
    stop(sprintf("No 'centroids' in @boundaries. Found: %s",
                 paste(names(cell_boundaries), collapse = ", ")))
  
  cell_mat <- cell_boundaries$centroids@coords  # numeric matrix, rows = cells
  
  if (nrow(cell_mat) != ncol(seurat_cells))
    stop(sprintf("Centroid rows (%d) != cells in Seurat object (%d).",
                 nrow(cell_mat), ncol(seurat_cells)))
  
  cells_dt <- data.table::data.table(
    barcode = colnames(seurat_cells),
    cx      = cell_mat[, "x"],
    cy      = cell_mat[, "y"]
  )
  if (verbose) message(sprintf("  %d cells.", nrow(cells_dt)))
  
  # ── 3. Extract bin centroids + decode grid indices from barcodes ─────────────
  if (verbose) message("Extracting bin centroids and decoding grid indices...")
  
  bin_coords <- tryCatch(
    {
      bc <- Seurat::GetTissueCoordinates(seurat_bins, image = bins_fov, which = "centroids")
      if (!"x"    %in% names(bc)) names(bc)[1] <- "x"
      if (!"y"    %in% names(bc)) names(bc)[2] <- "y"
      if (!"cell" %in% names(bc)) bc$cell <- colnames(seurat_bins)
      bc
    },
    error = function(e) {
      if (verbose) message("  GetTissueCoordinates failed; reading @boundaries directly.")
      mat <- seurat_bins@images[[bins_fov]]@boundaries$centroids@coords
      data.frame(x = mat[,"x"], y = mat[,"y"], cell = colnames(seurat_bins),
                 stringsAsFactors = FALSE)
    }
  )
  
  # Decode grid row/col from barcode: s_016um_RRRRR_CCCCC-1
  bc_vec  <- bin_coords$cell
  matches <- regmatches(bc_vec, regexpr("(\\d+)_(\\d+)(?=-\\d)", bc_vec, perl = TRUE))
  
  if (length(matches) != nrow(bin_coords) || any(matches == "")) {
    stop(paste(
      "Could not parse grid row/col from bin barcodes.",
      "Expected format: s_016um_RRRRR_CCCCC-1.",
      sprintf("Example barcode: %s", bc_vec[1])
    ))
  }
  
  split_idx   <- strsplit(matches, "_")
  grid_rows   <- as.integer(sapply(split_idx, `[`, 1))
  grid_cols   <- as.integer(sapply(split_idx, `[`, 2))
  
  bins_dt <- data.table::data.table(
    bin_barcode = bc_vec,
    bx          = as.numeric(bin_coords$x),
    by          = as.numeric(bin_coords$y),
    grid_row    = grid_rows,
    grid_col    = grid_cols
  )
  
  if (verbose) message(sprintf(
    "  %d bins. Grid row: %d–%d, col: %d–%d.",
    nrow(bins_dt),
    min(grid_rows), max(grid_rows),
    min(grid_cols), max(grid_cols)
  ))
  
  # ── 4. Infer pixel step size per grid unit ───────────────────────────────────
  # Use a single representative row and col (the median) to get clean steps,
  # avoiding gaps in tissue coverage that cause large spurious diffs.
  if (verbose) message("Inferring pixel-per-bin step size...")
  
  .infer_step <- function(bins_dt, fix_var, sort_var, coord_var) {
    med_fix  <- as.integer(median(bins_dt[[fix_var]]))
    sub      <- bins_dt[bins_dt[[fix_var]] == med_fix, ]
    # Try median row/col; fall back to adjacent rows/cols if fewer than 5 bins
    for (attempt in 0:10) {
      sub <- bins_dt[bins_dt[[fix_var]] == (med_fix + attempt), ]
      if (nrow(sub) >= 5) break
      sub <- bins_dt[bins_dt[[fix_var]] == (med_fix - attempt), ]
      if (nrow(sub) >= 5) break
    }
    sub  <- sub[order(sub[[sort_var]]), ]
    diffs <- diff(sub[[coord_var]])
    # Filter out large gaps (missing tissue bins)
    diffs <- diffs[abs(diffs) < 3 * median(abs(diffs))]
    median(diffs)  # signed: preserves axis direction
  }
  
  step_x <- .infer_step(bins_dt, "grid_row", "grid_col", "bx")  # px per col step
  step_y <- .infer_step(bins_dt, "grid_col", "grid_row", "by")  # px per row step
  
  if (verbose) {
    message(sprintf("  Pixel step: x = %.4f px/col, y = %.4f px/row", step_x, step_y))
    message(sprintf("  Scale:      %.4f µm/px (x), %.4f µm/px (y)",
                    bin_size_um / abs(step_x), bin_size_um / abs(step_y)))
    if (step_x < 0) message("  Note: X axis is flipped (x decreases as col increases).")
    if (step_y < 0) message("  Note: Y axis is flipped (y decreases as row increases).")
  }
  
  half_step_x <- abs(step_x) * max_offset_frac
  half_step_y <- abs(step_y) * max_offset_frac
  
  # ── 5. Snap each cell centroid to the nearest grid (row, col) ───────────────
  # Use the bin grid origin (row=0, col=0) as anchor.
  # Find the bin at row=0 (or min row) to get the origin coordinates.
  origin_col <- bins_dt[grid_row == min(grid_row) & grid_col == min(grid_col)]
  if (nrow(origin_col) == 0) {
    # Fall back: use the bin with min grid_row and min grid_col separately
    origin_bx <- bins_dt$bx[which.min(bins_dt$grid_col)]
    origin_by <- bins_dt$by[which.min(bins_dt$grid_row)]
  } else {
    origin_bx <- origin_col$bx[1]
    origin_by <- origin_col$by[1]
  }
  
  # Snap: grid_col = round((cx - origin_bx) / step_x)
  #        grid_row = round((cy - origin_by) / step_y)
  cells_dt[, snapped_col := round((cx - origin_bx) / step_x)]
  cells_dt[, snapped_row := round((cy - origin_by) / step_y)]
  
  # Clip to valid grid range (cells outside tissue get clamped but will fail
  # the offset check below)
  cells_dt[, snapped_col := pmax(min(bins_dt$grid_col),
                                 pmin(max(bins_dt$grid_col), snapped_col))]
  cells_dt[, snapped_row := pmax(min(bins_dt$grid_row),
                                 pmin(max(bins_dt$grid_row), snapped_row))]
  
  # ── 6. Join cells to bins on (grid_row, grid_col) ───────────────────────────
  if (verbose) message("Joining cells to bins on grid indices...")
  
  data.table::setkey(bins_dt,  grid_row, grid_col)
  data.table::setkey(cells_dt, snapped_row, snapped_col)
  
  # Right join: all cells kept; unmatched get NA bin columns
  matched_dt <- bins_dt[
    cells_dt,
    on = c(grid_row = "snapped_row", grid_col = "snapped_col")
  ]
  
  # ── 7. Offset check — discard cells not truly inside the matched bin ─────────
  matched_dt[, bin_offset_x_px := cx - bx]
  matched_dt[, bin_offset_y_px := cy - by]
  
  matched_dt[, in_bin := !is.na(bx) &
               (abs(bin_offset_x_px) <= half_step_x) &
               (abs(bin_offset_y_px) <= half_step_y)]
  
  n_outside <- sum(!matched_dt$in_bin & !is.na(matched_dt$bin_barcode), na.rm = TRUE)
  if (n_outside > 0 && verbose)
    message(sprintf("  %d cells snapped to a bin but failed offset check → NA.", n_outside))
  
  matched_dt[in_bin == FALSE, bin_barcode := NA_character_]
  
  # ── 8. Resolve annotation columns ───────────────────────────────────────────
  seurat_internals <- c(
    "orig.ident",
    "nCount_RNA",      "nFeature_RNA",
    "nCount_Spatial",  "nFeature_Spatial",
    "nCount_VisiumHD", "nFeature_VisiumHD",
    "nCount_Xenium",   "nFeature_Xenium"
  )
  
  if (is.null(annotation_cols)) {
    annotation_cols <- setdiff(names(seurat_bins@meta.data), seurat_internals)
    if (verbose) message(sprintf(
      "Transferring: %s", paste(annotation_cols, collapse = ", ")
    ))
  }
  
  missing_cols <- setdiff(annotation_cols, names(seurat_bins@meta.data))
  if (length(missing_cols))
    stop("Not in seurat_bins@meta.data: ", paste(missing_cols, collapse = ", "))
  
  # ── 9. Merge bin metadata onto matched cells ─────────────────────────────────
  bin_meta    <- seurat_bins@meta.data[, annotation_cols, drop = FALSE]
  bin_meta$bin_barcode <- rownames(bin_meta)
  bin_meta_dt <- data.table::as.data.table(bin_meta)
  
  data.table::setkey(bin_meta_dt, bin_barcode)
  data.table::setkey(matched_dt,  bin_barcode)
  result_dt <- bin_meta_dt[matched_dt, on = "bin_barcode"]
  
  # ── 10. Build metadata aligned to Seurat cell order ─────────────────────────
  result_df             <- as.data.frame(result_dt)
  rownames(result_df)   <- result_df$barcode
  result_df             <- result_df[colnames(seurat_cells), , drop = FALSE]
  result_df$matched_bin_id    <- result_df$bin_barcode
  result_df$matched_grid_row  <- result_df$grid_row
  result_df$matched_grid_col  <- result_df$grid_col
  
  cols_to_add <- c("matched_bin_id", "matched_grid_row", "matched_grid_col",
                   "bin_offset_x_px", "bin_offset_y_px", annotation_cols)
  meta_to_add <- result_df[, cols_to_add, drop = FALSE]
  
  # ── 11. Add to Seurat and summarise ─────────────────────────────────────────
  seurat_cells <- Seurat::AddMetaData(seurat_cells, metadata = meta_to_add)
  
  if (verbose) {
    n_matched <- sum(!is.na(meta_to_add$matched_bin_id))
    n_total   <- ncol(seurat_cells)
    message(sprintf(
      "\n=== Transfer complete ===\n  Matched   : %d / %d cells (%.1f%%)\n  Unmatched : %d cells (NA)\n  Columns   : %s",
      n_matched, n_total, 100 * n_matched / n_total,
      n_total - n_matched,
      paste(cols_to_add, collapse = ", ")
    ))
  }
  
  seurat_cells
}


# =============================================================================
# QC Helper
# =============================================================================

#' QC plots after transfer_visiumhd_point_in_bin()
#'
#' @param seurat_cells   Seurat object after transfer
#' @param annotation_col Optional: one transferred column for spatial plot
#' @param bin_step_px    Pixel step size per bin (default 58.33). Used to draw
#'                       the reference box on the offset scatter.
#'
#' @return Named list: $offset_scatter, $match_rate, $spatial (optional)

plot_point_in_bin_qc <- function(seurat_cells,
                                 annotation_col = NULL,
                                 bin_step_px    = 58.33) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("install.packages('ggplot2')")
  
  meta     <- seurat_cells@meta.data
  half_bin <- bin_step_px / 2
  
  if (!"bin_offset_x_px" %in% names(meta))
    stop("QC columns missing. Run transfer_visiumhd_point_in_bin() first.")
  
  plots <- list()
  
  # 1. Offset scatter — tight square cloud within ±half_bin box expected
  matched <- meta[!is.na(meta$bin_offset_x_px), ]
  plots$offset_scatter <- ggplot2::ggplot(
    matched, ggplot2::aes(x = bin_offset_x_px, y = bin_offset_y_px)
  ) +
    ggplot2::geom_hex(bins = 60) +
    ggplot2::annotate("rect",
                      xmin = -half_bin, xmax = half_bin,
                      ymin = -half_bin, ymax = half_bin,
                      fill = NA, color = "red", linewidth = 0.8, linetype = "dashed"
    ) +
    ggplot2::scale_fill_viridis_c(option = "magma") +
    ggplot2::labs(
      title    = "Cell centroid offset from matched bin center (pixels)",
      subtitle = sprintf("Red box = ±%.1f px (one bin half-width). Points should cluster inside.", half_bin),
      x = "Offset X (px)", y = "Offset Y (px)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12)
  
  # 2. Match rate
  match_df <- data.frame(
    status = c("Matched", "Unmatched"),
    count  = c(sum(!is.na(meta$matched_bin_id)), sum(is.na(meta$matched_bin_id)))
  )
  plots$match_rate <- ggplot2::ggplot(
    match_df, ggplot2::aes(x = status, y = count, fill = status)
  ) +
    ggplot2::geom_col(width = 0.5, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.4, size = 4) +
    ggplot2::scale_fill_manual(
      values = c("Matched" = "#4CAF50", "Unmatched" = "#E57373")
    ) +
    ggplot2::labs(title = "Point-in-bin match rate", x = NULL, y = "Cell count") +
    ggplot2::theme_minimal(base_size = 12)
  
  # 3. Spatial plot
  if (!is.null(annotation_col)) {
    if (!annotation_col %in% names(meta)) {
      warning(sprintf("'%s' not in meta.data; skipping spatial plot.", annotation_col))
    } else {
      plots$spatial <- Seurat::ImageDimPlot(
        seurat_cells, group.by = annotation_col, size = 0.5
      ) + ggplot2::ggtitle(sprintf("Transferred: %s", annotation_col))
    }
  }
  
  plots
}
