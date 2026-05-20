# ============================================================
# subset_infercnv_for_plot()
# ------------------------------------------------------------
# Subset an inferCNV object for plotting, with the option to
# move some reference cells into the observation panel.
#
# Args:
#   infercnv_obj         : an inferCNV S4 object (post-run)
#   ref_keep_barcodes    : character vector of cell barcodes to
#                          keep as REFERENCES. If NULL, keep all
#                          original reference cells.
#   ref_to_move_barcodes : character vector of cell barcodes
#                          originally in references that should
#                          be MOVED into the observation panel.
#                          If NULL, no cells are moved.
#   obs_keep_barcodes    : character vector of observation cells
#                          to keep. If NULL, keep all original
#                          observation cells.
#   split_moved_by_group : logical. If TRUE, moved ref cells are
#                          split into separate obs groups named
#                          "moved_<original_ref_group>". If FALSE,
#                          they all go into one group called
#                          "moved_from_reference".
#   verbose              : logical. Print summary if TRUE.
#
# Returns:
#   A modified inferCNV object and  infercnv::plot_cnv().
# ============================================================

subset_infercnv_for_plot <- function(infercnv_obj,
                                     ref_keep_barcodes    = NULL,
                                     ref_to_move_barcodes = NULL,
                                     obs_keep_barcodes    = NULL,
                                     split_moved_by_group = FALSE,
                                     verbose              = TRUE) {
  
  # ---- 0. Basic input validation ------------------------------------
  stopifnot(methods::is(infercnv_obj, "infercnv"))
  
  all_cells <- colnames(infercnv_obj@expr.data)
  
  # Helper: get barcodes for a given group from an index list
  get_barcodes <- function(index_list, group_name, colnames_vec) {
    colnames_vec[index_list[[group_name]]]
  }
  
  # Helper: get all barcodes belonging to any group in an index list
  all_group_barcodes <- function(index_list, colnames_vec) {
    unlist(lapply(index_list, function(idx) colnames_vec[idx]),
           use.names = FALSE)
  }
  
  orig_ref_cells <- all_group_barcodes(
    infercnv_obj@reference_grouped_cell_indices, all_cells
  )
  orig_obs_cells <- all_group_barcodes(
    infercnv_obj@observation_grouped_cell_indices, all_cells
  )
  
  # ---- 1. Resolve defaults -----------------------------------------
  # If user didn't specify, keep everything in that category
  if (is.null(ref_keep_barcodes))    ref_keep_barcodes    <- orig_ref_cells
  if (is.null(ref_to_move_barcodes)) ref_to_move_barcodes <- character(0)
  if (is.null(obs_keep_barcodes))    obs_keep_barcodes    <- orig_obs_cells
  
  # ---- 2. Validate barcodes ----------------------------------------
  # Drop anything not present in the object
  ref_keep_barcodes    <- intersect(ref_keep_barcodes,    orig_ref_cells)
  ref_to_move_barcodes <- intersect(ref_to_move_barcodes, orig_ref_cells)
  obs_keep_barcodes    <- intersect(obs_keep_barcodes,    orig_obs_cells)
  
  # Cells cannot be both kept-as-ref AND moved-to-obs
  conflict <- intersect(ref_keep_barcodes, ref_to_move_barcodes)
  if (length(conflict) > 0) {
    stop("These barcodes are in both ref_keep_barcodes and ",
         "ref_to_move_barcodes: ",
         paste(utils::head(conflict, 5), collapse = ", "),
         if (length(conflict) > 5) " ..." else "")
  }
  
  if (length(ref_keep_barcodes) == 0 && length(ref_to_move_barcodes) == 0) {
    warning("No reference cells will remain after subsetting. ",
            "plot_cnv() may behave unexpectedly without references.")
  }
  
  # ---- 3. Build final cell set and subset matrices -----------------
  final_cells <- unique(c(ref_keep_barcodes,
                          ref_to_move_barcodes,
                          obs_keep_barcodes))
  
  infercnv_sub <- infercnv_obj
  infercnv_sub@expr.data  <- infercnv_obj@expr.data[, final_cells, drop = FALSE]
  
  # Also subset @count.data — required, otherwise plot_cnv() errors
  if (ncol(infercnv_obj@count.data) > 0) {
    infercnv_sub@count.data <- infercnv_obj@count.data[, final_cells, drop = FALSE]
  }
  
  # Clear subcluster info (refers to old indices, will break plotting)
  if (.hasSlot(infercnv_sub, "tumor_subclusters")) {
    infercnv_sub@tumor_subclusters <- NULL
  }
  
  new_cols <- colnames(infercnv_sub@expr.data)
  
  # ---- 4. Rebuild index lists --------------------------------------
  # Rebuild references = only cells in ref_keep_barcodes
  new_ref_indices <- list()
  for (grp in names(infercnv_obj@reference_grouped_cell_indices)) {
    grp_barcodes <- get_barcodes(
      infercnv_obj@reference_grouped_cell_indices, grp, all_cells
    )
    kept <- intersect(grp_barcodes, ref_keep_barcodes)
    if (length(kept) > 0) {
      new_ref_indices[[grp]] <- match(kept, new_cols)
    }
  }
  
  # Rebuild observations = subset of original obs + moved ref cells
  new_obs_indices <- list()
  
  # (a) Original obs groups, subsetted
  for (grp in names(infercnv_obj@observation_grouped_cell_indices)) {
    grp_barcodes <- get_barcodes(
      infercnv_obj@observation_grouped_cell_indices, grp, all_cells
    )
    kept <- intersect(grp_barcodes, obs_keep_barcodes)
    if (length(kept) > 0) {
      new_obs_indices[[grp]] <- match(kept, new_cols)
    }
  }
  
  # (b) Moved reference cells -> new obs group(s)
  if (length(ref_to_move_barcodes) > 0) {
    if (split_moved_by_group) {
      for (grp in names(infercnv_obj@reference_grouped_cell_indices)) {
        grp_barcodes <- get_barcodes(
          infercnv_obj@reference_grouped_cell_indices, grp, all_cells
        )
        moved <- intersect(grp_barcodes, ref_to_move_barcodes)
        if (length(moved) > 0) {
          new_obs_indices[[paste0("moved_", grp)]] <- match(moved, new_cols)
        }
      }
    } else {
      new_obs_indices[["moved_from_reference"]] <-
        match(ref_to_move_barcodes, new_cols)
    }
  }
  
  infercnv_sub@reference_grouped_cell_indices   <- new_ref_indices
  infercnv_sub@observation_grouped_cell_indices <- new_obs_indices
  infercnv_sub@count.data <- infercnv_obj@count.data[, final_cells, drop = FALSE]
  infercnv_sub@tumor_subclusters <- NULL
  # ---- 5. Sanity checks --------------------------------------------
  n_indexed <- length(unlist(new_ref_indices)) +
    length(unlist(new_obs_indices))
  if (n_indexed != ncol(infercnv_sub@expr.data)) {
    warning("Indexed cells (", n_indexed,
            ") != ncol(expr.data) (", ncol(infercnv_sub@expr.data),
            "). Some cells may be unassigned.")
  }
  
  all_idx <- unlist(c(new_ref_indices, new_obs_indices))
  if (length(all_idx) > 0 &&
      (max(all_idx) > ncol(infercnv_sub@expr.data) || min(all_idx) < 1)) {
    stop("Internal error: indices out of range after rebuild.")
  }
  
  # Warn about singleton groups (will break hclust in plot_cnv)
  singletons <- names(new_obs_indices)[lengths(new_obs_indices) == 1]
  if (length(singletons) > 0) {
    warning("Observation groups with only 1 cell (may break clustering): ",
            paste(singletons, collapse = ", "),
            ". Consider cluster_by_groups = FALSE in plot_cnv().")
  }
  
  # ---- 6. Report ---------------------------------------------------
  if (verbose) {
    message(sprintf("References kept : %d cells in %d group(s)",
                    length(unlist(new_ref_indices)),
                    length(new_ref_indices)))
    message(sprintf("Observations    : %d cells in %d group(s)",
                    length(unlist(new_obs_indices)),
                    length(new_obs_indices)))
    message(sprintf("  (moved from ref: %d)", length(ref_to_move_barcodes)))
    message(sprintf("Total cells     : %d", ncol(infercnv_sub@expr.data)))
  }
  infercnv::plot_cnv(
    infercnv_obj      = infercnv_sub,
    out_dir           = ".",
    output_filename   = "infercnv_ref_moved",
    output_format     = "pdf",
    cluster_by_groups = TRUE,    # cluster within each group, preserves labels
    x.center          = 1,
    x.range           = "auto",
    title             = "InferCNV — partial refs, others moved to obs",
    write_expr_matrix = FALSE)
  return(infercnv_sub)
}