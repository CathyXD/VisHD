# apply_removal.R
# Drop contamination / junk subclusters flagged in removal_keep_map.csv
# from each per-compartment Seurat object, keeping only on-lineage clusters.
#
# Assumes: one Seurat object per (cell_type, method) whose subcluster IDs match
# the `cluster` column of the DEG CSVs, stored in a metadata column (edit
# `subcol` below to match yours, e.g. "seurat_clusters" or "sub.cluster").

library(Seurat)

map <- read.csv("removal_keep_map.csv", stringsAsFactors = FALSE)

# helper: parse "0,4,5" -> c(0,4,5)
parse_ids <- function(x) if (nchar(x) == 0) integer(0) else as.integer(strsplit(x, ",")[[1]])

# Example for one object -------------------------------------------------
# obj      : your Seurat object for, say, Stromal + rpca
# celltype : "Stromal"; meth : "rpca"; subcol : metadata column with subcluster id
clean_object <- function(obj, celltype, meth, subcol = "seurat_clusters",
                         drop_review = FALSE) {
  row  <- subset(map, cell_type == celltype & method == meth)
  keep <- parse_ids(row$keep_clusters)
  if (drop_review) {                      # by default REVIEW clusters are retained
    review <- parse_ids(row$review_clusters)
    keep   <- keep                        # (review already excluded from keep)
  } else {
    keep <- sort(unique(c(keep, parse_ids(row$review_clusters))))
  }
  ids <- as.integer(as.character(obj@meta.data[[subcol]]))
  obj$.keep <- ids %in% keep
  message(sprintf("%s/%s: keeping %d of %d cells across clusters {%s}",
                  celltype, meth, sum(obj$.keep), ncol(obj),
                  paste(keep, collapse = ",")))
  subset(obj, subset = .keep)
}

# usage:
# stromal_rpca_clean <- clean_object(stromal_rpca, "Stromal", "rpca",
#                                    subcol = "seurat_clusters")
