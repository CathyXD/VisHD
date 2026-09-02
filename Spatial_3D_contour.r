#' Spatial expression contour — 2D, liftable to 3D
#'
#' Grids a spatial feature's expression once, then renders it either as a
#' flat 2D contour map (`to3d = FALSE`) or as the SAME field lifted into a
#' 3D surface with the contour lines projected onto the floor (`to3d = TRUE`).
#' Workflow: build the contour first, inspect it, then flip `to3d = TRUE`.
#'
#' Contours describe a single scalar field, so this takes ONE feature. For
#' several features use spatial_expr_surface_3d() (stacked surfaces) or
#' spatial_expr_3d() (scatter).
#'
#' @param object      Seurat object OR data.frame.
#' @param features    A single gene / signature score / metadata column.
#' @param coords      NULL = GetTissueCoordinates(); 2 metadata col names;
#'                    a reduction name; or a coord matrix. data.frame: 2 cols.
#' @param assay,layer Assay to activate / expression layer (v5) or slot (v4).
#' @param method      "bin" | "interp" | "loess" — how scattered points become
#'                    a regular grid (bin = fast/HD-friendly; interp = smooth
#'                    for irregular data; loess = statistically smoothed).
#' @param grid.n      Grid resolution per axis.
#' @param fun.aggregate Aggregation for method="bin".
#' @param to3d        FALSE = 2D contour map; TRUE = lift to 3D surface with
#'                    contours projected onto the z floor.
#' @param palette_cont plotly colourscale.
#' @param coloring    2D contour style: "fill" | "heatmap" | "lines".
#' @param ncontours   Approximate number of contour levels.
#' @param show.labels Print level values on the contour lines (2D).
#' @param opacity     Trace transparency.
#' @param z.scale     z multiplier (3D only).
#' @param downsample  Optional integer: keep N cells before gridding.
#' @param smooth.span loess span (method="loess").
#' @param title       Plot title.
#' @param axis.labels Length 2-3 vector: x, y, (z) axis titles.
#' @param aspectmode  "data" keeps x/y at true spatial scale (locks the 2D
#'                    aspect ratio too); else "auto"/"cube" (3D).
#' @param save        Path. ".html" -> interactive; image ext -> flattened
#'                    static (needs python kaleido).
#' @param scale       Resolution multiplier for static export.
#' @param width,height Pixel dimensions.
#' @param seed        RNG seed for downsampling.
#'
#' @return A plotly object (also saved if `save` given).
spatial_expr_contour <- function(object,
                                 features,
                                 coords        = NULL,
                                 assay         = NULL,
                                 layer         = "data",
                                 method        = c("bin", "interp", "loess"),
                                 grid.n        = 150,
                                 fun.aggregate = mean,
                                 to3d          = FALSE,
                                 palette_cont  = "Viridis",
                                 coloring      = c("fill", "heatmap", "lines"),
                                 ncontours     = 15,
                                 show.labels   = FALSE,
                                 opacity       = 1,
                                 z.scale       = 1,
                                 downsample    = NULL,
                                 smooth.span   = 0.1,
                                 title         = NULL,
                                 axis.labels   = c("x", "y", "expression"),
                                 aspectmode    = "data",
                                 save          = NULL,
                                 scale         = 2,
                                 width         = 900,
                                 height        = 750,
                                 seed          = 1) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required: install.packages('plotly')")
  method   <- match.arg(method)
  coloring <- match.arg(coloring)
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (length(features) != 1)
    stop("spatial_expr_contour() takes a single feature. For several features ",
         "use spatial_expr_surface_3d() or spatial_expr_3d().")
  f <- features

  # ---- 1. Collect + optionally downsample -----------------------------------
  df <- .collect_spatial_data(object, features, coords, assay, layer, group.by = NULL)
  if (!is.null(downsample) && downsample < nrow(df)) {
    set.seed(seed)
    df <- df[sample.int(nrow(df), downsample), , drop = FALSE]
    message("Downsampled to ", downsample, " cells.")
  }

  # ---- 2. Grid the field ONCE (shared by 2D and 3D) -------------------------
  g <- .grid_feature(df$.x, df$.y, df[[f]], method, grid.n, fun.aggregate, smooth.span)

  # ---- 3a. 2D contour -------------------------------------------------------
  if (!to3d) {
    fig <- plotly::plot_ly(
      x = g$x, y = g$y, z = g$z, type = "contour",
      colorscale = palette_cont, ncontours = ncontours,
      opacity = opacity, connectgaps = FALSE,
      contours = list(coloring = coloring, showlabels = show.labels),
      colorbar = list(title = f),
      width = width, height = height
    )
    yax <- list(title = axis.labels[2])
    if (identical(aspectmode, "data")) {          # lock true spatial aspect
      yax$scaleanchor <- "x"; yax$scaleratio <- 1
    }
    fig <- plotly::layout(fig, title = title %||% f,
                          xaxis = list(title = axis.labels[1]), yaxis = yax)

  # ---- 3b. Same field lifted to 3D: surface + projected contours ------------
  } else {
    fig <- plotly::plot_ly(width = width, height = height)
    fig <- plotly::add_surface(
      fig, x = g$x, y = g$y, z = g$z * z.scale,
      colorscale = palette_cont, opacity = opacity, connectgaps = FALSE,
      colorbar = list(title = f),
      contours = list(z = list(show = TRUE, usecolormap = TRUE,
                               highlightcolor = "#ffffff",
                               project = list(z = TRUE)))
    )
    fig <- plotly::layout(
      fig, title = title %||% f,
      scene = list(
        xaxis      = list(title = axis.labels[1]),
        yaxis      = list(title = axis.labels[2]),
        zaxis      = list(title = if (length(axis.labels) >= 3) axis.labels[3] else "expression"),
        aspectmode = aspectmode
      )
    )
  }

  # ---- 4. Save --------------------------------------------------------------
  if (!is.null(save)) {
    ext <- tolower(tools::file_ext(save))
    if (ext == "html") {
      if (!requireNamespace("htmlwidgets", quietly = TRUE))
        stop("Saving .html needs 'htmlwidgets'.")
      htmlwidgets::saveWidget(fig, save, selfcontained = TRUE)
    } else {
      ok <- requireNamespace("reticulate", quietly = TRUE) &&
            reticulate::py_module_available("kaleido")
      if (!ok)
        stop("Flattened static export needs the python 'kaleido' package.\n",
             "  In your env:  pip install kaleido\n",
             "  then: reticulate::use_python('<path>')\n",
             "  (or save as .html for an interactive file)")
      plotly::save_image(fig, file = save, scale = scale,
                         width = width, height = height)
    }
    message("Saved: ", save)
  }

  fig
}


# ---------------- gridding + shared helpers --------------------------------- #
# Identical to spatial_expr_surface_3d.R — if you source both files, keep just
# one copy of these four (.grid_feature, .collect_spatial_data, .get_coords,
# .fetch).

.grid_feature <- function(x, y, v, method, grid.n, fun.aggregate, smooth.span) {
  keep <- is.finite(x) & is.finite(y) & is.finite(v)
  x <- x[keep]; y <- y[keep]; v <- v[keep]
  rx <- range(x); ry <- range(y)

  if (method == "bin") {
    bx <- cut(x, breaks = seq(rx[1], rx[2], length.out = grid.n + 1),
              include.lowest = TRUE, labels = FALSE)
    by <- cut(y, breaks = seq(ry[1], ry[2], length.out = grid.n + 1),
              include.lowest = TRUE, labels = FALSE)
    z <- tapply(v, list(factor(by, levels = seq_len(grid.n)),
                        factor(bx, levels = seq_len(grid.n))), fun.aggregate)
    z <- matrix(as.vector(z), nrow = grid.n, ncol = grid.n)   # [y, x], gaps = NA
    list(x = seq(rx[1], rx[2], length.out = grid.n),
         y = seq(ry[1], ry[2], length.out = grid.n), z = z)

  } else if (method == "interp") {
    if (!requireNamespace("interp", quietly = TRUE, lib.loc = "~/R_Library/4.5"))
      stop("method='interp' needs the 'interp' package: install.packages('interp')")
    ip <- interp::interp(x, y, v, nx = grid.n, ny = grid.n, duplicate = "mean")
    list(x = ip$x, y = ip$y, z = t(ip$z))                     # -> [y, x]

  } else { # loess
    lo <- loess(v ~ x + y, span = smooth.span,
                control = loess.control(surface = "direct"))
    xc <- seq(rx[1], rx[2], length.out = grid.n)
    yc <- seq(ry[1], ry[2], length.out = grid.n)
    gg <- expand.grid(x = xc, y = yc)
    zp <- matrix(predict(lo, newdata = gg), nrow = length(xc), ncol = length(yc))
    list(x = xc, y = yc, z = t(zp))                           # -> [y, x]
  }
}

.collect_spatial_data <- function(object, features, coords, assay, layer, group.by) {
  if (inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat is required for Seurat objects.")
    if (!is.null(assay)) Seurat::DefaultAssay(object) <- assay
    xy   <- .get_coords(object, coords)
    expr <- .fetch(object, features, layer)
    missing <- setdiff(features, colnames(expr))
    if (length(missing))
      stop("Feature(s) not found: ", paste(missing, collapse = ", "))
    expr <- expr[, features, drop = FALSE]
    df <- data.frame(.x = xy[, 1], .y = xy[, 2], expr, check.names = FALSE)
    if (!is.null(group.by))
      df$.group <- Seurat::FetchData(object, vars = group.by)[, 1]
  } else if (is.data.frame(object)) {
    if (length(coords) != 2 || !all(coords %in% names(object)))
      stop("For a data.frame, `coords` must be 2 existing column names.")
    if (!all(features %in% names(object)))
      stop("Feature(s) not found: ",
           paste(setdiff(features, names(object)), collapse = ", "))
    df <- data.frame(.x = object[[coords[1]]], .y = object[[coords[2]]],
                     object[, features, drop = FALSE], check.names = FALSE)
    if (!is.null(group.by)) df$.group <- object[[group.by]]
  } else {
    stop("`object` must be a Seurat object or a data.frame.")
  }
  df
}

.get_coords <- function(object, coords) {
  if (is.matrix(coords) || is.data.frame(coords)) return(as.matrix(coords)[, 1:2])
  if (is.null(coords)) {
    xy  <- Seurat::GetTissueCoordinates(object)
    num <- vapply(xy, is.numeric, logical(1))
    return(as.matrix(xy[, which(num)[1:2], drop = FALSE]))
  }
  if (length(coords) == 2 && all(coords %in% colnames(object[[]])))
    return(as.matrix(object[[coords]]))
  if (length(coords) == 1 && coords %in% Seurat::Reductions(object))
    return(Seurat::Embeddings(object, coords)[, 1:2])
  stop("Could not resolve `coords`.")
}

.fetch <- function(object, vars, layer) {
  tryCatch(
    Seurat::FetchData(object, vars = vars, layer = layer),
    error = function(e) Seurat::FetchData(object, vars = vars, layer = layer)
  )
}


# ------------------------------ examples ------------------------------------ #
if (FALSE) {

  # --- Step 1: build the 2D contour first ---
  spatial_expr_contour(obj, features = "FOLH1",coloring = "heatmap",  assay = "SpaNorm", opacity = 0.5, to3d = TRUE, method = "interp", grid.n = 200,save = "5.4.expression_proportion/3dfeature/folh1.html")          # filled
  spatial_expr_contour(tumour_anno_srt, features = "FOLH1",
                       coloring = "lines", show.labels = TRUE)         # labelled lines
  spatial_expr_contour(tumour_anno_srt, features = "AR",
                       coloring = "heatmap", ncontours = 20)

  # --- Step 2: transform the SAME field into 3D ---
  spatial_expr_contour(tumour_anno_srt, features = "FOLH1", to3d = TRUE)

  # smooth field for irregular data (Xenium/segmented), then lift
  spatial_expr_contour(obj, features = "FOLH1",
                       method = "interp", grid.n = 200, to3d = TRUE,save = "5.4.expression_proportion/3dfeature/folh1.html")

  # flatten to static images
  spatial_expr_contour(tumour_anno_srt, features = "FOLH1",
                       save = "folh1_contour2d.png", scale = 3)
  spatial_expr_contour(tumour_anno_srt, features = "FOLH1", to3d = TRUE,
                       save = "folh1_contour3d.png", scale = 3)
}