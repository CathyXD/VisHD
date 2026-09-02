#' 3D surface plot of spatial feature expression
#'
#' Renders spatial expression as a 3D SURFACE: x/y are cell coordinates,
#' z is a gridded/interpolated expression field. Unlike the scatter version,
#' scattered per-cell values must be aggregated onto a regular grid first
#' (see `method`), so the surface is a smoothed field rather than raw cells.
#' Non-tissue regions become gaps (NA), which traces the tissue outline.
#'
#' Colour modes (auto when `color.by = NULL`):
#'   - 1 feature  -> "expression": colourscale over z, with colourbar
#'   - >1 feature -> "feature": one solid-colour surface per feature, stacked
#'                    with transparency, labelled by name in the legend
#' (Discrete `group.by` colouring is NOT meaningful for a surface — use the
#'  scatter function `spatial_expr_3d()` for that.)
#'
#' @param object      Seurat object OR data.frame.
#' @param features    Gene(s), signature/module score(s), or metadata col(s).
#' @param coords      NULL = GetTissueCoordinates(); 2 metadata col names;
#'                    a reduction name; or a coord matrix. For a data.frame,
#'                    give 2 column names.
#' @param assay,layer Assay to activate / expression layer (v5) or slot (v4).
#' @param method      "bin" | "interp" | "loess" (how scattered points become
#'                    a grid). See notes above.
#' @param grid.n      Grid resolution per axis (higher = finer, slower).
#' @param fun.aggregate Aggregation for method="bin" (mean, median, sum, max).
#' @param color.by    "expression" | "feature". NULL = auto.
#' @param palette_cont plotly colourscale for the expression surface.
#' @param palette_disc Optional colours for the per-feature surfaces.
#' @param opacity     Surface transparency (lower it when overlaying features).
#' @param contours    Logical; project contour lines onto the z axis.
#' @param z.scale     Multiplier on z (exaggerate/compress the expression axis).
#' @param downsample  Optional integer: randomly keep N cells before gridding
#'                    (recommended for interp/loess on large HD samples).
#' @param smooth.span loess span (method="loess"); smaller = more local detail.
#' @param title       Plot title.
#' @param axis.labels Length 2-3 vector: x, y, (z) axis titles.
#' @param aspectmode  plotly scene aspect: "auto", "cube", or "data"
#'                    ("data" keeps x/y at true spatial scale).
#' @param save        Path. ".html" -> interactive; ".png"/".pdf"/".svg" ->
#'                    flattened static image (needs python kaleido).
#' @param scale       Resolution multiplier for static export.
#' @param width,height Pixel dimensions.
#' @param seed        RNG seed for downsampling.
#'
#' @return A plotly object (also saved if `save` given).
spatial_expr_surface_3d <- function(object,
                                    features,
                                    coords        = NULL,
                                    assay         = NULL,
                                    layer         = "data",
                                    method        = c("bin", "interp", "loess"),
                                    grid.n        = 150,
                                    fun.aggregate = mean,
                                    color.by      = NULL,
                                    palette_cont  = "Viridis",
                                    palette_disc  = NULL,
                                    opacity       = 1,
                                    contours      = FALSE,
                                    z.scale       = 1,
                                    downsample    = NULL,
                                    smooth.span   = 0.1,
                                    title         = NULL,
                                    axis.labels   = c("x", "y", "expression"),
                                    aspectmode    = "auto",
                                    save          = NULL,
                                    scale         = 2,
                                    width         = 900,
                                    height        = 750,
                                    seed          = 1) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required: install.packages('plotly')")
  method <- match.arg(method)
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- 1. Collect coordinates + expression ----------------------------------
  df <- .collect_spatial_data(object, features, coords, assay, layer,
                              group.by = NULL)   # groups not used for surfaces

  if (!is.null(downsample) && downsample < nrow(df)) {
    set.seed(seed)
    df <- df[sample.int(nrow(df), downsample), , drop = FALSE]
    message("Downsampled to ", downsample, " cells.")
  }

  n_feat <- length(features)

  # ---- 2. Colour mode -------------------------------------------------------
  if (is.null(color.by)) color.by <- if (n_feat > 1) "feature" else "expression"
  color.by <- match.arg(color.by, c("expression", "feature"))

  .disc_pal <- function(n, pal = NULL) {
    if (!is.null(pal)) return(rep(pal, length.out = n))
    okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
    if (n <= length(okabe)) okabe[seq_len(n)] else grDevices::hcl.colors(n, "Dynamic")
  }

  # ---- 3. Build the figure --------------------------------------------------
  fig <- plotly::plot_ly(width = width, height = height)
  contour_spec <- list(z = list(show = isTRUE(contours), usecolormap = TRUE,
                                 project = list(z = isTRUE(contours))))

  if (color.by == "feature") {
    cols <- .disc_pal(n_feat, palette_disc)
    for (i in seq_len(n_feat)) {
      f <- features[i]
      g <- .grid_feature(df$.x, df$.y, df[[f]], method, grid.n,
                         fun.aggregate, smooth.span)
      # solid single-hue surface (constant colourscale), no per-surface bar
      fig <- plotly::add_surface(
        fig, x = g$x, y = g$y, z = g$z * z.scale,
        opacity = opacity, showscale = FALSE, name = f,
        colorscale = list(c(0, cols[i]), c(1, cols[i])),
        contours = contour_spec
      )
      # dummy marker just to get a labelled legend entry per feature
      fig <- plotly::add_trace(
        fig, type = "scatter3d", mode = "markers",
        x = NA, y = NA, z = NA, name = f, showlegend = TRUE,
        marker = list(color = cols[i], size = 6), hoverinfo = "skip"
      )
    }
    z_lab <- "expression"

  } else { # expression
    f <- features[1]
    g <- .grid_feature(df$.x, df$.y, df[[f]], method, grid.n,
                       fun.aggregate, smooth.span)
    fig <- plotly::add_surface(
      fig, x = g$x, y = g$y, z = g$z * z.scale,
      colorscale = palette_cont, opacity = opacity,
      colorbar = list(title = f), contours = contour_spec
    )
    z_lab <- f
  }

  fig <- plotly::layout(
    fig,
    title = title %||% (if (n_feat > 1) "Spatial expression surface (3D)" else features[1]),
    scene = list(
      xaxis      = list(title = axis.labels[1]),
      yaxis      = list(title = axis.labels[2]),
      zaxis      = list(title = if (length(axis.labels) >= 3) axis.labels[3] else z_lab),
      aspectmode = aspectmode
    )
  )

  # ---- 4. Save (flatten to static image, or interactive html) ---------------
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
             "  then point reticulate at it: reticulate::use_python('<path>')\n",
             "  (or save as .html for an interactive file instead)")
      plotly::save_image(fig, file = save, scale = scale,
                         width = width, height = height)
    }
    message("Saved: ", save)
  }

  fig
}


# ---------------- gridding: scattered points -> regular grid ---------------- #

.grid_feature <- function(x, y, v, method, grid.n, fun.aggregate, smooth.span) {
  keep <- is.finite(x) & is.finite(y) & is.finite(v)
  x <- x[keep]; y <- y[keep]; v <- v[keep]
  rx <- range(x); ry <- range(y)

  if (method == "bin") {
    bx <- cut(x, breaks = seq(rx[1], rx[2], length.out = grid.n + 1),
              include.lowest = TRUE, labels = FALSE)
    by <- cut(y, breaks = seq(ry[1], ry[2], length.out = grid.n + 1),
              include.lowest = TRUE, labels = FALSE)
    # z[y, x]: rows = y bins, cols = x bins; empty cells -> NA (gaps)
    z <- tapply(v, list(factor(by, levels = seq_len(grid.n)),
                        factor(bx, levels = seq_len(grid.n))), fun.aggregate)
    z <- matrix(as.vector(z), nrow = grid.n, ncol = grid.n)
    list(x = seq(rx[1], rx[2], length.out = grid.n),
         y = seq(ry[1], ry[2], length.out = grid.n),
         z = z)

  } else if (method == "interp") {
    if (!requireNamespace("interp", quietly = TRUE))
      stop("method='interp' needs the 'interp' package: install.packages('interp')")
    ip <- interp::interp(x, y, v, nx = grid.n, ny = grid.n, duplicate = "mean")
    # interp z is [x, y] -> transpose to [y, x] for plotly
    list(x = ip$x, y = ip$y, z = t(ip$z))

  } else { # loess
    lo <- loess(v ~ x + y, span = smooth.span,
                control = loess.control(surface = "direct"))
    xc <- seq(rx[1], rx[2], length.out = grid.n)
    yc <- seq(ry[1], ry[2], length.out = grid.n)
    gg <- expand.grid(x = xc, y = yc)           # x varies fastest
    zp <- matrix(predict(lo, newdata = gg),
                 nrow = length(xc), ncol = length(yc))  # [x, y]
    list(x = xc, y = yc, z = t(zp))             # -> [y, x]
  }
}


# --------- shared helpers (identical to spatial_expr_3d.R) ------------------ #
# If you source both files, keep only one copy of these three.

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
  stop("Could not resolve `coords`: give 2 metadata columns, a reduction name, ",
       "a coordinate matrix, or NULL for tissue coordinates.")
}

.fetch <- function(object, vars, layer) {
  tryCatch(
    Seurat::FetchData(object, vars = vars, layer = layer),
    error = function(e) Seurat::FetchData(object, vars = vars, layer = layer)
  )
}


# ------------------------------ examples ------------------------------------ #
if (FALSE) {

  # 1) Single gene, binned surface (default) — great for Visium HD
  spatial_expr_surface_3d(obj, features = "FOLH1", assay = "SpaNorm", opacity = 0.5, save = "5.4.expression_proportion/3dfeature/folh1.html")

  # 2) Smooth interpolated surface for sparse/irregular data (Xenium/segmented)
  spatial_expr_surface_3d(tumour_anno_srt, features = "AR",
                          method = "interp", grid.n = 200)

  # 3) loess-smoothed field (downsample first on big samples)
  spatial_expr_surface_3d(tumour_anno_srt, features = "TACSTD2",
                          method = "loess", smooth.span = 0.15,
                          downsample = 30000)

  # 4) Add projected contour lines + keep true spatial aspect
  spatial_expr_surface_3d(tumour_anno_srt, features = "FOLH1",
                          contours = TRUE, aspectmode = "data")

  # 5) Two signatures as stacked translucent surfaces, coloured by name
  spatial_expr_surface_3d(tumour_anno_srt,
                          features = c("Hypoxia1", "AR_program1"),
                          opacity = 0.7)

  # 6) Flatten to a static PNG
  spatial_expr_surface_3d(tumour_anno_srt, features = "FOLH1",
                          save = "folh1_surface.png", scale = 3)
}