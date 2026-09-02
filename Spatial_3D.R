#' 3D spatial feature-expression plot
#'
#' Visualise spatial transcriptomics expression in 3D: x/y are cell
#' coordinates, z is feature expression. Handles a single feature or many
#' features on one plot, with three colouring modes, transparency, optional
#' downsampling, and a "flattened" static export.
#'
#' Colour logic (when `color.by = NULL`, auto-chosen):
#'   - >1 feature        -> colour by feature/signature NAME (one cloud each)
#'   - 1 feature + group -> colour by the discrete `group.by` variable
#'   - 1 feature only     -> colour by that feature's EXPRESSION (continuous)
#'
#' @param object      A Seurat object OR a data.frame.
#' @param features    Character vector: gene(s), module/signature score(s),
#'                    or metadata column(s). Anything FetchData() understands.
#' @param coords      How to get x/y. NULL = GetTissueCoordinates(); OR a
#'                    length-2 vector of metadata column names; OR a reduction
#'                    name (uses dims 1-2); OR a matrix/data.frame of coords.
#'                    For a plain data.frame `object`, give 2 column names.
#' @param assay       Optional assay to set as default before fetching.
#' @param layer       Expression layer/slot (Seurat v5 "layer", v4 "slot").
#' @param group.by    Optional metadata column for discrete colouring.
#' @param color.by    "expression" | "group" | "feature". NULL = auto.
#' @param palette_cont Named plotly colourscale for expression (e.g. "Viridis").
#' @param palette_disc Optional vector of colours for discrete modes.
#' @param opacity     Point transparency, 0-1.
#' @param point.size  Marker size.
#' @param z.scale     Multiplier on z (exaggerate/compress expression axis).
#' @param downsample  Optional integer: randomly keep this many cells (HD data
#'                    can be 100k+ points; 3D scatter gets heavy).
#' @param title       Plot title.
#' @param axis.labels Length-2 or 3 vector: x, y, (z) axis titles.
#' @param save        Optional path. ".html" -> interactive widget;
#'                    ".png"/".pdf"/".svg" -> flattened static image (kaleido).
#' @param scale       Resolution multiplier for static export.
#' @param width,height Pixel dimensions.
#' @param seed        RNG seed for downsampling reproducibility.
#'
#' @return A plotly object (also saved if `save` given).
spatial_expr_3d <- function(object,
                            features,
                            coords       = NULL,
                            assay        = NULL,
                            layer        = "data",
                            group.by     = NULL,
                            color.by     = NULL,
                            palette_cont = "Viridis",
                            palette_disc = NULL,
                            opacity      = 0.7,
                            point.size   = 3,
                            z.scale      = 1,
                            downsample   = NULL,
                            title        = NULL,
                            axis.labels  = c("x", "y", "expression"),
                            save         = NULL,
                            scale        = 2,
                            width        = 900,
                            height       = 750,
                            seed         = 1) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required: install.packages('plotly')")

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- 1. Collect coordinates + expression (+ group) into one data.frame ----
  df <- .collect_spatial_data(object, features, coords, assay, layer, group.by)
  # df columns: .x, .y, <one per feature>, and optionally .group

  # ---- 2. Optional downsample (HD data is large) ----------------------------
  if (!is.null(downsample) && downsample < nrow(df)) {
    set.seed(seed)
    df <- df[sample.int(nrow(df), downsample), , drop = FALSE]
    message("Downsampled to ", downsample, " cells.")
  }

  n_feat <- length(features)

  # ---- 3. Decide what colour encodes ----------------------------------------
  if (is.null(color.by)) {
    color.by <- if (n_feat > 1) "feature"
                else if (!is.null(group.by)) "group"
                else "expression"
  }
  color.by <- match.arg(color.by, c("expression", "group", "feature"))
  if (color.by == "group" && is.null(group.by))
    stop("color.by = 'group' needs `group.by`.")

  # discrete palette helper (colourblind-friendly Okabe-Ito, then HCL fallback)
  .disc_pal <- function(n, pal = NULL) {
    if (!is.null(pal)) return(rep(pal, length.out = n))
    okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
    if (n <= length(okabe)) okabe[seq_len(n)] else grDevices::hcl.colors(n, "Dynamic")
  }

  # ---- 4. Build the figure --------------------------------------------------
  fig <- plotly::plot_ly(width = width, height = height)

  if (color.by == "feature") {
    # One point cloud per feature, z = that feature's expression, coloured by name
    cols <- .disc_pal(n_feat, palette_disc)
    for (i in seq_len(n_feat)) {
      f <- features[i]
      fig <- plotly::add_trace(
        fig, type = "scatter3d", mode = "markers",
        x = df$.x, y = df$.y, z = df[[f]] * z.scale, name = f,
        marker = list(color = cols[i], size = point.size, opacity = opacity),
        hovertemplate = paste0(f, "<br>x:%{x:.1f}  y:%{y:.1f}  expr:%{z:.2f}<extra></extra>")
      )
    }
    z_lab <- "expression"

  } else if (color.by == "group") {
    # Single feature on z, one trace per group level (discrete colour + legend)
    f   <- features[1]
    grp <- factor(df$.group)
    lv  <- levels(grp)
    cols <- .disc_pal(length(lv), palette_disc)
    for (i in seq_along(lv)) {
      idx <- grp == lv[i]
      fig <- plotly::add_trace(
        fig, type = "scatter3d", mode = "markers",
        x = df$.x[idx], y = df$.y[idx], z = df[[f]][idx] * z.scale,
        name = lv[i],
        marker = list(color = cols[i], size = point.size, opacity = opacity),
        hovertemplate = paste0("%{text}<br>expr:%{z:.2f}<extra></extra>"),
        text = as.character(grp[idx])
      )
    }
    z_lab <- f

  } else { # "expression"
    f <- features[1]
    fig <- plotly::add_trace(
      fig, type = "scatter3d", mode = "markers",
      x = df$.x, y = df$.y, z = df[[f]] * z.scale,
      marker = list(
        color      = df[[f]],
        colorscale = palette_cont,
        showscale  = TRUE,
        size       = point.size,
        opacity    = opacity,
        colorbar   = list(title = f)
      ),
      hovertemplate = "x:%{x:.1f}  y:%{y:.1f}  expr:%{z:.2f}<extra></extra>"
    )
    z_lab <- f
  }

  fig <- plotly::layout(
    fig,
    title = title %||% (if (n_feat > 1) "Spatial expression (3D)" else features[1]),
    scene = list(
      xaxis = list(title = axis.labels[1]),
      yaxis = list(title = axis.labels[2]),
      zaxis = list(title = if (length(axis.labels) >= 3) axis.labels[3] else z_lab)
    )
  )

  # ---- 5. Save (flatten to static image, or interactive html) ---------------
  if (!is.null(save)) {
    reticulate::use_python("/scratch/pawsey1172/sweng/conda/envs/kaleido/bin/python", required = TRUE)
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


# ------------------------- internal helpers --------------------------------- #

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
  # explicit coordinate matrix / data.frame
  if (is.matrix(coords) || is.data.frame(coords)) return(as.matrix(coords)[, 1:2])

  # default: tissue coordinates (Seurat v5 returns x, y, cell -> take x, y)
  if (is.null(coords)) {
    xy <- Seurat::GetTissueCoordinates(object)
    num <- vapply(xy, is.numeric, logical(1))
    return(as.matrix(xy[, which(num)[1:2], drop = FALSE]))
  }

  # two metadata columns
  if (length(coords) == 2 && all(coords %in% colnames(object[[]])))
    return(as.matrix(object[[coords]]))

  # a reduction name -> first two dims
  if (length(coords) == 1 && coords %in% Seurat::Reductions(object))
    return(Seurat::Embeddings(object, coords)[, 1:2])

  stop("Could not resolve `coords`: give 2 metadata columns, a reduction name, ",
       "a coordinate matrix, or NULL for tissue coordinates.")
}

.fetch <- function(object, vars, layer) {
  # Legacy assays (e.g. SpaNorm, built via CreateAssayObject()) are old-style
  # `Assay` objects, not `Assay5`, and FetchData()'s `layer` arg doesn't
  # resolve against them cleanly. Convert to Assay5 so `layer` works; the old
  # v4 `slot` argument is defunct in SeuratObject 5.0.0, so no fallback there.
  assay_name <- Seurat::DefaultAssay(object)
  if (!inherits(object[[assay_name]], "Assay5"))
    object[[assay_name]] <- methods::as(object[[assay_name]], "Assay5")
  Seurat::FetchData(object, vars = vars, layer = layer)
}


# ------------------------------ examples ------------------------------------ #
if (FALSE) {
  
  # 1) Single gene, colour by continuous expression (the default)
  spatial_expr_3d(obj, features = "FOLH1", assay = "SpaNorm", opacity = 0.5, point.size = 2, save = "5.4.expression_proportion/3dfeature/folh1.html")

  # 2) Single feature, colour by a discrete group (e.g. curated tumour state)
  spatial_expr_3d(tumour_anno_srt, features = "FOLH1",
                  group.by = "tumour_state", opacity = 0.5)

  # 3) Multiple signatures on one plot, colour by signature name
  #    (module scores added via AddModuleScore, or genes directly)
  spatial_expr_3d(tumour_anno_srt,
                  features = c("Hypoxia1", "EMT1", "Stemness1"),
                  opacity  = 0.4, point.size = 2)

  # 4) Explicit coordinate columns + flatten to a static PNG
  spatial_expr_3d(tumour_anno_srt, features = "TACSTD2",
                  coords = c("x_centroid", "y_centroid"),
                  save   = "trop2_3d.png", scale = 3)

  # 5) Big HD sample: downsample for a responsive plot, save interactive html
  spatial_expr_3d(tumour_anno_srt, features = "AR",
                  downsample = 40000, save = "ar_3d.html")

  # 6) Plain data.frame input
  spatial_expr_3d(my_df, features = c("PSMA_score", "NE_score"),
                  coords = c("x", "y"))
}