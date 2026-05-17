#!/usr/bin/env Rscript
# plot_mean_variance.R
# Mean vs variance plot for genes in LUT-245-07/tumour/tumour_srt.qs2
# Uses SpaNorm log-normalised data; highlights SVGs and tumour marker genes.

suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(SpaNorm, lib.loc = "~/R_Library/4.5")
})

srt_path <- "~/VisHD/LUT-245-07/tumour/tumour_srt.qs2"
outdir   <- "~/VisHD/LUT-245-07/tumour/png"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Loading Seurat object ...")
srt <- qs_read(srt_path)

# ── Choose assay layer ─────────────────────────────────────────────────────────
# Prefer SpaNorm log-normalised; fall back to Spatial counts
if ("SpaNorm" %in% Assays(srt)) {
  mat <- LayerData(srt, assay = "SpaNorm", layer = "data")
  assay_used <- "SpaNorm (log-normalised)"
} else {
  mat <- LayerData(srt, assay = "Spatial", layer = "counts")
  assay_used <- "Spatial (raw counts)"
}
message("Using assay: ", assay_used)

# ── Compute per-gene mean and variance ────────────────────────────────────────
message("Computing gene-wise mean and variance (", nrow(mat), " genes × ", ncol(mat), " cells) ...")

if (inherits(mat, "sparseMatrix")) {
  gene_mean <- Matrix::rowMeans(mat)
  gene_var  <- sparseMatrixStats::rowVars(mat)   # fast sparse path
  if (is.null(gene_var)) {
    # fallback if sparseMatrixStats not available
    gene_var <- apply(mat, 1, var)
  }
} else {
  gene_mean <- rowMeans(mat)
  gene_var  <- apply(mat, 1, var)
}

df <- data.frame(
  gene = rownames(mat),
  mean = gene_mean,
  var  = gene_var,
  row.names = NULL
)

# ── Annotate SVGs and tumour markers ──────────────────────────────────────────
extract_svg_genes <- function(x) {
  if (is.character(x))        return(x)
  if (is.data.frame(x) || is.matrix(x)) {
    if ("spatially_variable" %in% colnames(x))
      return(rownames(x)[as.logical(x[["spatially_variable"]])])
    return(rownames(x))
  }
  character(0)
}

svg_genes <- character(0)
if ("SpaNorm" %in% Assays(srt)) {
  svg_meta <- srt[["SpaNorm"]]@misc[["SVGs"]]
  if (!is.null(svg_meta)) svg_genes <- extract_svg_genes(svg_meta)
}
if (length(svg_genes) == 0 && file.exists("~/VisHD/LUT-245-07/tumour/SVGs.Rds")) {
  svg_genes <- extract_svg_genes(readRDS("~/VisHD/LUT-245-07/tumour/SVGs.Rds"))
}

tumour_markers <- c("AR", "FOLH1", "KLK2", "KLK3", "KLK4",
                    "TMPRSS2", "NKX3-1", "HOXB13", "TRPM8")

df$group <- "Other"
df$group[df$gene %in% svg_genes]      <- "SVG"
df$group[df$gene %in% tumour_markers] <- "Tumour marker"
df$group <- factor(df$group, levels = c("Other", "SVG", "Tumour marker"))

df$label <- ifelse(df$gene %in% tumour_markers, df$gene, NA)

# ── Plot ───────────────────────────────────────────────────────────────────────
group_colors <- c("Other" = "grey70", "SVG" = "#4393C3", "Tumour marker" = "#D6604D")
group_alpha  <- c("Other" = 0.3,      "SVG" = 0.7,       "Tumour marker" = 1.0)
group_size   <- c("Other" = 0.4,      "SVG" = 0.6,       "Tumour marker" = 1.5)

p <- ggplot(df[order(df$group), ], aes(x = mean, y = var, colour = group,
                                        alpha = group, size = group)) +
  geom_point() +
  geom_text_repel(aes(label = label), colour = "#D6604D", size = 3,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = group_colors) +
  scale_alpha_manual(values = group_alpha) +
  scale_size_manual(values  = group_size) +
  scale_x_continuous(trans = "log1p", labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(trans = "log1p", labels = scales::label_number(accuracy = 0.01)) +
  labs(
    title    = "Gene mean vs variance — LUT-245-07 tumour",
    subtitle = paste0("Assay: ", assay_used, "  |  ", nrow(df), " genes"),
    x        = "Mean expression (log1p scale)",
    y        = "Variance (log1p scale)",
    colour   = NULL, alpha = NULL, size = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

out_png <- file.path(outdir, "mean_vs_variance.png")
ggsave(out_png, p, width = 8, height = 7, dpi = 200)
message("Saved: ", out_png)

# ── Save table ────────────────────────────────────────────────────────────────
write.csv(df[order(df$var, decreasing = TRUE), ],
          file.path(outdir, "gene_mean_variance.csv"), row.names = FALSE)
message("Done.")
