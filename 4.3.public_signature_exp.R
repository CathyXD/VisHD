library(Seurat)
library(patchwork)
library(ggplot2)
library(qs2)
library(readxl)
library(purrr)
library(scPearsonPCA, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")

# Read in gene sets ==========================
genesets <- readRDS("~/VisHD/public_signature/genesets2023.Rds")
sheetname <- excel_sheets("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx")
meta_programs <- set_names(sheetname, sheetname) |>
  map(~ read_excel("~/VisHD/public_signature/meta_programs_2025-01-29.xlsx", sheet = .x, col_names = TRUE)) |>
  # For each sheet (a tibble), split into a sublist keyed by column name,
  # dropping NA padding so each program is a clean character vector of genes.
  map(~ map(.x, ~ as.character(na.omit(.x))))
archetype_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(archetype_module) <- c("AR","Inflammation", "NE1","NE2", "Cycling","Glycolysis"  )

# Define file paths=========
args <- commandArgs(trailingOnly = TRUE)
arg  <- as.numeric(args[1])
paths <- system("realpath ~/VisHD/LUT-245-*/tumour", intern = T)
path <- paths[arg]
setwd(path)
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]
cat("working at", path, "\n")


out_qs <- "tumour_srt_with_public_signatures.qs2"
scores_present <- file.exists(out_qs)
if (scores_present) {
  cat("Loading existing object:", out_qs, "\n")
  srt <- qs_read(out_qs)
} else {
  srt <- qs_read("tumour_srt.qs2")
}

# ── Pearson residual PCA + UMAP ──────────────────────────────────────────────
if (!"pearsonumap" %in% Reductions(srt)) {
  srt <- do.pearson_pca(srt)
}

dir.create("png/public_signatures", showWarnings = FALSE, recursive = TRUE)

p <- DimPlot(srt, group.by = "category", cols = "polychrome") + DimPlot(srt,  group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25"))
p_pearson <- DimPlot(srt, group.by = "category", cols = "polychrome", reduction = "pearsonumap") + DimPlot(srt,  group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25"), reduction = "pearsonumap")
p2 <- ImageDimPlot(srt, group.by = "category", cols = "polychrome") + ImageDimPlot(srt,  group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25"))
ggsave(plot = p/p2, "png/public_signatures/category_tumouranno.png", width = 6, height = 8, dpi = 200)
ggsave(plot = p_pearson, "png/public_signatures/category_tumouranno_pearsonumap.png", width = 12, height = 4, dpi = 200)
# ── Output directories ────────────────────────────────────────────────────────
out_fp  <- file.path(path, "png", "public_signatures", "FeaturePlot")
out_ifp <- file.path(path, "png", "public_signatures", "ImageFeaturePlot")
dir.create(out_fp,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_ifp, showWarnings = FALSE, recursive = TRUE)

# ── Helper: score one parent group and save FeaturePlot + ImageFeaturePlot ───
# gene_sets : named list of character vectors (one vector per set)
# parent_name : string used as page title and filename prefix
plot_module_group <- function(srt, gene_sets, parent_name, out_fp, out_ifp, ncol = 3, add_scores = TRUE) {
  # Drop sets with fewer than 3 genes present in the object
  gene_sets <- Filter(function(g) length(intersect(g, rownames(srt))) >= 3, gene_sets)
  if (length(gene_sets) == 0) {
    message("Skipping '", parent_name, "': no sets with >=3 detectable genes")
    return(invisible(srt))
  }

  col_labels <- paste0(parent_name, " | ", names(gene_sets))

  if (add_scores) {
    # Unique internal prefix avoids column collisions across parent groups
    prefix <- paste0(".MS.", gsub("[^A-Za-z0-9]", ".", parent_name), ".")
    srt <- AddModuleScore(srt, features = gene_sets, name = prefix)

    # Rename numbered columns → "ParentName | SetName" for uniqueness
    raw_cols  <- paste0(prefix, seq_along(gene_sets))
    idx <- match(raw_cols, colnames(srt@meta.data))
    colnames(srt@meta.data)[idx] <- col_labels
  }

  # Paginate: max 20 panels per page
  max_per_page <- 20
  pages <- split(col_labels, ceiling(seq_along(col_labels) / max_per_page))

  for (pg in seq_along(pages)) {
    feats     <- pages[[pg]]
    titles    <- sub(paste0(parent_name, " | "), "", feats, fixed = TRUE)
    pg_suffix <- if (length(pages) > 1) paste0("_p", pg) else ""
    base_name <- gsub("[^A-Za-z0-9_-]", "_", parent_name)
    fname     <- paste0(base_name, pg_suffix, ".png")
    n_rows    <- ceiling(length(feats) / ncol)

    # FeaturePlot pages — one per reduction
    for (red in c("umap", "pearsonumap")) {
      if (!red %in% Reductions(srt)) next
      fp_list <- mapply(function(feat, ttl) {
        FeaturePlot(srt, feat, reduction = red) +
          scale_color_gradient2(low = "steelblue", mid = "gold", high = "indianred", name = ttl) +
          ggtitle(ttl) +
          theme(plot.title = element_text(size = 9), legend.position = "right")
      }, feats, titles, SIMPLIFY = FALSE)

      page_fp <- wrap_plots(fp_list, ncol = ncol) +
        plot_annotation(
          title = paste0(parent_name, " (", red, " ModuleScores)"),
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
      red_suffix <- if (red == "pearsonumap") "_pearsonumap" else ""
      fname_red  <- paste0(base_name, pg_suffix, red_suffix, ".png")
      ggsave(file.path(out_fp, fname_red), plot = page_fp,
             width = ncol * 4, height = n_rows * 3.5 + 0.8,
             dpi = 200, limitsize = FALSE)
    }

    # ImageFeaturePlot page
    ifp_list <- mapply(function(feat, ttl) {
      ImageFeaturePlot(srt, feat) +
        ggtitle(ttl) +
        labs(fill = ttl) +
        theme(plot.title = element_text(size = 9))
    }, feats, titles, SIMPLIFY = FALSE)

    page_ifp <- wrap_plots(ifp_list, ncol = ncol) +
      plot_annotation(
        title = parent_name,
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    ggsave(file.path(out_ifp, fname), plot = page_ifp,
           width = ncol * 4, height = n_rows * 4 + 0.8,
           dpi = 200, limitsize = FALSE)
  }

  message("Done: '", parent_name, "' (", length(gene_sets), " sets)")
  invisible(srt)
}

# ── 1. Archetype modules (flat named list → one page) ────────────────────────
srt <- plot_module_group(srt, archetype_module, "Archetype", out_fp, out_ifp, add_scores = !scores_present)

# ── 2. Meta-programs (nested: sheet → programs; one page per sheet) ───────────
for (sheet in names(meta_programs)) {
  srt <- plot_module_group(srt, meta_programs[[sheet]], sheet, out_fp, out_ifp, add_scores = !scores_present)
}

# ── 3. Public gene sets 2023 (auto-detect flat vs nested) ────────────────────
# If top-level values are character vectors → flat list, treat as one page.
# If top-level values are lists → nested; each top-level entry is one page.
is_nested <- is.list(genesets[[1]]) && !is.character(genesets[[1]])
if (is_nested) {
  for (parent in names(genesets)) {
    srt <- plot_module_group(srt, genesets[[parent]], parent, out_fp, out_ifp, add_scores = !scores_present)
  }
} else {
  srt <- plot_module_group(srt, genesets, "PublicSignatures2023", out_fp, out_ifp, add_scores = !scores_present)
}

message("\nAll module score plots saved to:\n  ", out_fp, "\n  ", out_ifp)

qs_save(srt, "tumour_anno_srt_with_public_signatures.qs2")

# ── 4. scMetabolism KEGG programs ────────────────────
# install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
# devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.0
# remotes::install_version("loe", "1.1", lib = "~/R_Library/4.5")
# install.packages(c('loe', "wesanderson"), lib = "~/R_Library/4.5")
# devtools::install_github("YosefLab/VISION@v2.1.0", lib = "~/R_Library/4.5") #Please note that the version would be v2.1.0
# devtools::install_github("wu-yc/scMetabolism", lib = "~/R_Library/4.5")

library(scMetabolism, lib.loc = "~/R_Library/4.5")
library(AUCell, lib.loc = "~/R_Library/4.5")
library(GSEABase, lib.loc = "~/R_Library/4.5")
library(wesanderson, lib.loc = "~/R_Library/4.5")

# scMetabolism hardcodes assay "RNA" — alias the raw Spatial counts so
# sc.metabolism.Seurat() can find them.
if (!"RNA" %in% Assays(srt)) {
  srt[["RNA"]] <- CreateAssayObject(counts = GetAssayData(srt, assay = "Spatial", layer = "counts"))
}

countexp.Seurat<-sc.metabolism.Seurat(obj = srt, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
colnames(metabolism.matrix) <- colnames(countexp.Seurat)

pathways <- rownames(metabolism.matrix)

qs_save(countexp.Seurat, "metabolism_srt.qs2")
write.csv(metabolism.matrix, "metabolism_KEGG_score.csv")

out_met <- file.path(path, "png", "public_signatures", "Metabolism")
dir.create(out_met, showWarnings = FALSE, recursive = TRUE)

met_ncol  <- 5
met_pages <- split(pathways, ceiling(seq_along(pathways) / 20))

DefaultAssay(countexp.Seurat) <- "RNA"

for (pg in seq_along(met_pages)) {
  feats  <- met_pages[[pg]]
  p_list <- lapply(feats, function(pw) {
    srt$score <- as.numeric(metabolism.matrix[pw, ])
    FeaturePlot(obj = srt, features = "score", reduction = "pearsonumap", order = TRUE) +
          scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred", name = "AUCell", midpoint = 0.1) +
          ggtitle(pw) +
          theme(plot.title = element_text(size = 9), legend.position = "right")
  })
  page <- wrap_plots(p_list, ncol = met_ncol) +
    plot_annotation(
      title = "KEGG metabolism (AUCell)",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  ggsave(file.path(out_met, paste0("metabolism_p", pg, ".png")), plot = page,
         width = met_ncol * 3, height = ceiling(length(feats) / met_ncol) * 3 + 0.8,
         dpi = 200, limitsize = FALSE)
}

message("Done: KEGG metabolism pathway plots saved to ", out_met)
