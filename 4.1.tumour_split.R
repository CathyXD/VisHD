library(Seurat)
library(dplyr)
library(SpaNorm,     lib.loc = "~/R_Library/4.5")
library(qs2)
library(leidenbase,  lib.loc = "~/R_Library/4.5")
library(UCell,       lib.loc = "~/R_Library/4.5")
library(ggplot2)
library(readxl)
library(purrr)
library(patchwork)
library(msigdbr)
library(stringr)

source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")
source("~/VisHD/celltype_annotation_function.R")

# ── Reference data ─────────────────────────────────────────────────────────
Hall <- readRDS("~/VisHD/Hall.Rds")
C6   <- readRDS("~/VisHD/C6.Rds")
C5   <- readRDS("~/VisHD/C5.Rds")
gene_sets <- list(Hallmark = Hall, C6 = C6, C5 = C5)

clean_module <- readRDS("~/VisHD/public_signature/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis")

# Gavish et al. meta-programs (sheet → program → genes)
meta_xlsx     <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname     <- excel_sheets(meta_xlsx)
meta_programs <- set_names(sheetname, sheetname) %>%
  map(~ read_excel(meta_xlsx, sheet = .x, col_names = TRUE)) %>%
  map(~ map(.x, ~ as.character(na.omit(.x))))

# Seminal Vesicle Epithelial Cell markers (used in both halves)
SVEC_marker <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                 "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")

# ── TME cell-type annotation panels ────────────────────────────────────────
# tme_markers: one vector per cell type — SVEC literature markers + union of
# all Gavish meta-programs within each cell-type sheet.
tme_markers <- c(
  list(SVEC = SVEC_marker),
  lapply(meta_programs, function(progs) unique(unlist(progs)))
)

# Map each tme_markers entry to a MSigDB C8 regex; unknown names fall back to
# the cell-type label uppercased so the pipeline still gets a non-empty hit set.
default_search <- function(nm) gsub("[^A-Z0-9]+", "_", toupper(nm))
known_search_terms <- list(
  SVEC          = "EPITHELIAL|SEMINAL_VESICLE",
  `B cells`     = "B_CELL",
  Endothelial   = "ENDOTHELIAL",
  Epithelial    = "EPITHELIAL|MALIGNANT",
  Fibroblasts   = "FIBROBLAST",
  Macrophages   = "MACROPHAGE",
  `CD4 T cells` = "CD4.*T_CELL",
  `CD8 T cells` = "CD8.*T_CELL"
)
tme_search_terms <- setNames(
  lapply(names(tme_markers), function(nm)
    if (!is.null(known_search_terms[[nm]])) known_search_terms[[nm]] else default_search(nm)),
  names(tme_markers)
)

# Cache MSigDB C8 (Cell Type Signatures) — network call on first run only.
c8_cache <- "~/VisHD/public_signature/c8_data_human.Rds"
c8_data <- if (file.exists(c8_cache)) {
  readRDS(c8_cache)
} else {
  cat("Fetching MSigDB C8 collection...\n")
  d <- msigdbr(species = "Homo sapiens", category = "C8")
  saveRDS(d, c8_cache)
  d
}

extract_c8_genes <- function(df, search_pattern) {
  df %>%
    filter(str_detect(gs_name, regex(search_pattern, ignore_case = TRUE))) %>%
    pull(gene_symbol) %>% unique()
}
c8_tme_markers <- map(tme_search_terms, ~ extract_c8_genes(c8_data, .x))

# ── CLI args ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || is.na(suppressWarnings(as.numeric(args[1])))) {
  stop("Usage: Rscript 4.1.tumour_normal_split.R <sample-index>")
}
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path  <- paths[arg]
i     <- basename(path)
cat("working at", path, "\n")

tumour_dir <- file.path(path, "tumour")
normal_dir <- file.path(path, "normal")
if (!dir.exists(tumour_dir)) dir.create(tumour_dir, recursive = TRUE)
if (!dir.exists(normal_dir)) dir.create(normal_dir, recursive = TRUE)

srt_cell_filtered <- qs_read(file.path(path, "tumour_subclone_srt.qs2"))


# ── Tumour cells ───────────────────────────────────────────────────────────
setwd(tumour_dir)
tumour_srt <- subset(srt_cell_filtered, subset = tumour_anno == "Tumour")
tumour_srt <- do.spanorm(tumour_srt)
tumour_srt <- do.pearson_pca(tumour_srt)
tumour_srt <- CellCycleScoring(tumour_srt,
                               s.features   = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes,
                               set.ident    = FALSE)
tumour_srt <- FindClusters(tumour_srt, resolution = 0.8, algorithm = 4)
spatial_plot(tumour_srt, outdir = "png/", name = "spanorm")

# ATAClone_cluster is absent from tumour_subclone_srt.qs2; prefer the refined
# `subclone` label and fall back to category so these diagnostics never crash.
clone_col <- if ("subclone" %in% colnames(tumour_srt@meta.data)) "subclone" else
             if ("ATAClone_cluster" %in% colnames(tumour_srt@meta.data)) "ATAClone_cluster" else "category"

g  <- DimPlot(tumour_srt, group.by = "category") +
      DimPlot(tumour_srt, group.by = clone_col, cols = "polychrome")
g2 <- ImageDimPlot(tumour_srt, group.by = "category") +
      ImageDimPlot(tumour_srt, group.by = clone_col, cols = "polychrome")
ggsave("png/spanorm_category_subclone.png", plot = g / g2,
       width = 10, height = 10, dpi = 350, create.dir = TRUE)

g_p <- DimPlot(tumour_srt, group.by = "category", reduction = "pearsonumap") +
       DimPlot(tumour_srt, group.by = clone_col, cols = "polychrome",
               reduction = "pearsonumap")
ggsave("png/pearson_category_subclone.png", plot = g_p,
       width = 10, height = 5, dpi = 350, create.dir = TRUE)

tumour_srt <- AddModuleScore(tumour_srt, features = clean_module, name = "Module")
colnames(tumour_srt@meta.data)[grep("Module", colnames(tumour_srt@meta.data), fixed = TRUE)] <-
  paste0(names(clean_module), "_Module")
mod_cols <- paste0(names(clean_module), "_Module")

g   <- FeaturePlot(tumour_srt, mod_cols, ncol = 3) &
       scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred")
ggsave("png/archetype_module_exp.pdf", plot = g, width = 9, height = 6)

g_p <- FeaturePlot(tumour_srt, mod_cols, ncol = 3, reduction = "pearsonumap") &
       scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred")
ggsave("png/pearson_archetype_module_exp.pdf", plot = g_p, width = 9, height = 6)

tumour_srt <- score_meta_programs(tumour_srt, meta_programs, out_dir = "png")

cat("save tumour object\n")
qs_save(tumour_srt, "tumour_srt.qs2")

tumour_srt <- SeuratWrappers::RunBanksy(tumour_srt, lambda = 0.2, verbose = TRUE,
                                        use_agf = TRUE, assay = "SpaNorm", slot = "data",
                                        k_geom = c(15), assay_name = "BANKSY_0.2")
tumour_srt <- RunPCA(tumour_srt, npcs = 30, features = rownames(tumour_srt),
                     reduction.name = "banksy0.2.pca")
tumour_srt <- RunUMAP(tumour_srt, dims = 1:20, reduction = "banksy0.2.pca",
                      reduction.name = "banksy0.2.umap")
tumour_srt <- FindNeighbors(tumour_srt, reduction = "banksy0.2.pca", dims = 1:20)
tumour_srt <- FindClusters(tumour_srt, resolution = 1, algorithm = 4)
qs_save(tumour_srt, "tumour_srt.qs2")
spatial_plot(tumour_srt, outdir = "png/", name = "Bansky_lam0.2")
cat("BANKSY DONE\n")

tryCatch(srt2anndata(tumour_srt, save_name = "tumour"),
         error = function(e) message("srt2anndata(tumour) failed: ", conditionMessage(e)))

if (file.exists("SVGs.Rds")) {
  SVGs <- readRDS("SVGs.Rds")
  top_svg <- as.data.frame(SVGs) %>%
    filter(svg.fdr < 0.05, !grepl("^MT-", symbol)) %>%
    arrange(desc(svg.F)) %>% dplyr::slice(1:20) %>% pull(symbol)
  if (length(top_svg) > 0) {
    ggsave("png/SVG_Featureplot.png",
           plot = FeaturePlot(tumour_srt, top_svg), width = 15, height = 12)
    ggsave("png/SVG_ImageFeatureplot.png",
           plot = ImageFeaturePlot(tumour_srt, top_svg), width = 15, height = 12)
  }
}

ggsave("png/SVEC_marker.png",
       plot = FeaturePlot(tumour_srt, SVEC_marker, reduction = "banksy0.2.umap"),
       width = 12, height = 12)
ggsave("png/pearson_SVEC_marker.png",
       plot = FeaturePlot(tumour_srt, SVEC_marker, reduction = "pearsonumap"),
       width = 12, height = 12)

cat("=============== Start DEG ===================\n")
DEG <- FindAllMarkers(tumour_srt, test.use = "MAST")
saveRDS(DEG, "deg_spanorm.Rds")
run_gsea_panel(DEG, gene_sets, "deg_enrich.Rds")
plot_top_deg(tumour_srt, DEG %>% filter(p_val_adj < 0.05),
             out_dir = "png", prefix = "spanorm_DEG")
cat("===================== tumour done ====================\n")


