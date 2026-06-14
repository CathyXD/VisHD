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

# ── Normal cells ───────────────────────────────────────────────────────────
setwd(normal_dir)
normal_srt <- subset(srt_cell_filtered, subset = tumour_anno == "Normal")
normal_srt <- do.spanorm(normal_srt)
normal_srt <- do.pearson_pca(normal_srt)
normal_srt <- FindClusters(normal_srt, resolution = 0.8, algorithm = 4)
spatial_plot(normal_srt, outdir = "png/", name = "spanorm")
qs_save(normal_srt, "normal_srt.qs2")
cat("SpaNorm Done\n")

normal_srt <- SeuratWrappers::RunBanksy(normal_srt, lambda = 0.2, verbose = TRUE,
                                        use_agf = TRUE, assay = "SpaNorm", slot = "data",
                                        k_geom = c(15), assay_name = "BANKSY_0.2")
normal_srt <- RunPCA(normal_srt, npcs = 30, features = rownames(normal_srt),
                     reduction.name = "banksy0.2.pca")
normal_srt <- RunUMAP(normal_srt, dims = 1:20, reduction = "banksy0.2.pca",
                      reduction.name = "banksy0.2.umap")
normal_srt <- FindNeighbors(normal_srt, reduction = "banksy0.2.pca", dims = 1:20)
normal_srt <- FindClusters(normal_srt, resolution = 1, algorithm = 4)
qs_save(normal_srt, "normal_srt.qs2")
spatial_plot(normal_srt, outdir = "png/", name = "Bansky_lam0.2")
cat("BANKSY DONE\n")

DefaultAssay(normal_srt) <- "SpaNorm"
DEG_n <- FindAllMarkers(normal_srt, test.use = "MAST")
saveRDS(DEG_n, "deg_spanorm.Rds")
run_gsea_panel(DEG_n, gene_sets, "deg_enrich.Rds")

normal_srt <- AddModuleScore(normal_srt, features = all_marker, name = "ct_")
added <- paste0("ct_", seq_along(all_marker))
colnames(normal_srt@meta.data)[match(added, colnames(normal_srt@meta.data))] <-
  names(all_marker)

ggsave("png/spanorm_ImageFeaturePlot_normal_score.png",
       plot = ImageFeaturePlot(normal_srt, names(all_marker), cols = c("white", "red")),
       width = 25, height = 15, dpi = 350)
ggsave("png/spanorm_FeaturePlot_normal_score.png",
       plot = FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red")),
       width = 25, height = 15, dpi = 350)
ggsave("png/pearson_FeaturePlot_normal_score.png",
       plot = FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"),
                          reduction = "pearsonumap"),
       width = 25, height = 15, dpi = 350)

# Gavish meta-programs (drop Malignant sheet for normal cells)
normal_srt <- score_meta_programs(
  normal_srt,
  meta_programs[setdiff(names(meta_programs), "Malignant")],
  out_dir = "png"
)

ggsave("png/SVEC_marker.png",
       plot = FeaturePlot(normal_srt, SVEC_marker, reduction = "banksy0.2.umap"),
       width = 12, height = 12)
ggsave("png/pearson_SVEC_marker.png",
       plot = FeaturePlot(normal_srt, SVEC_marker, reduction = "pearsonumap"),
       width = 12, height = 12)

plot_top_deg(normal_srt, DEG_n %>% filter(p_val_adj < 0.05),
             out_dir = "png", prefix = "spanorm_DEG")

# ── Cluster cell-type annotation ───────────────────────────────────────────
if (!"celltype_annotation" %in% colnames(normal_srt@meta.data)) {
  normal_srt <- tme_cluster_annotation_pipeline(
    obj                 = normal_srt,
    tme_markers         = tme_markers,
    secondary_genes     = c8_tme_markers,
    cluster_col         = "pearson_clusters",
    assay               = "SpaNorm",
    data_slot           = "data",
    expr_min_val        = 0,
    primary_expr_frac   = 0.05,
    secondary_expr_frac = 0.01,
    min_markers         = 3,
    conf_threshold      = 0.2,
    exclusivity_weight  = 0.30,
    detection_min       = 0.01,
    trim                = 0.10
  )
  qs_save(normal_srt, "normal_srt.qs2")
}

ggsave("png/cell_type_anno_Dimplot.pdf",
       plot = DimPlot(normal_srt, group.by = "celltype_annotation",
                      cols = "polychrome", reduction = "banksy0.2.umap"),
       width = 6, height = 4)
ggsave("png/cell_type_anno_ImageDimplot.png",
       plot = ImageDimPlot(normal_srt, group.by = "celltype_annotation",
                           cols = "polychrome"),
       width = 6, height = 4)
ggsave("png/cell_type_anno_QC.png",
       plot = FeaturePlot(normal_srt, "secondary_expr_frac",
                          reduction = "banksy0.2.umap") +
              VlnPlot(normal_srt, "nFeature_Spatial",
                      group.by = "celltype_annotation"),
       width = 6, height = 4)
ggsave("png/pearson_cell_type_anno_Dimplot.pdf",
       plot = DimPlot(normal_srt, group.by = "celltype_annotation",
                      cols = "polychrome", reduction = "pearsonumap"),
       width = 6, height = 4)

qs_save(normal_srt, "normal_srt.qs2")

tryCatch(srt2anndata(normal_srt, save_name = "normal"),
         error = function(e) message("srt2anndata(normal) failed: ", conditionMessage(e)))

cat("===================== normal done ====================\n")
