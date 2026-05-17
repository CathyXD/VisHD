library(Seurat)
library(qs, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")
library(dplyr)
library(SpaNorm, lib.loc = "~/R_Library/4.5")
library(qs2)
library(leidenbase, lib.loc = "~/R_Library/4.5")
library(UCell, lib.loc = "~/R_Library/4.5")
library(ggplot2)
library(patchwork)

library(msigdbr)
library(dplyr)
library(purrr)
library(stringr)


paths <- system("realpath ~/VisHD/LUT-*/normal/", intern = T)
source("~/VisHD/normal_markers.R")
SVEC_marker <-  c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4", "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")
# Select the specific sample directory based on the command line argument

## Try cell type annotation based on cluster using pre function====================
tme_markers <- list(
  Epithelial  = c("EPCAM","KRT8","KRT18","CDH1","MUC1"),
  T_Cells_Pan = c("CD3D","CD3E","CD3G","CD2"),
  CD8_T       = c("CD8A","CD8B","GZMK","GZMB"),
  CD4_T       = c("CD4","IL7R","LTB"),
  Tregs       = c("FOXP3","IL2RA","IKZF2","BATF"),
  B_Cells     = c("CD19","MS4A1","CD79A"),
  Plasma      = c("JCHAIN","MZB1","SDC1"),
  NK_Cells    = c("NCAM1","KLRB1","NKG7","GNLY"),
  Myeloid_Pan = c("CD14","LYZ","CD68","CSF1R"),
  TAMs        = c("CD163","MRC1","APOE","C1QA"),
  DCs         = c("HLA-DRA","CD1C","CLEC9A","THBD"),
  CAFs        = c("COL1A1","DCN","LUM","FAP","ACTA2"),
  Endothelial = c("PECAM1","VWF","ENG","CLDN5"),
  Pericytes   = c("RGS5","MCAM","CSPG4")
)


# ── Option B: secondary = named list (QC + fallback annotation) ────────────
secondary_markers <- list(
  Epithelial_broad = c("EPCAM","KRT5","KRT14","KRT17","KRT19","CLDN4","OCLN",
                       "TJP1","GJB2","GRHL2","ESRP1","ELF3"),
  Immune_broad     = c("PTPRC","HLA-A","HLA-B","HLA-C","B2M","LAPTM5",
                       "TYROBP","FCER1G","LST1","AIF1"),
  Stromal_broad    = c("VIM","FN1","THY1","S100A4","COL3A1","PDGFRA",
                       "PDGFRB","SPARC","POSTN","IGFBP3","IGFBP4"),
  Myeloid_broad    = c("CD14","ITGAM","ITGAX","CX3CR1","CSF1R","FCGR3A",
                       "S100A8","S100A9","VCAN","FCN1","LILRB2"),
  Lymphoid_broad   = c("CD3D","CD3E","CD19","MS4A1","NCAM1","IL2RA",
                       "CD27","CD38","SELL","CCR7","IL7R")
)

# Install required packages if you don't have them
# install.packages(c("msigdbr", "dplyr", "purrr", "stringr"))


# 1. Fetch the entire MSigDB C8 (Cell Type Signatures) collection for Human
cat("Fetching MSigDB C8 collection...\n")
c8_data <- msigdbr(species = "Homo sapiens", category = "C8")

# 2. Define our target cell types and their corresponding regex search terms
# Using regex allows us to catch variations (e.g., "T_CELL" vs "TCELL")
tme_search_terms <- list(
  Epithelial  = "EPITHELIAL|MALIGNANT",
  T_Cells_Pan = "T_CELL(?!.*(CD4|CD8|REGULATORY))", # T cells not explicitly CD4/CD8/Treg
  CD8_T       = "CD8.*T_CELL",
  CD4_T       = "CD4.*T_CELL",
  Tregs       = "REGULATORY_T_CELL|TREG",
  B_Cells     = "B_CELL",
  Plasma      = "PLASMA_CELL",
  NK_Cells    = "NK_CELL|NATURAL_KILLER",
  Myeloid_Pan = "MYELOID",
  TAMs        = "MACROPHAGE",
  DCs         = "DENDRITIC_CELL|DC",
  CAFs        = "FIBROBLAST",
  Endothelial = "ENDOTHELIAL",
  Pericytes   = "PERICYTE"
)

# 3. Create a helper function to extract and collapse genes for a given term
extract_c8_genes <- function(df, search_pattern) {
  df %>%
    # Filter for gene sets matching our regex pattern (case-insensitive)
    filter(str_detect(gs_name, regex(search_pattern, ignore_case = TRUE))) %>%
    # Pull the gene symbols
    pull(gene_symbol) %>%
    # Keep only unique genes to create a consensus list
    unique()
}

# 4. Map the function over our search terms to build the final list
c8_tme_markers <- map(tme_search_terms, ~extract_c8_genes(c8_data, .x))
source("~/VisHD/celltype_annotation_function.R")


metalist <- list()
for (path in paths){
  print(path)
  setwd(path)
  
  normal_srt <- qs_read("normal_srt.qs2")
  DefaultAssay(normal_srt) <- "SpaNorm"
  normal_srt <- AddModuleScore(normal_srt, all_marker)
  # Rename the resulting metadata columns (Cluster1, Cluster2...) to actual cell types
  colnames(normal_srt@meta.data)[colnames(normal_srt@meta.data) %in% paste0("Cluster", 1:13)] <- names(all_marker)
  
  g <- ImageFeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"))
  ggsave("png/spanorm_ImageFeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)
  g <- FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"))
  ggsave("png/spanorm_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)
  
  g <- FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"), reduction = "banksy0.2.umap")
  ggsave("png/Banksy_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)
  
  test <- FeaturePlot(normal_srt, SVEC_marker, reduction = "banksy0.2.umap")
  ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)
  
  normal_srt <- AddModuleScore(normal_srt, c8_tme_markers)
  colnames(normal_srt@meta.data)[colnames(normal_srt@meta.data) %in% paste0("Cluster", 1:14)] <- names(c8_tme_markers)
  g <- FeaturePlot(normal_srt, names(c8_tme_markers), cols = c("white", "red"), reduction = "banksy0.2.umap")
  ggsave("png/Banksy_FeaturePlot_c8_tme_markers.png", plot = g, width = 25, height = 15, dpi = 350)
  
  normal_srt <- tme_cluster_annotation_pipeline(
    obj                 = normal_srt,
    tme_markers         = tme_markers,
    secondary_genes     = c8_tme_markers,
    cluster_col         = "seurat_clusters",   # or e.g. "leiden_res0.5", "spatial_domain"
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
  
  g <- DimPlot(normal_srt, group.by= "celltype_annotation", cols = "polychrome", reduction= "banksy0.2.umap") +  DimPlot(normal_srt, group.by= "celltype_annotation", cols = "polychrome", reduction= "banksy0.2.umap", split.by = "category")
  ggsave(plot = g, "png/cell_type_anno_Dimplot.pdf", width = 12, height = 4)
  g <- ImageDimPlot(normal_srt, group.by= "celltype_annotation", cols = "polychrome")
  ggsave(plot = g, "png/cell_type_anno_ImageDimplot.png", width = 6, height = 4)
  g <- FeaturePlot(normal_srt, "secondary_expr_frac", reduction= "banksy0.2.umap") + VlnPlot(normal_srt, "nFeature_Spatial", group.by = "celltype_annotation")
  ggsave(plot = g, "png/cell_type_anno_QC.png", width = 6, height = 4)
  saveRDS(normal_srt@meta.data, "normal_meta.Rds")
  metalist[[strsplit(path, split = "/")[[1]][5]]] <- normal_srt@meta.data
}

df <- data.table::rbindlist(metalist, idcol = "sample", fill =T)
plotdf <- df %>% select(c("sample", "category", "celltype_annotation")) %>% group_by(sample, category, celltype_annotation) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))

ggplot(plotdf, aes(x = category, y = prop, fill = celltype_annotation)) +
  geom_bar(stat = "identity") +
  theme_classic()+
  facet_grid(.~ sample, space = "free", scale = "free")+
  scale_fill_manual(values = as.character(pals::polychrome(12)))
