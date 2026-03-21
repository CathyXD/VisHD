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

# ── Option A: secondary = plain vector (QC only, no fallback annotation) ──
secondary_qc_panel <- c(
  "ACTB","GAPDH","B2M","HSP90AB1","RPL13","RPS18","MALAT1","NEAT1",
  "EEF1A1","TPT1","FTL","FTH1","HSPA1A","HSPA1B",
  # add broad cell-type housekeeping genes here ...
  VariableFeatures(visHD_obj)   # or simply pass HVGs
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

library(msigdbr)
library(dplyr)
library(purrr)
library(stringr)

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

# Print a summary of how many genes were extracted per cell type
cat("Number of unique genes extracted per cell type from C8:\n")
print(lengths(c8_tme_markers))

# ── Run the pipeline ───────────────────────────────────────────────────────
# Use existing seurat_clusters (default)
visHD_annotated <- tme_cluster_annotation_pipeline(
  obj                 = normal_srt,
  tme_markers         = tme_markers,
  secondary_genes     = c8_tme_markers,
  cluster_col         = "seurat_clusters",   # or e.g. "leiden_res0.5", "spatial_domain"
  assay               = "SpaNorm",
  data_slot           = "data",
  expr_min_val        = 0,
  primary_expr_frac   = 0.05,
  secondary_expr_frac = 0.10,
  min_markers         = 3,
  conf_threshold      = 0.2,
  exclusivity_weight  = 0.30,
  detection_min       = 0.01,
  trim                = 0.10
)

DimPlot(visHD_annotated, group.by= "celltype_annotation", cols = "polychrome", reduction= "banksy0.2.umap")
FeaturePlot(visHD_annotated, "secondary_expr_frac", reduction= "banksy0.2.umap")
VlnPlot(visHD_annotated, "nFeature_Spatial", group.by = "celltype_annotation")

# 
# DotPlot(visHD_annotated, features = lapply(tme_markers,function(x) intersect(x, rownames(visHD_annotated))), group.by = "celltype_annotation")
# # Marker dot plot
# plot_marker_dotplot(visHD_annotated, tme_markers)
# 
# # Which clusters fell back to secondary?
# subset(visHD_annotated@misc$cluster_annotation,
#        annotation_source %in% c("secondary", "unknown"))