library(Seurat)
source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")
library(dplyr)
library(SpaNorm, lib.loc = "~/R_Library/4.5")
library(qs2)
library(leidenbase, lib.loc = "~/R_Library/4.5")
library(UCell, lib.loc = "~/R_Library/4.5")
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# Define file paths=========
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-245-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]
cat("working at", path, "\n")

srt_cell_filtered <- qs_read("tumour_anno_srt.qs2")


# tumour cells ----------
tumour_srt <- subset(srt_cell_filtered, cells = colnames(srt_cell_filtered)[srt_cell_filtered$tumour_anno == "Tumour"])
if (!file.exists("tumour")){
  dir.create("tumour")
}
setwd("tumour")
tumour_srt <- do.spanorm(tumour_srt)
#tumour_srt <- tumour_srt %>% RunUMAP(dims = 1:20) %>% 
#  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
 
 tumour_srt <- CellCycleScoring(tumour_srt, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
 FeaturePlot(tumour_srt, c("G2M.Score", "S.Score"))
 tumour_srt <- tumour_srt %>% FindClusters(resolution = 0.8,algorithm = 4)

spatial_plot(tumour_srt, outdir = "png/", name = "spanorm")

g <- DimPlot(tumour_srt, group.by= "category")+DimPlot(tumour_srt, group.by= "ATAClone_cluster", cols = "polychrome") 
g2 <-ImageDimPlot(tumour_srt, group.by= "category")+ImageDimPlot(tumour_srt, group.by= "ATAClone_cluster", cols = "polychrome") 
ggsave(plot = g/g2, "png/spanorm_category_subclone.png", width = 10, height = 10, dpi = 350, create.dir =T)

# Load Gene Signatures
clean_module <- readRDS("~/VisHD/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis") # Renaming tumour modules

tumour_srt <- AddModuleScore(tumour_srt, features = clean_module, name = "Module")
colnames(tumour_srt@meta.data)[grep("Module", fixed=T, colnames(tumour_srt@meta.data))] <- paste(names(clean_module), "Module")

g <- FeaturePlot(tumour_srt, paste(names(clean_module), "Module"), ncol = 3)&scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred")
ggsave(plot = g, "png/archetype_module_exp.pdf", width = 9 , height = 6)

cat("save object \n")
qs_save(tumour_srt, "tumour_srt.qs2")

tumour_srt <- SeuratWrappers::RunBanksy(tumour_srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
                                        assay = 'SpaNorm', slot = 'data',
                                        k_geom = c(15), assay_name = 'BANKSY_0.2')
tumour_srt <- RunPCA(tumour_srt, npcs = 30, features = rownames(tumour_srt), reduction.name = "banksy0.2.pca")
tumour_srt <- RunUMAP(tumour_srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name  = "banksy0.2.umap")
tumour_srt <-  FindNeighbors(tumour_srt, reduction = "banksy0.2.pca", dims = 1:20)
tumour_srt <-   FindClusters(tumour_srt,resolution = 1, algorithm = 4)
qs_save(tumour_srt, "tumour_srt.qs2")
spatial_plot(tumour_srt, "png/", "Bansky_lam0.2")
cat("BANKSY DONE \n")
srt2anndata(tumour_srt, save_name = 'tumour')

SVGs <- readRDS("~/VisHD/LUT-245-07/tumour/SVGs.Rds")
f <- FeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% dplyr::slice(1:20) %>% pull(symbol))
ggsave(plot = f, width = 15, height = 12, "png/SVG_Featureplot.png")
f <- ImageFeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% dplyr::slice(1:20) %>% pull(symbol))
ggsave(plot = f, width = 15, height = 12, "png/SVG_ImageFeatureplot.png")

 #Seminal Vesicle Epithelial Cells
SVEC_marker <-  c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4", "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")
test <- FeaturePlot(tumour_srt, SVEC_marker, reduction = "banksy0.2.umap")
  ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)

cat("===============Start DEG===================")
DEG <- FindAllMarkers(tumour_srt, test.used = "MAST")
saveRDS(DEG, "deg_spanorm.Rds")
  ## pathway enrichment
  Hall <- readRDS("~/VisHD/Hall.Rds")
  C6 <- readRDS("~/VisHD/C6.Rds")
  C5 <- readRDS("~/VisHD/C5.Rds")
  if (nrow(DEG) > 0) {
    DEG <- DEG %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
    
    if (nrow(DEG) > 0) {
      geneList <- lapply(split(DEG, DEG$cluster), function(x) setNames(x$avg_log2FC, x$gene))
      
      enrichlist <- lapply(list("Hallmark" = Hall, "C6" = C6, "C5" = C5), function(geneset) {
        enrich <- lapply(names(geneList), function(cluster_name) {
          x <- geneList[[cluster_name]]
          tryCatch({
            result <- clusterProfiler::GSEA(x, TERM2GENE = geneset)
            return(result)
          }, error = function(e) {
            message(sprintf("Skipping cluster '%s': %s", cluster_name, conditionMessage(e)))
            return(NULL)  # Return NULL for failed clusters
          }, warning = function(w) {
            message(sprintf("Warning in cluster '%s': %s", cluster_name, conditionMessage(w)))
            return(NULL)  # Return NULL for warned clusters, remove if warnings are acceptable
          })
        })
        names(enrich) <- names(geneList)
        return(enrich)
      })
      
      saveRDS(enrichlist, "deg_enrich.Rds")
      cat("ENRICHMENT DONE\n")
    }
  } else {
    cat("no DEG found\n")
  }

top5_genes <- DEG %>%
    filter(p_val_adj < 0.05) %>%
    filter(abs(pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 5)
  f3 <- ImageFeaturePlot(tumour_srt, top5_genes$gene, size = 0.3)
  ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)
  
  f3 <- FeaturePlot(tumour_srt, top5_genes$gene)
  ggsave(plot = f3, paste0("png/", "spanorm_DEG_FeaturePlot.png"), limitsize = FALSE, width = 15, height = 30, dpi = 350)
  cat("=====================tumour done====================", "\n")


# normal cells ========
setwd(path)
if (!file.exists("normal")){
  dir.create("normal")
}
setwd("normal")

tumour_srt <- subset(srt_cell_filtered, cells = colnames(srt_cell_filtered)[srt_cell_filtered$tumour_anno == "Normal"])
tumour_srt <- do.spanorm(tumour_srt)
tumour_srt <- tumour_srt %>% FindClusters(resolution = 0.8,algorithm = 4)
spatial_plot(tumour_srt, outdir = "png/", name = "spanorm")
qs_save(tumour_srt, "tumour_srt.qs2")
cat("SpaNorm Done \n")

tumour_srt <- SeuratWrappers::RunBanksy(tumour_srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
                                 assay = 'SpaNorm', slot = 'data',
                                 k_geom = c(15), assay_name = 'BANKSY_0.2')
tumour_srt <- RunPCA(tumour_srt, npcs = 30, features = rownames(tumour_srt), reduction.name = "banksy0.2.pca")
tumour_srt <- RunUMAP(tumour_srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name  = "banksy0.2.umap")
tumour_srt <-  FindNeighbors(tumour_srt, reduction = "banksy0.2.pca", dims = 1:20)
tumour_srt <-   FindClusters(tumour_srt,resolution = 1, algorithm = 4)
qs_save(tumour_srt, "tumour_srt.qs2")
spatial_plot(tumour_srt, "png/", "Bansky_lam0.2")
cat("BANKSY DONE \n")

DefaultAssay(tumour_srt) <- "SpaNorm"
DEG <- FindAllMarkers(tumour_srt, test.used = "MAST")
saveRDS(DEG, "deg_spanorm.Rds")
## pathway enrichment
Hall <- readRDS("~/VisHD/Hall.Rds")
C6 <- readRDS("~/VisHD/C6.Rds")
C5 <- readRDS("~/VisHD/C5.Rds")
if (nrow(DEG) > 0) {
  DEG <- DEG %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
  
  if (nrow(DEG) > 0) {
    geneList <- lapply(split(DEG, DEG$cluster), function(x) setNames(x$avg_log2FC, x$gene))
    
    enrichlist <- lapply(list("Hallmark" = Hall, "C6" = C6, "C5" = C5), function(geneset) {
      enrich <- lapply(names(geneList), function(cluster_name) {
        x <- geneList[[cluster_name]]
        tryCatch({
          result <- clusterProfiler::GSEA(x, TERM2GENE = geneset)
          return(result)
        }, error = function(e) {
          message(sprintf("Skipping cluster '%s': %s", cluster_name, conditionMessage(e)))
          return(NULL)  # Return NULL for failed clusters
        }, warning = function(w) {
          message(sprintf("Warning in cluster '%s': %s", cluster_name, conditionMessage(w)))
          return(NULL)  # Return NULL for warned clusters, remove if warnings are acceptable
        })
      })
      names(enrich) <- names(geneList)
      return(enrich)
    })
    
    saveRDS(enrichlist, "deg_enrich.Rds")
    cat("ENRICHMENT DONE\n")
  }
} else {
  cat("no DEG found\n")
}

tumour_srt <- AddModuleScore(tumour_srt, all_marker)
# Rename the resulting metadata columns (Cluster1, Cluster2...) to actual cell types
colnames(tumour_srt@meta.data)[colnames(tumour_srt@meta.data) %in% paste0("Cluster", 1:13)] <- names(all_marker)

g <- ImageFeaturePlot(tumour_srt, names(all_marker), cols = c("white", "red"))
ggsave("png/spanorm_ImageFeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)

g <- FeaturePlot(tumour_srt, names(all_marker), cols = c("white", "red"))
ggsave("png/spanorm_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)

test <- FeaturePlot(tumour_srt, SVEC_marker, reduction = "banksy0.2.umap")
ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)

top5_genes <- DEG %>%
    filter(p_val_adj < 0.05) %>%
    filter(abs(pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 5)
  f3 <- ImageFeaturePlot(tumour_srt, top5_genes$gene, size = 0.3)
  ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)
  
  f3 <- FeaturePlot(tumour_srt, top5_genes$gene)
  ggsave(plot = f3, paste0("png/", "spanorm_DEG_FeaturePlot.png"), limitsize = FALSE, width = 15, height = 30, dpi = 350)

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
source("~/VisHD/celltype_annotation_function.R")
tumour_srt <- tme_cluster_annotation_pipeline(
  obj                 = tumour_srt,
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

g <- DimPlot(tumour_srt, group.by= "celltype_annotation", cols = "polychrome", reduction= "banksy0.2.umap")
ggsave(plot = g, "png/cell_type_anno_Dimplot.pdf", width = 6, height = 4)
g <- ImageDimPlot(tumour_srt, group.by= "celltype_annotation", cols = "polychrome")
ggsave(plot = g, "png/cell_type_anno_ImageDimplot.png", width = 6, height = 4)
g <- FeaturePlot(tumour_srt, "secondary_expr_frac", reduction= "banksy0.2.umap") + VlnPlot(tumour_srt, "nFeature_Spatial", group.by = "celltype_annotation")
ggsave(plot = g, "png/cell_type_anno_QC.png", width = 6, height = 4)
qs_save(tumour_srt, "normal_srt.qs2")

srt2anndata(tumour_srt, save_name = 'normal')
