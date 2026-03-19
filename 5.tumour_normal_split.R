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
library(Hmisc, lib.loc = "~/R_Library/4.5")
library(corrplot)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# Define file paths=========
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

srt_cell <- qread("spanorm_srt.qs")
DefaultAssay(srt_cell) <- "SpaNorm"
srt_cell <-  AddMetaData(srt_cell, readRDS("tumour_anno.Rds"))

srt_cell <- UpdateSeuratObject(srt_cell)

# tumour cells ----------
tumour_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno == "Tumour"])
if (!file.exists("tumour")){
  dir.create("tumour")
}
setwd("tumour")
tumour_srt <- do.spanorm(tumour_srt)
#tumour_srt <- tumour_srt %>% RunUMAP(dims = 1:20) %>% 
#  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
 tumour_srt <- tumour_srt %>% FindClusters(resolution = 0.8,algorithm = 4)
 tumour_srt <- CellCycleScoring(tumour_srt, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
 FeaturePlot(tumour_srt, c("G2M.Score", "S.Score"))
 
qs_save(tumour_srt, "tumour_srt.qs2")
spatial_plot(tumour_srt, outdir = "png/", name = "spanorm")
g <- DimPlot(tumour_srt, group.by= "category")+DimPlot(tumour_srt, group.by= "infercnv_hc", cols = "polychrome") 
g2 <-ImageDimPlot(tumour_srt, group.by= "category")+ImageDimPlot(tumour_srt, group.by= "infercnv_hc", cols = "polychrome") 
ggsave(plot = g/g2, "png/spanorm_category_subclone.png", width = 10, height = 10, dpi = 350, create.dir =T)

# Load Gene Signatures
clean_module <- readRDS("~/VisHD/clean_module.Rds")
names(clean_module) <- c("AR", "Inflammation", "NE1", "NE2", "Cycling", "Glycolysis") # Renaming tumour modules
source("~/VisHD/normal_markers.R") # Loads 'all_marker' (Normal cells)

tumour_srt <- AddModuleScore(tumour_srt, features = clean_module, name = "Module")
colnames(tumour_srt@meta.data)[grep("Module", fixed=T, colnames(tumour_srt@meta.data))] <- paste(names(clean_module), "Module")

g <- FeaturePlot(tumour_srt, paste(names(clean_module), "Module"), ncol = 3)&scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred")
ggsave(plot = g, "png/archetype_module_exp.pdf", width = 9 , height = 6)

SVGs <- readRDS("~/VisHD/LUT-245-07/tumour/SVGs.Rds")
f <- FeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% slice(1:20) %>% pull(symbol))
ggsave(plot = f, width = 15, height = 12, "png/SVG_Featureplot.png")
f <- ImageFeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% slice(1:20) %>% pull(symbol))
ggsave(plot = f, width = 15, height = 12, "png/SVG_ImageFeatureplot.png")

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

  
# normal cells ========
setwd(path)
if (!file.exists("normal")){
  dir.create("normal")
}
setwd("normal")

normal_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno != "Tumour"])
normal_srt <- do.spanorm(normal_srt)
normal_srt <- normal_srt %>% FindClusters(resolution = 0.8,algorithm = 4)
spatial_plot(normal_srt, outdir = "png/", name = "spanorm")
qs_save(normal_srt, "normal_srt.qs2")

normal_srt <- SeuratWrappers::RunBanksy(normal_srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
                                 assay = 'SpaNorm', slot = 'data',
                                 k_geom = c(15), assay_name = 'BANKSY_0.2')
normal_srt <- RunPCA(normal_srt, npcs = 30, features = rownames(normal_srt), reduction.name = "banksy0.2.pca")
normal_srt <- RunUMAP(normal_srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name  = "banksy0.2.umap")
normal_srt <-  FindNeighbors(normal_srt, reduction = "banksy0.2.pca", dims = 1:20)
normal_srt <-   FindClusters(normal_srt,resolution = 1, algorithm = 4)
qs_save(normal_srt, "normal_srt.qs2")
spatial_plot(normal_srt, "png/", "Bansky_lam0.2")

DefaultAssay(normal_srt) <- "SpaNorm"
DEG <- FindAllMarkers(normal_srt, test.used = "MAST")
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

normal_srt <- AddModuleScore(normal_srt, all_marker)
# Rename the resulting metadata columns (Cluster1, Cluster2...) to actual cell types
colnames(normal_srt@meta.data)[colnames(normal_srt@meta.data) %in% paste0("Cluster", 1:13)] <- names(all_marker)

g <- ImageFeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"))
ggsave("png/spanorm_ImageFeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)

g <- FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"))
ggsave("png/spanorm_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)