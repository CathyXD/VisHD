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

path = "/home/sweng/VisHD/LUT-245-07"
setwd(path)
srt_cell <- qs_read("tumour_anno_srt.qs2")
DefaultAssay(srt_cell) <- "SpaNorm"

srt_cell <- UpdateSeuratObject(srt_cell)

setwd("normal")

normal_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno == "NoTumour"])
normal_srt <- do.spanorm(normal_srt)
normal_srt <- normal_srt %>% FindClusters(resolution = 0.8,algorithm = 4)
spatial_plot(normal_srt, outdir = "png/", name = "spanorm")
qs_save(normal_srt, "normal_srt.qs2")
cat("SpaNorm Done \n")

normal_srt <- SeuratWrappers::RunBanksy(normal_srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
                                        assay = 'SpaNorm', slot = 'data',
                                        k_geom = c(15), assay_name = 'BANKSY_0.2')
normal_srt <- RunPCA(normal_srt, npcs = 30, features = rownames(normal_srt), reduction.name = "banksy0.2.pca")
normal_srt <- RunUMAP(normal_srt, dims = 1:20, reduction = "banksy0.2.pca", reduction.name  = "banksy0.2.umap")
normal_srt <-  FindNeighbors(normal_srt, reduction = "banksy0.2.pca", dims = 1:20)
normal_srt <-   FindClusters(normal_srt,resolution = 1, algorithm = 4)
qs_save(normal_srt, "normal_srt.qs2")
spatial_plot(normal_srt, "png/", "Bansky_lam0.2")
cat("BANKSY DONE \n")

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

test <- FeaturePlot(normal_srt, SVEC_marker)
ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)


g <- FeaturePlot(normal_srt, names(all_marker), cols = c("white", "red"), reduction = )
ggsave("png/spanorm_FeaturePlot_normal_score.png", plot = g, width = 25, height = 15, dpi = 350)

test <- FeaturePlot(normal_srt, SVEC_marker)
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
