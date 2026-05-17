library(Seurat)
library(qs, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")
library(dplyr)
library(ggplot2)
library(patchwork)
library(SpaNorm, lib.loc = "~/R_Library/4.5")
library(qs2)



args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])
paths <- system("realpath ~/VisHD/LUT-*/normal", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
SVEC_marker <-  c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4", "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")
pre_tumour_marker <- c("AR", "FOLH1", "KLK2", "KLK3")
for(path in paths){
  setwd(path)
  setwd("tumour")
  
  tumour_srt <- qs_read("tumour_srt.qs2")
  g <- DimPlot(tumour_srt, group.by= "category")+DimPlot(tumour_srt, group.by= "infercnv_hc", cols = "polychrome") 
  g2 <-ImageDimPlot(tumour_srt, group.by= "category")+ImageDimPlot(tumour_srt, group.by= "infercnv_hc", cols = "polychrome") 
  ggsave(plot = g/g2, "png/spanorm_category_subclone.png", width = 10, height = 10, dpi = 350, create.dir =T)
  SVGs <- readRDS("~/VisHD/LUT-245-07/tumour/SVGs.Rds")
  f <- FeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% slice(1:20) %>% pull(symbol))
  ggsave(plot = f, width = 15, height = 12, "png/SVG_Featureplot.png")
  f <- ImageFeaturePlot(tumour_srt, as.data.frame(SVGs) %>% filter(svg.fdr <0.05) %>% arrange(desc(svg.F)) %>% filter(!grepl("MT-", symbol)) %>% slice(1:20) %>% pull(symbol))
  ggsave(plot = f, width = 15, height = 12, "png/SVG_ImageFeatureplot.png")
  
  #Seminal Vesicle Epithelial Cells
  test <- FeaturePlot(tumour_srt, SVEC_marker)
  ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)
  

  
  test <- FeaturePlot(tumour_srt, pre_tumour_marker)|VlnPlot(tumour_srt, pre_tumour_marker, ncol = 2)
  ggsave(plot = test, "png/tumour_marker.png", width = 10, height= 6)
  
  test <-VlnPlot(tumour_srt, pre_tumour_marker, group.by = "SR_Cluster", ncol = 2)
  ggsave(plot = test, "png/tumour_marker_SR_Cluster.png", width = 6, height= 6)
  
  deg_spanorm <- readRDS("~/VisHD/LUT-245-07/tumour/deg_spanorm.Rds")
  top5_genes <- deg_spanorm %>%
    filter(p_val_adj < 0.05) %>%
    filter(abs(pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5)
  f3 <- ImageFeaturePlot(tumour_srt, top5_genes$gene, size = 0.3)
  ggsave(plot = f3 + plot_layout(ncol = 5), paste0("png/", "spanorm_DEG_ImageFeaturePlot.png"), limitsize = FALSE, width = 15, height = 50, dpi = 350)
  
  f3 <- FeaturePlot(tumour_srt, top5_genes$gene)
  ggsave(plot = f3, paste0("png/", "spanorm_DEG_FeaturePlot.png"), limitsize = FALSE, width = 15, height = 30, dpi = 350)
  cat(path, "\n")
}

## Tumour further clearification: fix anno 

# DT_cluster <- which.max(table(tumour_srt$category, tumour_srt$seurat_clusters)["DT", ])
# C7_deg<- deg_spanorm %>%
#   filter(p_val_adj < 0.05) %>%
#   filter(abs(pct.1 - pct.2) > 0.1) %>%
#   filter(!grepl("MT-", gene)) %>%
#   filter(cluster == DT_cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 20)
# f <- FeaturePlot(tumour_srt, C7_deg$gene)
# ggsave(plot = f, paste0("png/", "DT_DEG_FeaturePlot_top20.png"), limitsize = FALSE, width = 15, height = 12, dpi = 350)
# deg_enrich <- readRDS("~/VisHD/LUT-245-07/tumour/deg_enrich.Rds")
# names(deg_enrich[["Hallmark"]])
# pathwayenrich_plot(top_n = 10, gsea_result = deg_enrich[["C5"]][["7"]], save.path = "png/DTvsCB_")

  for(path in paths){
    setwd(path)
    setwd("normal")
    
    normal_srt <- qs_read("normal_srt.qs2")
    DefaultAssay(normal_srt) <- "SpaNorm"
    test <- FeaturePlot(normal_srt, SVEC_marker)
    ggsave(plot = test, "png/SVEC_marker.png", width = 12, height= 12)
    
    test <- FeaturePlot(normal_srt, pre_tumour_marker)|VlnPlot(normal_srt, pre_tumour_marker, ncol = 2)
    ggsave(plot = test, "png/tumour_marker.png", width = 10, height= 6)
    
    test <-VlnPlot(normal_srt, pre_tumour_marker, group.by = "SR_Cluster", ncol = 2)
    ggsave(plot = test, "png/tumour_marker_SR_Cluster.png", width = 6, height= 6)
    cat(path, "\n")
  }