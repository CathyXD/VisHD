library(Seurat)
library(qs, "~/R_Library/4.5")
library(fastCNV)
source("~/VisHD/functions.R")
library(fastCNV)
library(dplyr)
library(enrichplot)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# Define file paths=========
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

srt_cell <- qread("banksy_srt.qs")
DefaultAssay(srt_cell) <- "SpaNorm"

srt_cell <- AddMetaData(srt_cell, readRDS("tumour_anno.Rds"))

# Tumour vs non Tumour differential expression analysis==========

DEG <- FindMarkers(srt_cell, ident.1 = "Tumour", ident.2 = "NoTumour", group.by = "tumour_anno", test.use = "MAST", min.diff.pct  = 0.5)
saveRDS(DEG, "deg_tumourvsnontumour.Rds")

# Hall <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
#   dplyr::select(gs_name, gene_symbol)
# 
# C6 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C6") %>%
#   dplyr::select(gs_name, gene_symbol)
# 
# C5 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5") %>%
#   dplyr::select(gs_name, gene_symbol)
# saveRDS(Hall, "~/VisHD/Hall.Rds")
# saveRDS(C6, "~/VisHD/C6.Rds")

## pathway enrichment ==========
Hall <- readRDS("~/VisHD/Hall.Rds")
C6 <- readRDS("~/VisHD/C6.Rds")
C5 <- readRDS("~/VisHD/C5.Rds")

geneList <- DEG %>% filter(p_val_adj <0.05) %>%  arrange(desc(avg_log2FC)) %>% select(avg_log2FC)
geneList <- setNames(geneList$avg_log2FC, rownames(geneList))
enrichlist <- lapply(list(Hall, C6, C5),function(geneset){
  enrich <- clusterProfiler::GSEA(geneList, TERM2GENE = geneset)
  ridgeplot(gsea_result, showCategory = 20, fill = "p.adjust") +
    scale_fill_gradient(low = "#d73027", high = "#4575b4") +
    ggtitle("GSEA Ridge Plot") +
    theme_bw(base_size = 12) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  pathwayenrich_plot(top_n = 8, gsea_result = enrich, save.path = "png/tumourvsnontumour_")
  return(enrich)
} )
saveRDS(enrichlist, "enrichlist_tumourvsnontumour.Rds")

# Tumour vs non Tumour differential expression analysis==========
deglist <- list()
enrichlist_list <- list()
for(cb in grep(" Tumour", fixed = T, value = T, grep("CB", value = T, unique(srt_cell$tumour_fullanno)))){
  DEG <- FindMarkers(srt_cell, ident.1 = "DT Tumour", ident.2 = cb, group.by = "tumour_fullanno", test.use = "MAST", min.diff.pct = 0.5)
  
  geneList <- DEG %>% filter(p_val_adj <0.05) %>%  arrange(desc(avg_log2FC)) %>% select(avg_log2FC)
  geneList <- setNames(geneList$avg_log2FC, rownames(geneList))
  enrichlist <- lapply(list(Hall, C6,C5),function(geneset){
    enrich <- clusterProfiler::GSEA(geneList, TERM2GENE = geneset)
    ridgeplot(gsea_result, showCategory = 20, fill = "p.adjust") +
      scale_fill_gradient(low = "#d73027", high = "#4575b4") +
      ggtitle("GSEA Ridge Plot") +
      theme_bw(base_size = 12) +
      theme(
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 9)
      )
    return(enrich)
    pathwayenrich_plot(top_n = 8, gsea_result = enrich, save.path = paste0("png/DTvs", cb))
  } )
  enrichlist_list[[paste0("DTvs", cb)]] <- enrichlist
  deglist[[paste0("DTvs", cb)]] <- DEG %>% filter(p_val_adj <0.05)
}

deglist <- data.table::rbindlist(deglist, fill = T, idcol = "comp")
saveRDS(deglist, "deg_DTvsCBtumour.Rds")
saveRDS(deglist, "deg_DTvsCBtumour.Rds")
