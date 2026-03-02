library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)
library(mLLMCelltype)
library(pals)
library(ggplot2)
library(patchwork)

# install.packages("mLLMCelltype")

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# 4. Define file paths
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
path <- paste0(path, "/analysis/segmented_outputs/")
srt <-  qread(paste0(path, "banksy_srt.qs"))
setwd(path)
# SingleR ==========
# ref <- HumanPrimaryCellAtlasData()
# qsave(ref, "ref.qs")
# ref <- qread("~/VisHD/ref.qs")
# pred.main <- SingleR(test = GetAssayData(srt, assay = "Spatial", layer = "counts"), 
#                      ref = ref, 
#                      labels = ref$label.main)
# saveRDS(pred.main, paste0(path, "SingleR_pred.Rds"))

pred.main <- readRDS("SingleR_pred.Rds")
CB <- read.csv("CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
srt$category <- NA
srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
srt$category[srt$category == "CB NA"] <- "DT"

srt$SingleR.labels <- pred.main$labels
# DimPlot(srt, cols = "polychrome")
DefaultAssay(srt) <- "BANKSY_0.2"
srt <- FindNeighbors(srt, reduction = "banksy0.2.pca", dims = 1:20) %>%
  FindClusters(resolution = 2, algorithm = 4)
cell_col <- setNames( as.character(polychrome(36)), unique(pred.main$labels))
df <- data.frame(cluster = srt$BANKSY_0.2_snn_res.2, SR.pred = srt$SingleR.labels) %>% group_by(cluster, SR.pred) %>%
  summarise(n = n()) %>% group_by(cluster) %>% mutate(prop = n/sum(n))
f <- ggplot(df, aes(x = cluster, y = prop, fill = SR.pred))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values =cell_col)
g1 <- DimPlot(srt, cols = "polychrome",label = T, group.by = "BANKSY_0.2_snn_res.2", reduction = "banksy0.2.umap")+NoLegend()
g <- DimPlot(srt, group.by = "SingleR.labels", cols =cell_col,label = T, reduction = "banksy0.2.umap")+NoLegend()

ggsave("png/singleR_DimPlot.png", plot = f +g1+g, width = 15, height = 4, dpi = 350)
df2 <- data.frame(cluster = srt$BANKSY_0.2_snn_res.2,FOLH1 = srt[["Spatial"]]$counts["FOLH1", ]) 
labdf <- df2 %>% group_by(cluster) %>% summarise(mean = signif(mean(FOLH1), 3) , y = max(FOLH1)+1) 

p <- ggplot(df2, aes(x = cluster , y = FOLH1)) + 
  geom_jitter(alpha = 0.7)+
  geom_text(data = labdf, aes(x = cluster, label = mean, y = y), color = "red")+
  theme_minimal()
ggsave("png/cluster_FOLH1exp.png", plot = p, width = 15, height = 4, dpi = 350)

g <- ImageDimPlot(srt, group.by = "SingleR.labels", cols = cell_col)+NoLegend()
ggsave("png/singleR_ImageDimPlot.png", plot = g, width = 5, height = 4, dpi = 350)

# Found potential normal clsuters, prepare input for Infercnv
DefaultAssay(srt) <- "SpaNorm"

cellcut <- df %>% dplyr::filter(SR.pred == "Epithelial_cells"&prop<quantile(df$prop[df$SR.pred == "Epithelial_cells"], 0.1)) %>% pull(cluster) %>% as.character()
normal_clust <- labdf %>% dplyr::filter(mean <= signif(min(labdf$mean[labdf$cluster %in% cellcut]),1)) %>% pull(cluster) %>% as.character()

saveRDS(normal_clust, "normal_clust.Rds")

anno <- data.frame(colnames(srt), srt$BANKSY_0.2_snn_res.2)
write.table(anno, "anno.txt", sep = "\t", quote = F, row.names = F, col.names = F)

count.mat <- GetAssayData(srt, assay = "Spatial", layer = "counts")
qsave(count.mat, "countmat.qs")

# gene order files -----------------------
# gtf_data <- rtracklayer::import("~/SpaceRanger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")
# gtf_data <- gtf_data[gtf_data$type == "gene"]
# gene_ord <- data.frame(gene = gtf_data$gene_name, chr =as.character(seqnames(gtf_data)), start =start(gtf_data), end = end(gtf_data))
# saveRDS(gene_ord, "~/VisHD/gene_ord.Rds")
# write.table(gene_ord, "~/VisHD/gene_ord.txt", sep = "\t", col.names = F, row.names = F, quote = F)

# mllmcelltype ==========
# Set up cache directory to speed up processing
cache_dir <- "./mllmcelltype_cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)


deg <- FindAllMarkers(srt, test.used = "MAST") # Using MAST test for single-cell DE
deg <- deg %>% filter(p_val_adj < 0.05)
saveRDS(deg, "deg.res=2.Rds")

# top15_genes <- deg %>%
#   filter(p_val_adj < 0.05) %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 15)
# mat <- sapply(split(top15_genes, top15_genes$cluster), function(x) x$gene)
# rownames(mat) <- paste0("Gene", 1:15)
# write.table(mat, "~/VisHD/LUT-245-11/analysis/segmented_outputs/top15_deg.csv", row.names = T, col.names = T, quote = F, sep = "\t" )
# 
# 
# consensus_results <- interactive_consensus_annotation(
#   input = deg,
#   tissue_name = "Prostate primary tumour",  # provide tissue context
#   models = c(
#     "gemini-3-pro-preview",         
#     "gemini-2.5-pro"           # Google
#   ),
#   api_keys = list(
#     gemini = "AIzaSyCAofVWil8J4Fqat3Uge7hyogG-nx-tYUc"
#   ),
#   top_gene_count = 20,
#   controversy_threshold = 0.8,
#   entropy_threshold = 1.0,
#   cache_dir = cache_dir
# )
# 
# # Print structure of results to understand the data
# print("Available fields in consensus_results:")
# print(names(consensus_results))
# saveRDS(consensus_results, "gemini_consensus_anno.Rds")
# # Add annotations to Seurat object
# # Get cell type annotations from consensus_results$final_annotations
# cluster_to_celltype_map <- consensus_results$final_annotations
# 
# 
# # Create new cell type identifier column
# cell_types <- as.character(Idents(srt))
# for (cluster_id in names(cluster_to_celltype_map)) {
#   cell_types[cell_types == cluster_id] <- cluster_to_celltype_map[[cluster_id]]
# }
# 
# # Add cell type annotations to Seurat object
# srt$gemini_cell_type <- cell_types
# DimPlot(srt, group.by = "gemini_cell_type", cols = "polychrome", split.by = "category")
# ImageDimPlot(srt, group.by = "gemini_cell_type", cols = "polychrome")

# Run cell type annotation with a single LLM model
single_model_results <- annotate_cell_types(
  input = deg,
  tissue_name = "Prostate primary tumour",  # provide tissue context
  model = "gemini-3-pro-preview",  # specify a single model (Claude 4 Opus)
  api_key = "AIzaSyCAofVWil8J4Fqat3Uge7hyogG-nx-tYUc",  # provide the API key directly
  top_gene_count = 25
)
saveRDS(single_model_results, "gemini_3_anno.Rds")