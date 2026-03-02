library(qs, lib.loc = "~/R_Library/4.5")
library(Seurat)
library(dplyr)
library(ggplot2)
library(pals)
library(patchwork)
library(data.table)
library(scPearsonPCA, lib.loc =  "~/R_Library/4.5")

clean_module <- readRDS("~/CosMx/clean_module.Rds")
#clean_module <- lapply(clean_module, function(x) intersect(x, rownames(srt)))
names(clean_module) <- c("AR","Inflammation","NE1" , "NE2" , "Cycling" ,"Glycolysis"  )

source("~/CosMx/normal_markers.R")

paths <- list.files("~/VisHD", pattern = "LUT", full.names = T)
setwd(paths[1])
srt <- qread("banksy_srt.qs")


CB <- read.csv("CB.csv")
CB$cell <- as.character(as.numeric( sub(".*_(.*)-.*", "\\1", CB$Barcode) ))
srt$group <- "DT"
srt$group <- CB$CB[match(colnames(srt), CB$cell)]
srt$group[is.na(srt$group)] <- "DT"


tc <- colSums(GetAssayData(srt, assay = "Spatial", layer = "counts"))
genefreq <- scPearsonPCA::gene_frequency(GetAssayData(srt, assay = "Spatial", layer = "counts")) ## gene frequency (across all cells)
sum(genefreq)==1 # TRUE

srt <- Seurat::FindVariableFeatures(srt, nfeatures = 5000, assay ="Spatial")
hvgs <- VariableFeatures(srt, assay = "Spatial")

pcaobj <- 
  sparse_quasipoisson_pca_seurat(GetAssayData(srt, assay = "Spatial", layer = "counts")[hvgs,]
                                 ,totalcounts = tc
                                 ,grate = genefreq[hvgs]
                                 ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                                 ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                                 ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
  )
umapobj <- scPearsonPCA::make_umap(pcaobj)
srt[["pearsonpca"]] <- pcaobj$reduction.data
srt[["pearsonumap"]] <- umapobj$ump  ## umap

DimPlot(srt, reduction = "pearsonumap", split.by = "group") + DimPlot(srt, reduction = "banksy0.2.umap", split.by = "group")

FeaturePlot(srt, "FOLH1",reduction = "pearsonumap" )
p <- DimPlot(srt, group.by = "BANKSY_0.2_snn_res.1", reduction =  "banksy0.2.umap", cols = "polychrome", label = T) +
  coord_fixed()+NoLegend()
DimPlot(srt, group.by = "group", reduction =  "banksy0.2.umap") +
  ImageDimPlot(srt, group.by = "group")
ggsave("png/DTvsCB.png", width = 8, height = 4, dpi = 350)

DefaultAssay(srt) <- "SpaNorm"
srt <- RunUMAP(srt, metric = "correlation", reduction.name = "correlation.umap", dims = 1:20)
# cronanocal marker check--------
tumour_modules <- intersect(unlist(clean_module), rownames(srt))
normal_modules <- intersect(unlist(all_marker), rownames(srt))
srt <- AddModuleScore(srt, list(tumour_modules))
srt$tumour_score <- srt$Cluster1
srt <- AddModuleScore(srt, list(normal_modules))
srt$normal_score <- srt$Cluster1

FeaturePlot(srt, c("tumour_score","normal_score"), cols = c("white", "red"), reduction =  "banksy0.2.umap")
ggsave("png/Bansky0.2_FeaturePlot_TN_score.png", width = 10, height = 5, dpi = 350)

srt <- AddModuleScore(srt, all_marker)
colnames(srt@meta.data)[colnames(srt@meta.data) %in% paste0("Cluster", 1:13)] <- names(all_marker)
g <- FeaturePlot(srt, names(all_marker), cols = c("white", "red"), reduction =  "pearsonumap")&coord_fixed()
ggsave("png/Bansky0.2_FeaturePlot_normal_score.png", plot = g+p, width = 25, height = 15, dpi = 350)

srt <- AddModuleScore(srt, clean_module)
colnames(srt@meta.data)[colnames(srt@meta.data) %in% paste0("Cluster", 1:6)] <- paste(names(clean_module), "Module")
g <- FeaturePlot(srt, paste(names(clean_module), "Module"), cols = c("white", "red"), reduction =  "banksy0.2.umap", ncol = 3)&coord_fixed()
ggsave("png/Bansky0.2_FeaturePlot_tumour_score.png", plot = g/p+plot_layout(heights = c(2,1)), width = 10, height = 8, dpi = 350)



## prepare for Infercnv
annodf <- data.frame(cell = colnames(srt), cluster = srt$BANKSY_0.2_snn_res.1)
write.table(annodf, "anno.txt", quote = F, sep = "\t", col.names = F, row.names = F)

countmat <- GetAssayData(srt, assay = "Spatial", layer = "counts")

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl")
# Define which attributes you want to extract
attributes <- c("ensembl_gene_id", 
                "hgnc_symbol", 
                "chromosome_name", 
                "start_position", 
                "end_position", 
                "gene_biotype")

# Fetch the data
ref_hg38 <- getBM(attributes = attributes, 
                  mart = ensembl)

# Filter for standard chromosomes only (1-22, X, Y, MT)
ref_hg38 <- ref_hg38[ref_hg38$chromosome_name %in% c(1:22, "X", "Y"), ]
ref_hg38  <- ref_hg38[, 2:5] %>%filter(hgnc_symbol %in% rownames(srt)) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = c(1:22, "X", "Y"))) %>% 
  arrange(chromosome_name, start_position)

write.table(hg38,"~/CosMx/ref/hg38.txt", quote = F, sep = "\t", col.names = F, row.names = F)

countmat <- countmat[rownames(countmat) %in% hg38$`Gene name`,]
saveRDS(countmat, "countmat.Rds")
