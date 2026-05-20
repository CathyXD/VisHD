library(data.table)
library(qs)
library(Seurat)
library(readr)

# LUT-245-07
srt  <- qread("~/VisHD/LUT-245-07/analysis/segmented_outputs/spanorm_srt.qs")
CB <- read_csv("LUT-245-07/analysis/segmented_outputs/CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
srt$category <- NA
srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
srt$category[srt$category == "CB CB"] <- "CB 2"
srt$category[srt$category == "CB NA"] <- "DT"
ImageDimPlot(srt, group.by = "category")
qsave(srt, "~/VisHD/LUT-245-07/analysis/segmented_outputs/spanorm_srt.qs")

# LUT-245-11
srt  <- qread("~/VisHD/LUT-245-11/analysis/segmented_outputs/spanorm_srt.qs")
CB <- read_csv("LUT-245-11/analysis/segmented_outputs/CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
srt$category <- NA
srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
g <- ImageDimPlot(srt, group.by = "category")
CBcells <- CellSelector(g)
CBcells <- gsub("centroids_", "", g@data[CBcells, "cell"])
CB <- data.frame(Barcode = sprintf("cellid_%09d-1", as.numeric(CBcells)), CB = 1)
write.table(CB, "LUT-245-11/analysis/segmented_outputs/CB.csv", col.names = T, row.names = F, quote = F, sep ="\t")


srt$category <- NA
srt$category[which(as.character(srt$cell) %in% CBcells)] <- "CB 1"
srt$category[is.na(srt$category)] <- "DT"
ImageDimPlot(srt, group.by = "category")
qsave(srt, "~/VisHD/LUT-245-11/analysis/segmented_outputs/spanorm_srt.qs")


# LUT-245-15
srt  <- qread("~/VisHD/LUT-245-15/analysis/segmented_outputs/spanorm_srt.qs")
CB <- read_csv("LUT-245-15/analysis/segmented_outputs/CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
srt$category <- NA
srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
srt$category[srt$category == "CB NA"] <- "DT"

ImageDimPlot(srt, group.by = "category")
qsave(srt, "~/VisHD/LUT-245-15/analysis/segmented_outputs/spanorm_srt.qs")

# LUT-245-20
srt  <- qread("~/VisHD/LUT-245-20/analysis/segmented_outputs/spanorm_srt.qs")
CB <- read_csv("CB.csv")
CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
srt$category <- NA
srt$category <- paste("CB", CB$CB[match(srt$cell, CB$cellid)])
srt$category[srt$category == "CB NA"] <- "DT"

ImageDimPlot(srt, cells =  colnames(srt)[srt$nFeature_Spatial >50])

ImageDimPlot(srt, group.by = "category")
qsave(srt, "~/VisHD/LUT-245-20/analysis/segmented_outputs/spanorm_srt.qs")
