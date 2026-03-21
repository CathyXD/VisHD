library(Seurat)
library(qs, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")
library(dplyr)
library(SpaNorm, lib.loc = "~/R_Library/4.5")
library(qs2)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# Define file paths=========
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

srt_cell <- qread("spanorm_srt.qs")
DefaultAssay(srt_cell) <- "Spatial"
srt_cell <-  AddMetaData(srt_cell, readRDS("tumour_anno.Rds"))

srt_cell <- UpdateSeuratObject(srt_cell)
tumour_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno == "Tumour"])
if (!file.exists("tumour")){
  dir.create("tumour")
}
setwd("tumour")

tumour_srt <- do.spanorm(tumour_srt)
qs_save(tumour_srt, "tumour_srt.qs2")
spatial_plot(srt, outdir = "png/", name = "spanorm")

setwd(path)
if (!file.exists("normal")){
  dir.create("normal")
}
setwd("normal")

normal_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno != "Tumour"])
normal_srt <- do.spanorm(normal_srt)
normal_srt <- SeuratWrappers::RunBanksy(normal_srt, lambda = 0.2, verbose = TRUE, use_agf = TRUE,
                                 assay = 'SpaNorm', slot = 'data',
                                 k_geom = c(15), assay_name = 'BANKSY_0.2')