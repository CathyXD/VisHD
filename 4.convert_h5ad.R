library(Seurat)
library(qs, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")
library(dplyr)
library(anndataR, lib.loc = "~/R_Library/4.5")

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])


# Define file paths=========
# Finds all directories matching 'LUT-*' to create a list of samples.
paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

srt_cell <- qread("banksy_srt.qs")
DefaultAssay(srt_cell) <- "Spatial"

slimsrt <- DietSeurat(srt_cell, layers = "counts", assays = "Spatial")
slimsrt <-  AddMetaData(slimsrt, readRDS("tumour_anno.Rds"))
adata <- as_AnnData(slimsrt)
write_h5ad(adata, "vishd_counts.h5ad", mode = "w")
