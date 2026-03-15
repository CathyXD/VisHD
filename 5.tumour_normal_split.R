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

srt_cell <- qread("banksy_srt.qs")
DefaultAssay(srt_cell) <- "Spatial"
srt_cell <-  AddMetaData(srt_cell, readRDS("tumour_anno.Rds"))

srt_cell <- UpdateSeuratObject(srt_cell)
tumour_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno == "Tumour"])

tumour_srt <- do.spanorm(tumour_srt)
qs_save(tumour_srt, "tumour_srt.qs2")

normal_srt <- subset(srt_cell, cells = colnames(srt_cell)[srt_cell$tumour_anno != "Tumour"])
normal_srt <- do.spanorm(normal_srt)
normal_srt <- do.bank