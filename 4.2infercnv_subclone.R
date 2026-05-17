# remotes::install_github("must-bioinfo/fastCNV", lib = "~/R_Library/4.5")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
library(qs, lib.loc = "~/R_Library/4.5")  
library(UCell, lib.loc = "~/R_Library/4.5")  
library(pals)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
cat("working on", path , "\n")
samples <- sapply(strsplit(paths, split = "/"), '[', 5)
names(paths) <- samples
i = samples[arg]


setwd(paths[[i]])
setwd("bined_ouput")
srt <- qs_read("subclone_srt.qs2")
# srt$subclone[!is.na(srt$subclone)] <- paste(srt$subclone[!is.na(srt$subclone)], "subclone")
# ref_cluster <-  unique(srt$ATAClone_cluster[is.na(srt$subclone)])
srt$subclone[is.na(srt$subclone)] <- "Normal"
gene_ord2 <- readRDS("~/VisHD/gene_ord2.Rds")
generef <- readRDS("~/VisHD/proberef.Rds")
generef <- rbind(generef[, 1:3], gene_ord2[setdiff(rownames(srt), rownames(generef)), ])
generef <- generef[order(generef$chr, generef$start), ]

# srt$ATAClone_cluster[srt$ATAClone_cluster %in% normal_clusters] <- "Normal"


library(infercnv)
infercnvobject = CreateInfercnvObject(raw_counts_matrix=as.matrix(GetAssayData(srt, assay = "Spatial.016um", layer = "counts")),
                                      annotations_file=as.data.frame(srt$subclone),
                                      delim="\t",
                                      min_max_counts_per_cell = c(50, Inf), 
                                      chr_exclude = c('chrY', 'chrM'), 
                                      gene_order_file= generef,
                                      ref_group_names="Normal")

if (!file.exists("infercnv_subclone")){
  dir.create("infercnv_subclone")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_subclone",
                               analysis_mode='cells',
                               cluster_by_groups= T,
                               cluster_references = T, 
                               denoise=T,
                               HMM=T,
                               save_rds = T, 
                               plot_steps = F, 
                               write_phylo = T, 
                               write_expr_matrix = F,
                               output_format = "pdf",
                               num_threads = 6, 
                               resume_mode = F)
infercnv::plot_per_group(infercnvobject, on_observations = F, out_dir ="infercnv_subclone",  write_expr_matrix = F, save_objects =T)