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
srt <- qs_read("srt_infercnv.qs2")
leiden_clusters <- readRDS("leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters

tumour_cluster <- list(c(4, 5, 7, 10), #07
                       c(4, 5, 7, 8),#09
                       c(3, 5), #10
                       c(1), #11
                       c(4, 6, 8, 9, 10), #15
                       c(3, 4),  #16, 
                       c(2, 6), #17
                       c(5, 9, 13, 10, 14)#20
                       
                       )
normal_cluster <- list(c(1,2,3, 6, 8), 
                       c(1,2,3, 6, 10), 
                       c(1,2,4, 6, 8),
                       c(2, 3), 
                       c(1,2, 3, 5, 7), 
                       c(1, 2), 
                       c(1, 3, 4, 5), 
                       c(1:4,6:8, 11)
                       )
normal_clusters <- normal_cluster[[arg]]
# include_clusters <- c(tumour_cluster[[arg]], normal_cluster[[arg]])
# srt <- subset(srt, cells = colnames(srt)[srt$ATAClone_cluster %in% include_clusters])
# srt$ATAClone_cluster <- factor(srt$ATAClone_cluster)
# srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[arg]], "Tumour", "Normal")
# s1 = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90"))
# s2 = ImageDimPlot(srt, split.by = "tumour_anno", group.by = "ATAClone_cluster", cols = "polychrome")
# s = s1 + s2 + plot_layout(widths = c(1, 2))
# ggsave("ATAClone_tumour_annotation.png", plot = s, width = 10, height = 3)
# qs_save(srt, "srt_ATAClone.qs2")
srt <- qs_read("srt_ATAClone.qs2")
gene_ord2 <- readRDS("~/VisHD/gene_ord2.Rds")
generef <- readRDS("~/VisHD/proberef.Rds")
generef <- rbind(generef[, 1:3], gene_ord2[setdiff(rownames(srt), rownames(generef)), ])
generef <- generef[order(generef$chr, generef$start), ]

# srt$ATAClone_cluster[srt$ATAClone_cluster %in% normal_clusters] <- "Normal"


library(infercnv)
infercnvobject = CreateInfercnvObject(raw_counts_matrix=as.matrix(GetAssayData(srt, assay = "Spatial.016um", layer = "counts")),
                                      annotations_file=as.data.frame(srt$ATAClone_cluster),
                                      delim="\t",
                                      min_max_counts_per_cell = c(50, Inf), 
                                      chr_exclude = c('chrY', 'chrM'), 
                                      gene_order_file= generef,
                                      ref_group_names=as.character(normal_clusters))

if (!file.exists("infercnv_ATAClone")){
  dir.create("infercnv_ATAClone")
}

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_ATAClone",
                               analysis_mode='cells',
                               cluster_by_groups= T,
                               cluster_references = T, 
                               denoise=F,
                               HMM=F,
                               save_rds = T, 
                               plot_steps = F, 
                               write_phylo = T, 
                               write_expr_matrix = F,
                               output_format = "pdf",
                               num_threads = 6, 
                               resume_mode = T)