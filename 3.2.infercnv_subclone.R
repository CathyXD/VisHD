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

# Downsample Normal cells to avoid hclust timeout (~7h+ for 77k cells)
normal_cells <- colnames(srt)[srt$subclone == "Normal"]
max_normal   <- 5000
if (length(normal_cells) > max_normal) {
  set.seed(42)
  keep_normal <- sample(normal_cells, max_normal)
  srt <- srt[, c(colnames(srt)[srt$subclone != "Normal"], keep_normal)]
  message("Downsampled Normal cells from ", length(normal_cells), " to ", max_normal)
}

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

if (dir.exists("infercnv_subclone")) {
  unlink("infercnv_subclone", recursive = TRUE)
}
dir.create("infercnv_subclone")

infercnvobject = infercnv::run(infercnvobject,
                               cutoff=0.01, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir="infercnv_subclone",
                               analysis_mode='cells',
                               cluster_by_groups= T,
                               cluster_references = T, 
                               denoise=F,
                               HMM=F,
                               save_rds= TRUE,
                              plot_steps         = FALSE,
                              write_phylo        = TRUE,
                              write_expr_matrix  = FALSE,
                              no_plot            = TRUE,
                              num_threads        = 8,
                              resume_mode        = TRUE
                               )
plot_per_group(
  infercnvobject,
  on_references     = TRUE,
  on_observations   = TRUE,
  base_filename     = "infercnv_per_group",
  output_format     = "pdf",
  write_expr_matrix = FALSE,
  png_res           = 300,
  dynamic_resize    = 0,
  useRaster         = TRUE,
  out_dir           = "infercnv_subclone"
)