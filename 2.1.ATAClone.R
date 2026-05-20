library("ATAClone")
library("MatrixGenerics")
library("ggplot2")
library(qs2)
library(dplyr)
library(Matrix)
library(Seurat)

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


# #make Granges for genes
# #this is not in same order as count matrix for some reason... but has gene coordinates
# vis_probes <- data.table::fread("https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.bed")
# #this is in same order as count matrix... but does not have gene coordinates
# vis_probes2 <- data.table::fread("https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv")
# #same dimensions as count matrix
# length(unique(sapply(strsplit(as.character(vis_probes2$probe_id), "\\|"), `[`, 2)[vis_probes2$included])) == nrow(srt)
# 
# #get included probes first
# included_hashes <- sapply(strsplit(as.character(vis_probes2$probe_id), "\\|"), `[`, 3)[vis_probes2$included]
# 
# vis_genes <- sapply(strsplit(as.character(vis_probes$V4), "\\|"), `[`, 2)
# vis_hash <- sapply(strsplit(as.character(vis_probes$V4), "\\|"), `[`, 3)
# names(vis_genes) <- vis_hash
# vis_genes <- vis_genes[included_hashes]
# chr_split <- split(setNames(vis_probes$V1, vis_hash)[included_hashes], vis_genes)
# start_split <- split(setNames(vis_probes$V2, vis_hash)[included_hashes], vis_genes)
# end_split <- split(setNames(vis_probes$V3, vis_hash)[included_hashes], vis_genes)
# 
# gene_coords <- setNames(paste0(sapply(chr_split, `[`, 1), ":", sapply(start_split, min), "-", sapply(end_split, max)), names(chr_split))
# gene_coords <- gene_coords[rownames(srt)[rownames(srt) %in% names(chr_split)]]
# gene_gr <- GenomicRanges::GRanges(gene_coords)
# 
# generef <- as.data.frame(gene_gr)
# rownames(generef) <- names(gene_gr)
# colnames(generef) <- c("chr", "start", "end", "width", "strand")
# saveRDS(generef, "~/VisHD/proberef.Rds")
generef <- readRDS("~/VisHD/proberef.Rds")

bin_width <- 10e6 
gene_coords <- generef %>% 
  filter(!chr %in% c("chrM", "chrY")) %>%
  mutate(bin = paste0(chr, "_", floor(start / bin_width) * bin_width)
  )
gene_coords$gene <- rownames(gene_coords)
gene_bin_map <- gene_coords %>%
  filter(gene %in% rownames(srt)) %>%
  select(gene , bin) %>%
  distinct(gene, .keep_all = TRUE)

gene_bin_map <- gene_bin_map %>% filter(bin %in% names(which(table(gene_bin_map$bin)>5 ))) # filter bins with less than 5 genes 
mat_sub <- GetAssayData(srt, layer = "counts")[gene_bin_map$gene, ]

bin_matrix <- function(mat, gene_bin_map){
  # Build a bin aggregation matrix (genes x bins)
  bins       <- unique(gene_bin_map$bin)
  gene_order <- gene_bin_map$gene
  
  agg_mat <- Matrix(0, nrow = length(bins), ncol = length(gene_order), sparse = TRUE)
  rownames(agg_mat) <- bins
  colnames(agg_mat) <- gene_order
  
  for (b in bins) {
    genes_in_bin <- gene_bin_map$gene[gene_bin_map$bin == b]
    agg_mat[b, genes_in_bin] <- 1
  }
  
  # bin x cell = (bin x gene) %*% (gene x cell)
  bin_mat <- agg_mat %*% mat[gene_order, ]
  # Result: sparse bin x cell matrix
  return(bin_mat)
}
bin_mat <- bin_matrix(mat_sub, gene_bin_map)
stable_counts_filtered <- bin_mat[, apply(bin_mat, 2, function(x) sum(x>0))>10] 

### Sample 10000 cells from Breast Cancer VisHD data set as external normal reference ==============
# library(fastCNVdata)
# library(fastCNV)
# HDBreast <- load_HDBreast()
# HDBreast <- annotations8umTo16um(HDBreast, referenceVar = "annots_8um")
# refmat <- GetAssayData(HDBreast, layer = "counts", assay = "Spatial.016um")
# gene_bin_map_ref <- gene_bin_map %>% filter(gene %in% rownames(refmat))
# temp <- names(which(HDBreast$projected_annots_8um == "NoTumor"))
# temp <- temp[temp %in% colnames(refmat)]
# stable_counts_ref <- refmat[gene_bin_map_ref$gene, sample(temp, 10000, replace = F) ]
# stable_counts_ref <- bin_matrix(stable_counts_ref, gene_bin_map_ref)
# colnames(stable_counts_ref) <- paste0("HDB_", colnames(stable_counts_ref))
# saveRDS(stable_counts_ref, "~/VisHD/stable_counts_ref.Rds")
stable_counts_ref <- readRDS("~/VisHD/stable_counts_ref.Rds")

#0.01 by default. Usually set in range of 0.005-0.03. 0 = Poisson distributed data.
nb_overdispersion <- 0.005
#pseudo-count for vst or power transformation
pseudo_count <- 0.5
#0.4 by default. Usually set in range of 0.3 (high overdispersion) to 0.5 (low overdispersion). Setting higher than optimal can help with removing
#technical variation associated with transposition efficiency when confounded with biological variation.
lambda <- 0.4 #can also try estimate with find_lambda(stable_counts_filtered, nb_overdispersion, pseudo_count)
discard_pcs <- c(1)
npcs <- 25

loess_mean_var <- fit_sim_mean_var(stable_counts_filtered, nb_overdispersion, pseudo_count, lambda)

stable_counts_filtered_norm <- normalise_counts(stable_counts_filtered, nb_overdispersion, pseudo_count, lambda)
stable_counts_filtered_norm_cor <- correct_normalised_counts(stable_counts_filtered_norm, discard_pcs)

norm_df <- data.frame(mean = rowMeans(stable_counts_filtered_norm_cor), var = rowVars(stable_counts_filtered_norm_cor))
f <- ggplot(norm_df, aes(x = mean, y = var)) + geom_point() + stat_function(fun = function(x){predict(loess_mean_var, x)})


atac_stable_total_filtered <- colSums(stable_counts_filtered)

stable_counts_filtered_norm_pca <- get_pca(stable_counts_filtered_norm, npcs)

stable_total_cor <- cor(stable_counts_filtered_norm_pca$x, sqrt(atac_stable_total_filtered))[,1]

pca_df <- data.frame(pc = 1:npcs, stable_total_cor)

g <- ggplot(pca_df, aes(x = pc, y = stable_total_cor, fill = abs(stable_total_cor))) + geom_col() + 
  ylim(c(-1,1)) + scale_x_continuous(breaks = 1:npcs) + scale_fill_gradient(low = "blue", high = "yellow", limits = c(0,1))

knn_k <- 20
set.seed(123)
leiden_clusters <- iterative_cluster_sim(stable_counts_filtered, nb_overdispersion, pseudo_count, lambda, npcs, discard_pcs, knn_k, start = 0,iter.limit = 2)
umap_embeddings <- get_umap(stable_counts_filtered_norm_pca$x, npcs, discard_pcs)
cluster_outliers <- get_cluster_outliers(stable_counts_filtered_norm_pca$x, leiden_clusters, npcs, discard_pcs, knn_k, 0.2)

cluster_df <- data.frame(umap_embeddings, cluster = leiden_clusters, stable_total = atac_stable_total_filtered,  is_outlier = cluster_outliers)
p1 <- ggplot(cluster_df, aes(x = UMAP_1, y= UMAP_2, color = cluster)) + geom_point() + scale_color_manual(values = as.character(pals::polychrome(20)))
p2 <- ggplot(cluster_df, aes(x = UMAP_1, y= UMAP_2, color = stable_total)) + geom_point() + scale_color_gradient(low = "blue", high = "yellow", trans = "log10")
p3 <- ggplot(cluster_df, aes(x = UMAP_1, y= UMAP_2, color = is_outlier)) + geom_point()




qs_save(stable_counts_filtered, "binned_10Mbp_counts.qs2")
saveRDS(leiden_clusters, "leiden_clusters_10Mbp_iter2.Rds")
saveRDS(umap_embeddings, "umap_embeddings_10Mbp.Rds")
saveRDS(cluster_outliers, "cluster_outliers_10Mbp_iter2.Rds")

get_copy_number <- function(count_matrix, clusters, ref.cluster, joint, is_ref_female, ploidy.correction.method = "prop", external.ref = NULL){
  #require("MatrixGenerics")
  #removed bias correction - need to investigate whether it improves calling. Scaling to mean in presence of zeroes may be a problem.
  cor.list <- list()
  if (is.null(ref.cluster) & !is.null(external.ref)){
    ref.cluster <- "external"
    if (ploidy.correction.method == "prop"){
      cor.list[["external"]] <- rowMeans(external.ref) / sum(external.ref)
    }
  }
  for (i in levels(clusters)){
    cor.list[[i]] <- rowMeans(count_matrix[,as.integer(clusters) == i])
  }
  cn.mat <- count_matrix
  if (ploidy.correction.method == "prop"){
    for (i in seq_along(cor.list)){
      cor.list[[i]] <- cor.list[[i]] / sum(cor.list[[i]])
    }
    if (!joint){
      cn.mat <- t(t(count_matrix) / colSums(count_matrix))
    }
  }
  if (joint){
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- rep(2 * as.numeric(cor.list[[i]]) / as.numeric(cor.list[[ref.cluster]]), sum(as.integer(clusters) == i))
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- ifelse(is_ref_female, 1, 0.5) * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  } else {
    for (i in levels(clusters)){
      cn.mat[,as.integer(clusters) == i] <- 2 * cn.mat[,as.integer(clusters) == i] / as.numeric(cor.list[[ref.cluster]])
      cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i] <- ifelse(is_ref_female, 1, 0.5) * cn.mat[grep("^chrX|^chrY", rownames(cn.mat)),as.integer(clusters) == i]
    }
  }
  cn.mat[is.infinite(cn.mat)] <- NA
  cn.mat
}
test <- lapply(split(names(leiden_clusters), leiden_clusters), function(x) as.data.frame(table(srt$tumour_normal[x])))
test <- data.table::rbindlist(test, idcol = "cluster")
ref_cluster <- test$cluster[test$Var1 == "0 1"][which.max(test$Freq[test$Var1 == "0 1"])]
cat("Reference cluster is", ref_cluster, "\n")
single_cn_estimates <- get_copy_number(stable_counts_filtered, leiden_clusters, ref.cluster = ref_cluster, joint = F, external.ref = NULL, is_ref_female = F)
qs_save(single_cn_estimates, "cn_mat_10Mbp_iter2.qs2")
color_fun <- circlize::colorRamp2(breaks = c(0,2,4), colors = c("blue", "white", "red"))
library(ComplexHeatmap)
chrs <- gsub("chr", "", sapply(strsplit(rownames(single_cn_estimates),split = "_"), '[', 1))
chrs <- factor(chrs, levels = c(1:22, "X"))
ht <- Heatmap(as.matrix(t(single_cn_estimates)), col = color_fun, name = "CN est", 
        cluster_rows = F, cluster_columns = F, 
        row_split = leiden_clusters, column_split = chrs, 
        show_column_names = F, show_row_names = F)

cluster_means <- t(apply(single_cn_estimates, 1, function(gene_expr) {
  tapply(gene_expr, leiden_clusters, mean)
}))

diff_mat <- cluster_means -
  cluster_means[, ref_cluster]
ht2 <- Heatmap(as.matrix(t(diff_mat)), name = "CN diff",
               cluster_columns = F, 
               column_split = chrs, 
               show_column_names = F)

srt$ATAClone_cluster <- NA
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
s1 <- ImageDimPlot(srt, group.by = "ATAClone_cluster", cols = "polychrome") + ggtitle(paste("ref=", ref_cluster))
ggsave("ATAClone_cluster_10Mbp_iter2.png", plot = s1)



pdf("ATAClone_report_10Mbp_iter2.pdf", width = 6, height = 5)
print(f)
print(g)
print(p1)
print(p2)
print(p3)
draw(ht)
draw(ht2)
dev.off()

# merged_mat <- cbind(stable_counts_filtered, stable_counts_ref)
# qs_save(merged_mat, "merged_mat.qs2")
# anno <- data.frame(srt$tumour_normal)[colnames(stable_counts_filtered), , drop = F]
# anno <- rbind(anno, data.frame( srt.tumour_normal = rep("Breast", ncol(stable_counts_ref)), row.names = colnames(stable_counts_ref)))
# saveRDS(anno, "infercnv.anno.Rds")


 
## genome reference for the bin =========
# --- Define chromosome sizes (hg38) ---
# chrom_sizes <- data.frame(
#   chr = paste0("chr", c(1:22, "X", "Y")),
#   size = c(
#     248956422, 242193529, 198295559, 190214555, 181538259,
#     170805979, 159345973, 145138636, 138394717, 133797422,
#     135086622, 133275309, 114364328, 107043718, 101991189,
#     90338345,  83257441,  80373285,  58617616,  64444167,
#     46709983,  50818468,  156040895, 57227415
#   )
# )
# 
# # --- Generate bins ---
# bin_width <- 10e6  # 10 Mbp
# 
# bin_ref <- chrom_sizes %>%
#   rowwise() %>%
#   reframe(
#     chr   = chr,
#     start = seq(0, size - 1, by = bin_width),
#     end   = pmin(start + bin_width, size)
#   ) %>%
#   mutate(bin_name = paste0(chr, "_", start)) %>%
#   as.data.frame()
# 
# rownames(bin_ref) <- bin_ref$bin_name
# bin_ref <- bin_ref %>% select(chr, start, end)
# all(gene_bin_map$bin %in% rownames(bin_ref))
# saveRDS(bin_ref, "~/VisHD/infercnv_10Mbp_genomeref.Rds")
