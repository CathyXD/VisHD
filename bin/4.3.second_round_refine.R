library(patchwork)
library(infercnv)
library(qs2)
library(dendextend, lib.loc = "~/R_Library/4.5")
library(Seurat)
library(ggplot2)
library(qs, lib.loc = "~/R_Library/4.5")
source("~/VisHD/functions.R")

paths <- system("realpath ~/VisHD/LUT-245-*/bined_ouput/infercnv_subclone", intern = T)
for (path in paths[5:6]){
  setwd(path)
  infercnvobj <- readRDS("preliminary.infercnv_obj")
  saveRDS(infercnvobj@tumor_subclusters$hc, "all_dendrogram.Rds")
  infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir =".",  write_expr_matrix = F, save_objects =F, output_format = "pdf", k_obs_groups = 2)
}

bintocell <- function(srt_bin, sample_name){
  require(dplyr)
  require(pals)
  srt_cell <- qread("raw_srt.qs")
  cell_types <- read.csv(paste0("~/VisHD/raw/",sample_name,"/outs/segmented_outputs/cell_types/Azimuth/cell_types.csv"))
  simple_ids <- as.integer(gsub("cellid_|-1", "", cell_types$barcode))
  cell_types$cell_id <- simple_ids
  
  CB <- read.csv("category.csv")
  CB$cellid <- as.integer(gsub("cellid_|-1", "", CB$Barcode))
  srt_cell$category <- NA
  srt_cell$category <- CB$category[match(srt_cell$cell, CB$cellid)]
  srt_cell@meta.data <- cbind(srt_cell@meta.data, cell_types[match(srt_cell$cell, CB$cellid), ])
  srt_cell <- transfer_visiumhd_to_cells(
    srt_cell,
    srt,
    annotation_cols = c("tumour_score_UCell", "normal_score_UCell", "tumour_normal", "ATAClone_cluster", "tumour_anno", "subclone")
  )
  srt_cell <- UpdateSeuratObject(srt_cell)
  srt_cell_filtered <- subset(srt_cell, cells  = colnames(srt_cell)[srt_cell$nFeature_Spatial>20& srt_cell$cell_area > quantile(srt_cell$cell_area, 0.05)&srt_cell$cell_area < quantile(srt_cell$cell_area, 0.99) & srt_cell$tumour_anno != "Removed"])
  srt_cell_filtered <- filter_artefacts_knn(srt_cell_filtered, min_neighbours  = 5)
  s1 = ImageDimPlot(srt_cell_filtered, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) 
  s2 = ImageDimPlot(srt_cell_filtered,  group.by = "ATAClone_cluster", cols = "polychrome")
  s3 = ImageDimPlot(srt_cell_filtered,  group.by = "subclone")
  labs <-table(srt_cell_filtered$tumour_anno)
  print(s1 + s2 + s3 + plot_annotation(title = sample_name, 
                                       subtitle = paste("n Normal = ", labs["Normal"], "|n Tumour =", labs["Tumour"])))
  ggsave("png/subclone_anno2.png", width = 12, height = 4)
  qs_save(srt_cell_filtered, "tumour_anno_srt.qs2")
  prop_tumour <- srt_cell_filtered@meta.data %>%
    group_by(category, tumour_anno) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(category) %>%
    mutate(proportion = n / sum(n))
  
  g1 <- ggplot(prop_tumour, aes(x = category, y = proportion, fill = tumour_anno)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Tumour Annotation",
      x = "Category",
      y = "Proportion",
      fill = "Tumour Anno"
    ) +
    theme_classic() +
    scale_fill_manual(values = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  # ----- Plot 2: Proportion of subclone by category -----
  
  prop_subclone <-  srt_cell_filtered@meta.data %>%
    filter(!is.na(subclone)) %>% 
    mutate(subclone = as.factor(subclone)) %>%
    group_by(category, subclone) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(category) %>%
    mutate(proportion = n / sum(n))
  
  g2 <- ggplot(prop_subclone, aes(x = category, y = proportion, fill = subclone)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Subclone",
      x = "Category",
      y = "Proportion",
      fill = "Subclone"
    ) +
    theme_classic() +
    scale_fill_manual(values = brewer.accent(5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g1 + g2)
  ggsave("png/tumour_anno_proportion2.png", width = 6, height = 4)
  return(srt_cell_filtered)
}

#LUT07
setwd("~/VisHD/LUT-245-07")
infercnvobj <- readRDS("bined_ouput/infercnv_subclone/preliminary.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_subclone/",  write_expr_matrix = F, save_objects =F, output_format = "pdf", k_obs_groups = 1)
saveRDS(infercnvobj@tumor_subclusters$hc, "bined_ouput/infercnv_subclone/all_dendrogram.Rds")
hctree <- infercnvobj@tumor_subclusters$hc$Normal
tree <- cutree(hctree, k = 8, order_clusters_as_data = F)
png("test.png")
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 8, border = c("red"))
dev.off()

table(tree)
tumourcells <- names(tree[tree == 4])
srt <- qs_read("bined_ouput/subclone_srt.qs2")
srt$tumour_anno <- ifelse(colnames(srt) %in% tumourcells, "Tumour", srt$tumour_anno)
srt$subclone <- ifelse(colnames(srt) %in% tumourcells, 1, srt$subclone)
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-07")
ggsave("bined_ouput/subclone_anno2.png", plot = p,  width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, "LUT-245-07")
qs_save(srt, "bined_ouput/subclone_srt.qs2")

#LUT10
setwd("~/VisHD/LUT-245-10")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt <- qs_read("bined_ouput/subclone_srt.qs2")
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-10")

all_dendrogram <- readRDS("~/VisHD/LUT-245-10/bined_ouput/infercnv_subclone/all_dendrogram.Rds")
hctree <- all_dendrogram$Normal
plot(hctree, labels = FALSE)
rect.hclust(hctree, h = 5, border = c("red"))
tree <- cutree(hctree,h = 5, order_clusters_as_data = F)
tumourcells <- names(tree[tree == 16])

hctree <- all_dendrogram$`1`
plot(hctree, labels = FALSE)
tree <- cutree(hctree,k = 3, order_clusters_as_data = F)
subclone2 <- names(tree[tree == 3])
srt$subclone[subclone2] <- 2
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-10")
ggsave("bined_ouput/subclone_anno2.png", width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, "LUT-245-10")
qs_save(srt, "bined_ouput/subclone_srt.qs2")

#LUT11
setwd("~/VisHD/LUT-245-11")
# Reference clsutering unclear showing tumour cells so skip 


#LUT15
setwd("~/VisHD/LUT-245-15")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
all_dendrogram <- readRDS("bined_ouput/infercnv_subclone/all_dendrogram.Rds")
hctree <- all_dendrogram$Normal
tree <- cutree(hctree,k = 2, order_clusters_as_data = F)
unkown <- names(tree[tree ==1])

#LUT16
setwd("~/VisHD/LUT-245-16")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt <- qs_read("bined_ouput/subclone_srt.qs2")
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-16")

all_dendrogram <- readRDS("bined_ouput/infercnv_subclone/all_dendrogram.Rds")
ref_tree <- all_dendrogram$Normal
ref_tree <- cutree(ref_tree, 3, order_clusters_as_data = F)
ref_tree[ref_tree %in% c(1,3)] <- NA
ref_tree[ref_tree == 2] <- 1

new_subclone <- all_dendrogram$`2`
new_subclone <- cutree(new_subclone, k = 2)
new_subclone[new_subclone == 2] <- 4
new_subclone[new_subclone == 1] <- 2

tumourcells1 <- cutree(all_dendrogram$`1`, k = 1)
tumourcells3 <- setNames(rep(3, length(labels(all_dendrogram$`3`))),labels(all_dendrogram$`3`)) 

allcells <- c(ref_tree, new_subclone, tumourcells1, tumourcells3)
srt$tumour_anno[names(which(is.na(allcells)))] <- "Normal"
srt$subclone[names(allcells)] <- allcells
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-16")

srt_cell_filtered <- bintocell(srt, "LUT-245-16")
qs_save(srt, "bined_ouput/subclone_srt.qs2")

#LUT20
setwd("~/VisHD/LUT-245-20")
infercnvobj <- readRDS("bined_ouput/infercnv_subclone/preliminary.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_subclone/",  write_expr_matrix = F, save_objects =F, output_format = "pdf", k_obs_groups = 1)
saveRDS(infercnvobj@tumor_subclusters$hc, "bined_ouput/infercnv_subclone/all_dendrogram.Rds")
hctree <- infercnvobj@tumor_subclusters$hc$Normal
tree <- cutree(hctree, k = 8, order_clusters_as_data = F)
png("test.png")
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 8, border = c("red"))
dev.off()

table(tree)
tumourcells <- names(tree[tree == 4])
srt <- qs_read("bined_ouput/subclone_srt.qs2")
srt$tumour_anno <- ifelse(colnames(srt) %in% tumourcells, "Tumour", srt$tumour_anno)
srt$subclone <- ifelse(colnames(srt) %in% tumourcells, 1, srt$subclone)
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-20")
ggsave("bined_ouput/subclone_anno2.png", plot = p,  width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, "LUT-245-20")
qs_save(srt, "bined_ouput/subclone_srt.qs2")