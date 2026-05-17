library(dplyr)
library(dendextend)
library(qs2)
library(qs)
library(ggplot2)
library(Seurat)
library(patchwork)
source("~/VisHD/functions.R")
library(phylogram)
library(pals)
tumour_cluster <- list(c(4, 5, 7, 10), #07
                       c(4, 5, 7, 8, 10),#09
                       c(3, 5), #10
                       c(1), #11
                       c(4, 6, 8, 9, 10), #15
                       c(3, 4),  #16, 
                       c(2, 6), #17
                       c(5, 9, 13, 10, 14)#20
                       
)

normal_cluster <- list(c(1,2,3, 6, 8), 
                       c(1,2,3, 6), 
                       c(1,2,4, 6, 8),
                       c(2, 3), 
                       c(1,2, 3, 5, 7), 
                       c(1, 2), 
                       c(1, 3, 4, 5), 
                       c(1:4,6:8, 11)
)


# LUT07 ============
setwd("~/VisHD/LUT-245-07")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[1]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[1]], "Normal", "Removed"))
srt$ATAClone_cluster <- factor(srt$ATAClone_cluster)
srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% c(4,5)] <- 1
srt$subclone[srt$ATAClone_cluster %in% c(7, 10)] <- 2
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-07")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

bintocell <- function(srt_bin, sample_name){
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
  ggsave("png/subclone_anno.png", width = 12, height = 4)
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
  ggsave("png/tumour_anno_proportion.png", width = 6, height = 4)
  return(srt_cell_filtered)
}

bintocell(srt, sample_name = "LUT-245-07")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")

# LUT09 ============
setwd("~/VisHD/LUT-245-09")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% c(4,7,8), "Tumour", ifelse(srt$ATAClone_cluster %in% c(1:3, 5,6,9:11), "Normal", "Removed"))

srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% c(4)] <- 1
srt$subclone[srt$ATAClone_cluster %in% c(7)] <- 2
srt$subclone[srt$ATAClone_cluster %in% c(8)] <- 3

p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-09")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-09")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")
# LUT10 ============
setwd("~/VisHD/LUT-245-10")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
srt$ATAClone_cluster <- factor(srt$ATAClone_cluster)
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[1]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[1]], "Normal", "Removed"))
tree <- ape::read.tree("bined_ouput/infercnv_ATAClone/infercnv.observations_dendrogram.txt")
dend_list <- lapply(tree, as.dendrogram)
plot(color_branches(dend_list[[1]], k = 4))
tree_cluster <- cutree(dend_list[[1]], 4)
mergingcells <- names(tree_cluster)[tree_cluster %in% c(2, 3)]

srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% c(3)] <- 1
srt$subclone[srt$ATAClone_cluster %in% c(5)] <- 2
srt$subclone[mergingcells] <- 2

p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-10")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

bintocell(srt, sample_name = "LUT-245-10")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")

# LUT11 ============
setwd("~/VisHD/LUT-245-11")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
#need to load infercnv obj 
infercnvobj <- readRDS("bined_ouput/infercnv_ATAClone/run.final.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_ATAClone/",  write_expr_matrix = F, save_objects =T)
hctree <- infercnvobj@tumor_subclusters$hc$`2` # reference 
plot(hctree, labels = FALSE)
tree <- cutree(hctree, k = 3, order_clusters_as_data = F)
tumourcells <- names(tree[tree == 2])

hctree <- infercnvobj@tumor_subclusters$hc$`1`
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 3, border = c("red", "blue", "green"))
tree <- cutree(hctree, k = 3, order_clusters_as_data = F)
tumourcells <- c(tumourcells, names(tree)[tree == 3])
junkcells <- names(tree)[tree %in% c(1, 2)]

srt$tumour_anno <- ifelse(colnames(srt) %in% tumourcells, "Tumour", "Normal")
srt$tumour_anno[junkcells] <- "Removed"
srt$subclone <- NA
srt$subclone[tumourcells] <- 1

p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-11")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-11")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")
#LUT 15 ==========
# need to plot the reference by group when loaded infercnv
setwd("~/VisHD/LUT-245-15")
infercnvobj <- readRDS("bined_ouput/infercnv_ATAClone/run.final.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_ATAClone/",  write_expr_matrix = F, save_objects =T)
#clear 
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
srt$ATAClone_cluster <- factor(srt$ATAClone_cluster)
srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[5]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[5]], "Normal", "Removed"))
srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% c(6, 8)] <- 1
srt$subclone[srt$ATAClone_cluster %in% c(4, 9, 10)] <- 2
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-15")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-15")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")

#LUT 16 ==========
setwd("~/VisHD/LUT-245-16")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
#need to load infercnv obj 
infercnvobj <- readRDS("bined_ouput/infercnv_ATAClone/run.final.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_ATAClone/",  write_expr_matrix = F, save_objects =T)
hctree <- infercnvobj@tumor_subclusters$hc$`1` # reference 
plot(hctree, labels = FALSE)
tree <- cutree(hctree, k = 3, order_clusters_as_data = F)
tumourcells <- names(tree[tree == 2])

hctree <- infercnvobj@tumor_subclusters$hc$`3` # reference 
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 5, border = c("red" ))
tree <- cutree(hctree, k = 5, order_clusters_as_data = F)
subclcells <- names(tree[tree == 1])

srt$tumour_anno <- ifelse(colnames(srt) %in% tumourcells|srt$ATAClone_cluster %in% tumour_cluster[[6]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[6]], "Normal", "Removed"))
srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% c(4)] <- 3
srt$subclone[srt$ATAClone_cluster %in% c(3)] <- 1
srt$subclone[c(tumourcells, subclcells)] <- 2
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-16")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)

srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-16")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")

#LUT 17 ==========
setwd("~/VisHD/LUT-245-17")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[7]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[7]], "Normal", "Removed"))
srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% tumour_cluster[[7]]] <- 1

p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-17")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)
srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-17")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")

#LUT 20 ==========
setwd("~/VisHD/LUT-245-20")
srt <- qs_read("bined_ouput/srt_infercnv.qs2")
leiden_clusters <- readRDS("bined_ouput/leiden_clusters_10Mbp_iter2.Rds")
srt$ATAClone_cluster <- 0
srt$ATAClone_cluster[names(leiden_clusters)] <- leiden_clusters
ImageDimPlot(srt, group.by= "ATAClone_cluster", cols = "polychrome")
infercnvobj <- readRDS("bined_ouput/infercnv_ATAClone/run.final.infercnv_obj")
infercnv::plot_per_group(infercnvobj, on_observations = F, out_dir ="bined_ouput/infercnv_ATAClone/",  write_expr_matrix = F, save_objects =T)
hctree <- infercnvobj@tumor_subclusters$hc$`5` # reference 
plot(hctree, labels = FALSE)
tree <- cutree(hctree, k = 2, order_clusters_as_data = F)
normalcells <- names(tree[tree == 2])

srt$tumour_anno <- ifelse(srt$ATAClone_cluster %in% tumour_cluster[[8]], "Tumour", ifelse(srt$ATAClone_cluster %in% normal_cluster[[8]]|colnames(srt) %in% normalcells, "Normal", "Removed"))
srt$subclone <- NA
srt$subclone[srt$ATAClone_cluster %in% tumour_cluster[[8]]] <- 1
srt$subclone[srt$ATAClone_cluster %in% c(13)] <- 2
srt$subclone[srt$ATAClone_cluster %in% c(14, 10)] <- 3
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
  ImageDimPlot(srt, group.by = "subclone")
p + plot_annotation(title = "LUT-245-20")
ggsave("bined_ouput/subclone_anno.png", width = 8, height = 3)
srt_cell_filtered <- bintocell(srt, sample_name = "LUT-245-20")
srt_filtered <- subset(srt, cells = colnames(srt)[srt$tumour_anno != "Removed"])
srt_filtered$anno <- paste(srt_filtered$tumour_anno, ifelse(is.na(srt_filtered$subclone), "", srt_filtered$subclone))
srt_filtered <- filter_artefacts_knn(srt_filtered, min_neighbours  = 5)
qs_save(srt_filtered, "bined_ouput/subclone_srt.qs2")