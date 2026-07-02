library(dplyr)
library(dendextend, lib.loc = "~/R_Library/4.5")
library(qs2)
library(ggplot2)
library(Seurat)
library(patchwork)
library(ape)
source("~/VisHD/functions.R")
library(phylogram, lib.loc = "~/R_Library/4.5")
library(pals)
library(UCell, lib.loc = "~/R_Library/4.5")

subclone_cols <- c(setNames(brewer.set1(5), as.character(1:5)), "Normal" = "grey90", "Removed" = "grey25")
SVEC_marker <- c("SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
                 "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8")
source("~/VisHD/normal_markers.R")
tumour_markers <- c("KLK2", "KLK3", "KLK4", "TMPRSS2", "FOLH1", "NKX3-1", "HOXB13", "TRPM8")
normal_modules <- unlist(all_marker)
ucell_names <- c("tumour_score", "normal_score")
ucell_cols  <- paste(ucell_names, "UCell", sep = "_")
make_ucell_fp <- function(reduction) {
  lapply(seq_along(ucell_cols), function(j) {
    mid <- FetchData(srt, vars = ucell_cols[j])[, 1]
    mid <-max(mid)/2
    FeaturePlot(srt, features = ucell_cols[j], reduction = reduction, order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = mid, name = ucell_names[j]) +
      ggtitle(ucell_names[j]) +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))
  })
}


out_dir = "png"
test_counts <- function(k, var = "nCount_Spatial"){
    srt$test <- NA
    srt$test[names(cutree(hctree, k = k, order_clusters_as_data = FALSE))] <- cutree(hctree, k = k, order_clusters_as_data = FALSE)
    t <- VlnPlot(srt, features = var, group.by = "test", cols = setNames(polychrome(10), as.character(1:10)))+ylim(0, 500)
     t2<- VlnPlot(srt, features = var, group.by = "test", cols = setNames(polychrome(10), as.character(1:10)))
    ggsave("test.png", t/ t2, width = 5, height = 5)
}
plot_k <- function(k){
    png("test.png")
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = k)
dev.off()
}

highlight <- function(cells){
  t <- DimPlot(srt, cells.highlight = cells, cols.highlight = "red", reduction = "pearsonumap")
ggsave("test.png", t, width = 4, height = 4)
}


# LUT07 ============
setwd("~/VisHD/LUT-245-07")
srt <- qs_read("tumour_anno_srt2.qs2")
srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(2,3,6), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(2, 6), "1", ifelse(srt$pearson_clusters == 3, "2", "Normal"))
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-07", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT09 ============
setwd("~/VisHD/LUT-245-09")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

tree <- read.tree("second_round_check/infercnv/infercnv.preliminary.observations_dendrogram.txt")
hctree<-as.hclust(tree[[1]])
png("test.png")
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 10)
dev.off()
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 10, order_clusters_as_data = FALSE) %in% c(9, 10, 5))
cells <- labs[idx]
t <- DimPlot(srt, cells.highlight = cells, cols.highlight = "red", reduction = "pearsonumap")
ggsave("test.png", t, width = 4, height = 4)

srt$test <- NA
srt$test[names(cutree(hctree, k = 10, order_clusters_as_data = FALSE))] <- cutree(hctree, k = 10, order_clusters_as_data = FALSE)
t <- DimPlot(srt, group.by = "test", cols = setNames(pals::alphabet(10), as.character(1:10)), reduction = "pearsonumap")
t <- ImageDimPlot(srt, group.by = "test", cols = setNames(pals::polychrome(10), as.character(1:10)))
t <- VlnPlot(srt, features = c("tumour_score_UCell", "normal_score_UCell"), group.by = "test", cols = setNames(polychrome(10), as.character(1:10)))
ggsave("test.png", t, width = 5, height = 5)


srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(1,4,5), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(1,4), "1", ifelse(srt$pearson_clusters == 5, "2", "Normal"))
srt$tumour_anno[as.character(cells)] <- "Removed" 
srt$subclone[as.character(cells)] <- "Removed" 
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-09", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT10 ============
setwd("~/VisHD/LUT-245-10")
srt <- qs_read("tumour_anno_srt2.qs2")
tree <- read.tree("second_round_check/infercnv/infercnv.preliminary.observations_dendrogram.txt")
hctree<-as.hclust(tree[[1]])
png("test.png")
plot(hctree, labels = FALSE)
rect.hclust(hctree, k = 3)
dev.off()
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
cells1 <- cells
hctree<-as.hclust(tree[[1]])
srt$test <- NA
srt$test[names(cutree(hctree, k = 3, order_clusters_as_data = FALSE))] <- cutree(hctree, k = 3, order_clusters_as_data = FALSE)
t <- VlnPlot(srt, features = "nCount_Spatial", group.by = "test", cols = setNames(polychrome(10), as.character(1:10)))
ggsave("test.png", t, width = 5, height = 5)

hctree<-as.hclust(tree[[2]])
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
cells2 <- cells

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(4,6), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c( 6), "1", ifelse(srt$pearson_clusters == 4, "2", "Normal"))
srt$tumour_anno[as.character(c(cells1, cells2))] <- "Removed" 
srt$subclone[as.character(c(cells1, cells2))] <- "Removed" 
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-10", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT11 ============
setwd("~/VisHD/LUT-245-11")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

srt <- AddModuleScore_UCell(srt,
    features = list(SVEC = SVEC_marker), assay = "Spatial")
t <- FeaturePlot(srt, features = "SVEC_UCell", reduction = "pearsonumap", order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = max(FetchData(srt, vars = "SVEC_UCell")[, 1])/2, name = "SVEC_UCell") +
      ggtitle("SVEC_UCell") +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))    
ggsave(file.path(out_dir, "SVEC.png"),
       width =  3, height = 3, dpi = 400)

t <- VlnPlot(srt, features = "nCount_Spatial", group.by = "pearson_clusters", cols = setNames(polychrome(10), as.character(1:10)))
ggsave("test.png", t, width = 5, height = 5)


srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(0,1,4), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(0, 1,4), "1", "Normal")

  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-11", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT15 ============
setwd("~/VisHD/LUT-245-15")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

srt <- AddModuleScore_UCell(srt,
    features = list(SVEC = SVEC_marker), assay = "Spatial")
t <- FeaturePlot(srt, features = "SVEC_UCell", reduction = "pearsonumap", order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = max(FetchData(srt, vars = "SVEC_UCell")[, 1])/2, name = "SVEC_UCell") +
      ggtitle("SVEC_UCell") +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))    
ggsave(file.path(out_dir, "SVEC.png"),
       width =  3, height = 3, dpi = 400)

tree <- read.tree("second_round_check/infercnv/infercnv.preliminary.observations_dendrogram.txt")
hctree<-as.hclust(tree[[1]])
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 5, order_clusters_as_data = FALSE) %in% c(3))
cells <- labs[idx]
cells1 <- cells
hctree<-as.hclust(tree[[1]])


hctree<-as.hclust(tree[[2]])
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 6, order_clusters_as_data = FALSE) %in% c(3))
cells <- labs[idx]
cells2 <- cells

hctree<-as.hclust(tree[[3]])
test_counts(8)
plot_k(8)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 8, order_clusters_as_data = FALSE) %in% c(3))
cells <- labs[idx]
cells3 <- cells

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(0, 4, 6, 8), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(0, 6, 8), "1", ifelse(srt$pearson_clusters == 4, "2", "Normal"))
srt$tumour_anno[as.character(c(cells1, cells2, cells3))] <- "Removed" 
srt$subclone[as.character(c(cells1, cells2, cells3))] <- "Removed" 
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-15", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")


# LUT16 ============
setwd("~/VisHD/LUT-245-16")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

srt <- AddModuleScore_UCell(srt,
    features = list(SVEC = SVEC_marker), assay = "Spatial")
t <- FeaturePlot(srt, features = "SVEC_UCell", reduction = "pearsonumap", order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = max(FetchData(srt, vars = "SVEC_UCell")[, 1])/2, name = "SVEC_UCell") +
      ggtitle("SVEC_UCell") +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))    
ggsave(file.path(out_dir, "SVEC.png"),
       width =  3, height = 3, dpi = 400)

tree <- read.tree("second_round_check/infercnv/infercnv.preliminary.observations_dendrogram.txt")
hctree<-as.hclust(tree[[3]])
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(2))
cells <- labs[idx]
removecells <- cells


hctree<-as.hclust(tree[[4]])
test_counts(2)
plot_k(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(removecells, cells)

hctree<-as.hclust(tree[[6]])
test_counts(3)
plot_k(3)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(removecells, cells)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(2))
cells <- labs[idx]
tumourcells <- cells

hctree<-as.hclust(tree[[7]])
test_counts(6)
plot_k(6)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 6, order_clusters_as_data = FALSE) %in% c(4))
cells <- labs[idx]
removecells <- c(removecells, cells)
idx   <- which(cutree(hctree, k = 6, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
tumourcells2 <-  cells

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(1,2, 4,5, 6, 9), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(1,2,5, 6, 9), "1", ifelse(srt$pearson_clusters == 4, "2", "Normal"))
srt$tumour_anno[as.character(removecells)] <- "Removed" 
srt$subclone[as.character(removecells)] <- "Removed" 
srt$tumour_anno[as.character(c(tumourcells, tumourcells2))] <- "Tumour"
srt$subclone[as.character(tumourcells)] <- "1"
srt$subclone[as.character(tumourcells2)] <- "2"
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-16", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT17 ============
setwd("~/VisHD/LUT-245-17")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

srt <- AddModuleScore_UCell(srt,
    features = list(SVEC = SVEC_marker), assay = "Spatial")
t <- FeaturePlot(srt, features = "SVEC_UCell", reduction = "pearsonumap", order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = max(FetchData(srt, vars = "SVEC_UCell")[, 1])/2, name = "SVEC_UCell") +
      ggtitle("SVEC_UCell") +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))    
ggsave(file.path(out_dir, "SVEC.png"),
       width =  3, height = 3, dpi = 400)

tree <- read.tree("second_round_check/infercnv/infercnv.preliminary.observations_dendrogram.txt")
hctree<-as.hclust(tree[[2]])
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- cells


hctree<-as.hclust(tree[[3]])
test_counts(3)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(removecells, cells)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(2))
cells <- labs[idx]
tumourcells <- cells

hctree<-as.hclust(tree[[4]])
test_counts(3)
plot_k(3)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(removecells, cells)

hctree<-as.hclust(tree[[5]])
test_counts(3)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 3, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(removecells, cells)

hctree<-as.hclust(tree[[6]])
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(2))
cells <- labs[idx]
removecells <- c(removecells, cells)

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(1,2, 4,5, 9), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(1,2, 4,5, 9), "1",  "Normal")
srt$tumour_anno[as.character(removecells)] <- "Removed" 
srt$subclone[as.character(removecells)] <- "Removed" 
srt$tumour_anno[as.character(tumourcells)] <- "Tumour"
srt$subclone[as.character(tumourcells)] <- "1"
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-17", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")

# LUT20 ============
setwd("~/VisHD/LUT-245-20")
srt <- qs_read("tumour_anno_srt2.qs2")
srt <- AddModuleScore_UCell(srt,
    features = list(tumour_score = tumour_markers, normal_score = normal_modules), assay = "Spatial")
ggsave(file.path(out_dir, "FeaturePlot_tumour_normal.png"),
       wrap_plots(make_ucell_fp("pearsonumap"), nrow = 1) +
         plot_annotation(title = "tumour/normal scores (UCell)"),
       width = length(ucell_cols) * 4, height = 3, dpi = 400)

srt <- AddModuleScore_UCell(srt,
    features = list(SVEC = SVEC_marker), assay = "Spatial")
t <- FeaturePlot(srt, features = "SVEC_UCell", reduction = "pearsonumap", order = TRUE, pt.size = 0.01) +
      scale_color_gradient2(low = "navy", mid = "white", high = "darkred",
                            midpoint = max(FetchData(srt, vars = "SVEC_UCell")[, 1])/2, name = "SVEC_UCell") +
      ggtitle("SVEC_UCell") +
      theme(plot.title = element_text(size = 9), legend.key.width = unit(0.3, "cm"))    
ggsave(file.path(out_dir, "SVEC.png"),
       width =  3, height = 3, dpi = 400)

infercnvobj <- readRDS("second_round_check/infercnv/run.final.infercnv_obj")
tree <- infercnvobj@tumor_subclusters$hc
names(tree)

hctree <- tree[["10"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- cells

hctree <- tree[["9"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(2))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["8"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["6"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["5"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["4"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["2"]]
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["0"]]
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["7"]]
test_counts(2)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)
plot_k(4)
idx   <- which(cutree(hctree, k = 4, order_clusters_as_data = FALSE) %in% c(3))
cells <- labs[idx]
tumourcells <- cells

hctree <- tree[["3"]]
test_counts(4, "tumour_score_UCell")
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

hctree <- tree[["1"]]
test_counts(4)
labs <- labels(hctree)
idx   <- which(cutree(hctree, k = 2, order_clusters_as_data = FALSE) %in% c(1))
cells <- labs[idx]
removecells <- c(cells, removecells)

srt$tumour_anno <- ifelse(srt$pearson_clusters %in% c(1,3), "Tumour", "Normal")
srt$subclone <- ifelse(srt$pearson_clusters %in% c(1), "1",  ifelse(srt$pearson_clusters %in% c(3), "2", "Normal"))
srt$tumour_anno[as.character(removecells)] <- "Removed" 
srt$subclone[as.character(removecells)] <- "Removed" 
srt$tumour_anno[as.character(tumourcells)] <- "Tumour"
srt$subclone[as.character(tumourcells)] <- "1"
  n_tab <- table(srt$tumour_anno)
  subtitle <- paste0("Tumour: ", n_tab["Tumour"], "   Normal: ", n_tab["Normal"], "   Removed: ", n_tab["Removed"])
p = ImageDimPlot(srt, group.by = "tumour_anno", cols = c("Tumour" = "red", "Normal" = "grey90", "Removed" = "grey25")) + 
ImageDimPlot(srt, group.by = "subclone", cols = subclone_cols)
p + plot_annotation(title = "LUT-245-20", subtitle = subtitle)
ggsave("png/subclone_anno.png", width = 8, height = 3)
qs_save(srt, "tumour_subclone_srt.qs2")