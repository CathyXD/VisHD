library(Seurat)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

paths <- system("realpath ~/VisHD/LUT-*/", intern = T)

# Select the specific sample directory based on the command line argument
path <- paths[arg]
setwd(path)

infercnvobj <- readRDS("bined_ouput/infercnv_uncluster/run.final.infercnv_obj")

# 1. Extract the matrix (genes as rows, cells as columns)
# Note: infercnv_obj@expr.data contains the smoothed/normalized values
cnv_matrix <- infercnvobj@expr.data

# 2. Create a temporary Seurat object for clustering

cnv_seurat <- CreateSeuratObject(counts = cnv_matrix, data = cnv_matrix)

# 3. Standardize and Run PCA
# We skip normalization as infercnv data is already processed
cnv_seurat <- ScaleData(cnv_seurat, features = rownames(cnv_seurat), do.scale = F, do.center = F)
cnv_seurat <- RunPCA(cnv_seurat, features = rownames(cnv_seurat), approx = FALSE)

# 4. Graph-based Clustering (Louvain/Leiden)
cnv_seurat <- FindNeighbors(cnv_seurat, dims = 1:20)
# Use algorithm = 1 for Louvain, algorithm = 4 for Leiden (requires leiden pocket)
cnv_seurat <- FindClusters(cnv_seurat, resolution = 0.1, algorithm = 1)
subclone_assignments <- Idents(cnv_seurat)
saveRDS(subclone_assignments, "Infercnv_subclone.Rds")


# 1. Ensure metadata is aligned
# Assuming cnv_seurat has the subclone labels in Idents()
# and srt_cell has your BANKSY clusters
srt_cell$cnv_subclone <- Idents(cnv_seurat)

# 2. Create the nested grouping variable
srt_cell$nested_group <- paste0("Sub:", srt_cell$cnv_subclone, 
                                "_Bnk:", srt_cell$BANKSY_0.2_snn_res.1)

# 3. Aggregate Expression (Average)
# Use the @expr.data (CNV ratios) slot
# We use 'SCT' or 'RNA' assay depending on where you stored the cnv_matrix
avg_cnv <- AverageExpression(srt_cell, 
                             group.by = "nested_group", 
                             features = rownames(cnv_matrix),
                             slot = "data")$RNA # Change to appropriate assay name

# 4. Prepare Metadata for Heatmap Annotation
# Extract the original labels back from the column names of avg_cnv
meta_split <- data.frame(row.names = colnames(avg_cnv)) %>%
  mutate(
    group_name = colnames(avg_cnv),
    subclone = sub("_Bnk:.*", "", sub("Sub:", "", group_name)),
    banksy_cluster = sub(".*_Bnk:", "", group_name)
  )

# Define Colors
col_subclone <- setNames(RColorBrewer::brewer.pal(8, "Set1"), unique(meta_split$subclone))
col_banksy <- setNames(RColorBrewer::brewer.pal(12, "Paired"), unique(meta_split$banksy_cluster))

# Create Top Annotation
top_anno <- HeatmapAnnotation(
  Subclone = meta_split$subclone,
  Banksy = meta_split$banksy_cluster,
  col = list(Subclone = col_subclone, Banksy = col_banksy),
  show_legend = TRUE
)

# 5. Generate Heatmap
# We center the color scale at 1 (neutral CNV)
ht <- Heatmap(as.matrix(avg_cnv),
        name = "Avg CNV Ratio",
        column_title = "Subclones nested with BANKSY Clusters",
        top_annotation = top_anno,
        
        # Clustering & Splitting
        cluster_columns = TRUE, 
        column_split = meta_split$subclone, # Groups by Subclone visually
        cluster_rows = FALSE,               # Usually keep rows in genomic order
        
        # Aesthetics
        show_column_names = TRUE,
        show_row_names = FALSE,
        col = colorRamp2(c(0.8, 1, 1.2), c("#313695", "white", "#a50026")),
        
        # Rasterize if the number of genes is very high
        use_raster = TRUE)
pdf("Infercnv_heatmap.pdf", width =10, height = 7)
draw(ht)
dev.off()