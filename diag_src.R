library(qs2); library(Seurat)

cat("######## NORMAL SOURCE: 9.4.1 normal_srt_final_anno.qs2 ########\n")
n <- qs_read(path.expand("~/VisHD/8.4.final_clear_normal_integration/normal_srt_final_anno.qs2"))
nm <- n@meta.data
cat("cells:", nrow(nm), "\n")
cat("has x_centroid col:", "x_centroid" %in% colnames(nm), " y:", "y_centroid" %in% colnames(nm), "\n")
if ("x_centroid" %in% colnames(nm)) {
  cat("NA x_centroid:", sum(is.na(nm$x_centroid)), "/", nrow(nm), "\n")
  cat("NA by slide:\n"); print(tapply(is.na(nm$x_centroid), nm$slide, sum))
  cat("total by slide:\n"); print(table(nm$slide))
}
rm(n); gc()

cat("\n######## TUMOUR SOURCE: LUT-245-10/tumour/tumour_srt.qs2 ########\n")
t <- qs_read(path.expand("~/VisHD/LUT-245-10/tumour/tumour_srt.qs2"))
t <- UpdateSeuratObject(t)
cat("tumour cells:", ncol(t), "\n")
co <- GetTissueCoordinates(t, which = "centroids")
cat("GetTissueCoordinates: nrow =", nrow(co), " cols:", paste(colnames(co), collapse=","), "\n")
cat("rownames(co) head:", paste(head(rownames(co)), collapse=" | "), "\n")
cat("co$cell head:", paste(head(co$cell), collapse=" | "), "\n")
cat("colnames(t) head:", paste(head(colnames(t)), collapse=" | "), "\n")
cat("colnames in co$cell:", sum(colnames(t) %in% co$cell), "/", ncol(t), "\n")
cat("colnames in rownames(co):", sum(colnames(t) %in% rownames(co)), "/", ncol(t), "\n")
cat("any duplicated co$cell:", anyDuplicated(co$cell), "\n")
