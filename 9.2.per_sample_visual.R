suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SpaNorm,    lib.loc = "~/R_Library/4.5")
  library(qs2)
  library(leidenbase, lib.loc = "~/R_Library/4.5")
  library(UCell,      lib.loc = "~/R_Library/4.5")
  library(ggplot2)
  library(readxl)
  library(purrr)
  library(patchwork)
  library(stringr)
  library(pals)
})

source("~/VisHD/functions.R")
source("~/VisHD/normal_markers.R")   # all_marker

paths <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
heat.shock_mod <- c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "NR4A1", "DUSP1", "ZFP36", "HSPA1A","HSPA1B", "IER3","IER2", "BTG2", "SOCS3")
ifp_save <- function(srt, feats, outfile, ncol = 3, dpi = 300) {
  feats <- feats[feats %in% colnames(srt@meta.data)]
  if (length(feats) == 0) return(invisible(NULL))
  nrows <- ceiling(length(feats) / ncol)
  plots <- lapply(feats, function(f) {
    v   <- srt@meta.data[[f]]
    mid <- (max(v, na.rm = TRUE) + min(v, na.rm = TRUE)) *0.3
    ImageFeaturePlot(srt,  features = f, size = 0.3) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = mid)+
      scale_colour_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = mid)
  })
  p <- patchwork::wrap_plots(plots, ncol = ncol)
  ggsave(outfile, p,
         width = ncol * 8, height = nrows * 7,
         dpi = dpi, limitsize = FALSE)
}
meta_xlsx     <- "~/VisHD/public_signature/meta_programs_2025-01-29.xlsx"
sheetname     <- excel_sheets(meta_xlsx)
Imagescore_meta_programs <- function(obj, sheetname, meta_cols, out_dir,
                                tag = "") {
  require(Seurat)
  require(ggplot2)

  for (sheet in sheetname) {
    # Columns for this sheet: ucell_score named them _meta<SheetName>.<ProgramName>_UCell
    sheet_pat  <- paste0("^", make.names(sheet))
    sheet_cols <- grep(sheet_pat, meta_cols, value = TRUE)
    sheet_cols <- sheet_cols[sheet_cols %in% colnames(obj@meta.data)]
    if (length(sheet_cols) == 0) next

    ncol_g <- min(3L, length(sheet_cols))
    nrow_g <- ceiling(length(sheet_cols) / ncol_g)
    base   <- make.names(sheet)
    plots_s <- lapply(sheet_cols, function(f) {
      v   <- obj@meta.data[[f]]
      mid <- (max(v, na.rm = TRUE) + min(v, na.rm = TRUE)) * 0.3
      tt <- gsub(paste0(sheet, "."), fix = T, "", f)
      tt <- gsub(".", " ",fix = T, tt)
      tt <- gsub("_meta_UCell", "", tt)
      ImageFeaturePlot(obj, features = f) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                             midpoint = mid)+
        scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                              midpoint = mid) + ggtitle(tt)
    })
    g <- patchwork::wrap_plots(plots_s, ncol = ncol_g)
    ggsave(file.path(out_dir, paste0("meta_", base, tag, ".png")),
           plot = g, width = ncol_g * 8, height = nrow_g * 7,
           dpi = 300, limitsize = FALSE)
  }
}


for (arg in 1:8){
  path  <- paths[arg]
  i     <- basename(path)
  setwd(path)
  srt            <- qs_read("tumour_normal_anno_srt.qs2")
  spatial_dir  <- file.path(path, "final_png", "spatial")
  srt <- AddModuleScore_UCell(srt, features = list(heat.shock_mod), name = "_heat_shock_mod", slot = "data", assay = "SpaNorm")

  arch_mod_cols <- grep("_arch", colnames(srt@meta.data), value = TRUE)
  mod_score_cols<- grep("_gd", colnames(srt@meta.data), value = TRUE)
  gs23_cols <-  grep("_gs23", colnames(srt@meta.data), value = TRUE)
  meta_cols <-  grep("_meta", colnames(srt@meta.data), value = TRUE)
  tn_cols <- c("tumour_score_UCell", "normal_score_UCell")



  ifp_save(srt, feats = arch_mod_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_archetype_modules.png"))
  ifp_save(srt, feats = mod_score_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_G123_scores.png"))
  ifp_save(srt, feats = gs23_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_genesets2023.png"))
  ifp_save(srt, feats = tn_cols, outfile = file.path(spatial_dir, "ImageFeaturePlot_tumour_normal_scores.png"))
  Imagescore_meta_programs(srt, sheetname = sheetname, meta_cols = meta_cols, out_dir = spatial_dir)

  ggsave(file.path(spatial_dir, "ImageFeaturePlot_heat_shock_mod.png"),
         plot = ImageFeaturePlot(srt, features = "signature_1_heat_shock_mod") +
           scale_fill_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                midpoint = 0.3) +
           scale_color_gradient2(low = "steelblue", mid = "white", high = "indianred",
                                 midpoint = 0.3),
         width = 8, height = 7, dpi = 300, limitsize = FALSE)
}
