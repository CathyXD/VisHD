library(Seurat)
library(qs2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
source("~/VisHD/functions.R")

args <- commandArgs(trailingOnly = TRUE)
arg  <- as.numeric(args[1])

paths   <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
path    <- paths[arg]
samples <- sapply(strsplit(paths, "/"), "[", 5)
names(paths) <- samples
sample  <- samples[arg]
cat("Working on", sample, "\n")

setwd(path)

# ── Gene sets ─────────────────────────────────────────────────────────────────
Hall <- readRDS("~/VisHD/Hall.Rds")
C5   <- readRDS("~/VisHD/C5.Rds")
C6   <- readRDS("~/VisHD/C6.Rds")

# ── Load tumour object ────────────────────────────────────────────────────────
srt <- qs_read("tumour/tumour_srt.qs2")
DefaultAssay(srt) <- "SpaNorm"

png_dir <- "tumour/png"
if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)

# Pool CB 0 / CB 1 → "CB"
srt$category_bin <- ifelse(srt$category == "DT", "DT",
                     ifelse(grepl("^CB", srt$category), "CB", NA))

# ── Helpers ───────────────────────────────────────────────────────────────────
run_gsea <- function(deg_df) {
  sig <- deg_df %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
  if (nrow(sig) == 0) return(NULL)
  geneList <- setNames(sig$avg_log2FC, sig$gene)
  lapply(list(Hallmark = Hall, C6 = C6, C5 = C5), function(gs) {
    tryCatch(clusterProfiler::GSEA(geneList, TERM2GENE = gs, verbose = FALSE),
             error = function(e) NULL, warning = function(w) NULL)
  })
}

volcano_plot <- function(deg_df, title) {
  df <- deg_df %>%
    mutate(
      sig   = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25,
      label = ifelse(sig & rank(-abs(avg_log2FC)) <= 20, gene, NA_character_)
    )
  ggplot(df, aes(avg_log2FC, -log10(p_val_adj + 1e-300), colour = sig, label = label)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    scale_colour_manual(values = c("grey70", "firebrick")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = title, x = "avg log2FC  (positive = DT)", y = "-log10(adj p)") +
    theme_classic() + theme(legend.position = "none")
}

save_enrich_plots <- function(enrich_list, prefix) {
  for (nm in names(enrich_list)) {
    res <- enrich_list[[nm]]
    if (is.null(res) || nrow(res@result) == 0) next
    sig_n <- sum(res@result$p.adjust < 0.05)
    if (sig_n == 0) next
    p <- pathwayenrich_plot(top_n = min(10, sig_n), gsea_result = res)
    ggsave(file.path(png_dir, paste0(prefix, "_", nm, "_gsea.pdf")),
           plot = p, width = 8, height = 8)
  }
}

run_comparison <- function(cells, label) {
  n_DT <- sum(srt$category_bin[cells] == "DT", na.rm = TRUE)
  n_CB <- sum(srt$category_bin[cells] == "CB", na.rm = TRUE)
  cat("  DT:", n_DT, "| CB:", n_CB, "\n")
  if (n_DT < 3 || n_CB < 3) { cat("  Skipping: insufficient cells\n"); return(NULL) }

  deg <- FindMarkers(srt,
                     ident.1         = "DT",
                     ident.2         = "CB",
                     group.by        = "category_bin",
                     cells           = cells,
                     test.use        = "MAST",
                     latent.vars     = "nCount_Spatial",
                     logfc.threshold = 0,
                     min.pct         = 0.05)
  deg$gene <- rownames(deg)

  p_vol <- volcano_plot(deg, paste(sample, gsub("_", " ", label)))
  ggsave(file.path(png_dir, paste0(label, "_volcano.png")),
         p_vol, width = 6, height = 5, dpi = 200)

  enrich <- run_gsea(deg)
  save_enrich_plots(enrich, label)

  list(deg = deg, enrich = enrich)
}

future::plan("multisession", workers = 4)

# ══════════════════════════════════════════════════════════════════════════════
# 1.  Overall DT vs CB
# ══════════════════════════════════════════════════════════════════════════════
cat("── Overall DT vs CB ──\n")
cells_all <- colnames(srt)[!is.na(srt$category_bin)]
res_overall <- run_comparison(cells_all, "DT_vs_CB_overall")
if (!is.null(res_overall)) {
  saveRDS(res_overall$deg,    "tumour/deg_DTvsCB_overall.Rds")
  saveRDS(res_overall$enrich, "tumour/enrich_DTvsCB_overall.Rds")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2.  Per-subclone DT vs CB  (only when >1 subclone)
# ══════════════════════════════════════════════════════════════════════════════
subclones <- sort(unique(na.omit(srt$subclone)))
cat("Subclones found:", paste(subclones, collapse = ", "), "\n")

if (length(subclones) > 1) {
  deg_list    <- list()
  enrich_list <- list()
  vol_plots   <- list()

  for (sc in subclones) {
    sc_label <- paste0("DT_vs_CB_subclone_", sc)
    cat("── Subclone", sc, "──\n")
    sc_cells <- colnames(srt)[!is.na(srt$subclone) & srt$subclone == sc &
                                !is.na(srt$category_bin)]
    res <- run_comparison(sc_cells, sc_label)
    if (is.null(res)) next
    deg_list[[sc_label]]    <- res$deg
    enrich_list[[sc_label]] <- res$enrich
    vol_plots[[sc_label]]   <- volcano_plot(res$deg,
                                 paste("DT vs CB — subclone", sc))
  }

  saveRDS(deg_list,    "tumour/deg_DTvsCB_per_subclone.Rds")
  saveRDS(enrich_list, "tumour/enrich_DTvsCB_per_subclone.Rds")

  if (length(vol_plots) > 0) {
    ncols  <- min(3, length(vol_plots))
    p_grid <- wrap_plots(vol_plots, ncol = ncols)
    ggsave(file.path(png_dir, "DT_vs_CB_subclones_volcano_grid.png"),
           p_grid, width = 6 * ncols, height = 5 * ceiling(length(vol_plots) / ncols),
           dpi = 200)
  }
} else {
  cat("Single subclone — skipping per-subclone comparison\n")
}

cat("Done:", sample, "\n")
