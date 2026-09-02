#!/usr/bin/env Rscript
# 10.6.aggregate_colocalization.R   (run-once)
# Aggregates the per-sample Cross-K curves + AUC from
# 10.5.SPIAT_colocalization.R into one grid figure: rows = cell type,
# columns = sample, each panel showing the observed vs. theoretical Cross-K
# curve for that (sample, cell type) pair with its Cross-K AUC annotated.
#
# Reads:  ~/VisHD/10.5.tumour_colocalization/per_sample_tables/<slide>_crossK_summary.csv
#         ~/VisHD/10.5.tumour_colocalization/per_sample_crossK/<slide>_<cell_type>_crossK.csv
# Writes: ~/VisHD/10.6.aggregate_colocalization/
#           combined_tables/all_samples_crossK_summary.csv
#           combined_tables/all_samples_crossK_curves.csv
#           crossK_grid.png
#
#   Rscript 10.6.aggregate_colocalization.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

src_summary_dir <- "~/VisHD/10.5.tumour_colocalization/per_sample_tables"
src_crossk_dir  <- "~/VisHD/10.5.tumour_colocalization/per_sample_crossK"

out_dir <- "~/VisHD/10.6.aggregate_colocalization"
tbl_dir <- file.path(out_dir, "combined_tables")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

paths   <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
samples <- basename(paths)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Combine per-sample Cross-K AUC summaries
# ══════════════════════════════════════════════════════════════════════════════
read_summary <- function(s) {
  f <- file.path(src_summary_dir, paste0(s, "_crossK_summary.csv"))
  if (!file.exists(f)) {
    message(s, ": no crossK_summary.csv found, skipping")
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}
all_auc <- do.call(rbind, lapply(samples, read_summary))
if (is.null(all_auc) || !nrow(all_auc))
  stop("No crossK_summary.csv files found under ", src_summary_dir,
       " — run 10.5.SPIAT_colocalization.R first.")
write.csv(all_auc, file.path(tbl_dir, "all_samples_crossK_summary.csv"), row.names = FALSE)
cat("Combined AUC summary:", nrow(all_auc), "(slide, cell_type) pairs across",
    length(unique(all_auc$slide)), "samples\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Combine the per-pair Cross-K curves (r/theo/border), joined to the AUC table
# ══════════════════════════════════════════════════════════════════════════════
read_curve <- function(slide, cell_type) {
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", cell_type)
  f <- file.path(src_crossk_dir, sprintf("%s_%s_crossK.csv", slide, ct_safe))
  if (!file.exists(f)) {
    message(slide, "/", cell_type, ": no Cross-K curve CSV found, skipping")
    return(NULL)
  }
  cbind(slide = slide, cell_type = cell_type, read.csv(f, stringsAsFactors = FALSE))
}
all_curves <- do.call(rbind, Map(read_curve, all_auc$slide, all_auc$cell_type))
if (is.null(all_curves) || !nrow(all_curves))
  stop("No Cross-K curve CSVs found under ", src_crossk_dir)
write.csv(all_curves, file.path(tbl_dir, "all_samples_crossK_curves.csv"), row.names = FALSE)
cat("Combined curves:", nrow(all_curves), "rows\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Grid figure: rows = cell type, columns = sample; AUC annotated per panel
# ══════════════════════════════════════════════════════════════════════════════
slide_lvls <- sort(unique(all_curves$slide))
ct_lvls    <- sort(unique(all_curves$cell_type))
all_curves$slide     <- factor(all_curves$slide, levels = slide_lvls)
all_curves$cell_type <- factor(all_curves$cell_type, levels = ct_lvls)

auc_labels <- all_auc %>%
  mutate(slide = factor(slide, levels = slide_lvls),
         cell_type = factor(cell_type, levels = ct_lvls),
         label = sprintf("AUC=%.3f", crossK_AUC))

p <- ggplot(all_curves, aes(x = r)) +
  geom_line(aes(y = theo, linetype = "Theoretical (CSR)"), colour = "grey40") +
  geom_line(aes(y = border, linetype = "Observed"), colour = "firebrick") +
  geom_text(data = auc_labels, aes(label = label), x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.4, size = 3, inherit.aes = FALSE) +
  facet_grid(cell_type ~ slide, scales = "free_y") +
  scale_linetype_manual(values = c("Theoretical (CSR)" = "dashed", "Observed" = "solid"), name = NULL) +
  labs(x = "r (px)", y = "Cross-K(r)",
       title = "Tumour vs. normal cell-type colocalization — Cross-K function",
       subtitle = "AUC > 0: clustering; AUC < 0: dispersion") +
  theme_bw(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")

ggsave(file.path(out_dir, "crossK_grid.png"), p,
       width = max(10, length(slide_lvls) * 2.4), height = max(8, length(ct_lvls) * 1.9),
       dpi = 200, limitsize = FALSE, bg = "white")

cat("\n==================== 10.6.aggregate_colocalization done ====================\n")
