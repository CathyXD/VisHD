#!/usr/bin/env Rscript
# 10.4.aggregate_tumour_border.R   (run-once)
# Aggregates the per-sample cell-position tables from 10.3.SPIAT_tumour_border.R
# into cross-sample final_annotation composition per tumour-border Structure
# category (Internal margin / External margin / Inside / Outside), split by
# category_bin (CB vs DT).
#
# Reads:  ~/VisHD/10.3.tumour_border/per_sample_tables/<slide>_cell_position.csv
# Writes: ~/VisHD/10.4.tumour_border_aggregate/
#           combined_tables/all_samples_cell_position.csv
#           aggregate_barplots/persample_barplot.png (+.csv)
#           aggregate_barplots/average_barplot.png   (+.csv)
#           spatial_plots/all_samples_structure_spatial.png
#
#   Rscript 10.4.aggregate_tumour_border.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(pals)
})

out_dir    <- "~/VisHD/10.4.tumour_border_aggregate"
tbl_dir    <- file.path(out_dir, "combined_tables")
barplt_dir <- file.path(out_dir, "aggregate_barplots")
spatial_dir <- file.path(out_dir, "spatial_plots")
for (d in c(tbl_dir, barplt_dir, spatial_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

src_dir <- "~/VisHD/10.3.tumour_border/per_sample_tables"
paths   <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
samples <- basename(paths)

group_pal    <- c("Neg" = "lightgrey", "G1" = "red", "G2" = "gold", "G3" = "royalblue",
                  "G2/G1" = "orange", "G3/G1" = "purple", "G3/G2" = "green",
                  "G3/G2/G1" = "darkgrey")


read_one <- function(s) {
  f <- file.path(src_dir, paste0(s, "_cell_position.csv"))
  if (!file.exists(f)) {
    message(s, ": no cell-position table found, skipping")
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}

all_pos <- do.call(rbind, lapply(samples, read_one))
write.csv(all_pos, file.path(tbl_dir, "all_samples_cell_position.csv"), row.names = FALSE)
cat("Combined table:", nrow(all_pos), "cells across", length(unique(all_pos$slide)), "samples\n")

cols = c("#8DD3C7","#FFFFB3","#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
         "#FCCDE5", "#fd64b3", "#888787")
anno_pal <- setNames(cols[seq_along(unique(all_pos$cell_type))], unique(all_pos$cell_type))
anno_pal <- c(anno_pal, group_pal)

struct_lvls <- c("Inside",  "Infiltrated.CoI", "Internal.margin", "Internal.margin.CoI", "Border", "External.margin","External.margin.CoI", "Stromal.CoI" ,  "Outside")
all_pos$Structure <- factor(all_pos$Structure,
                            levels = intersect(struct_lvls, unique(all_pos$Structure)))

# SPIAT's plot_cell_categories() default Structure palette
# (https://trigosteam.github.io/SPIAT/articles/tissue-structure.html)
struct_pal <- c("Border" = "black", "Inside" = "pink", "Infiltrated.CoI" = "purple",
                "Outside" = "yellow", "Stromal.CoI" = "orange",
                "Internal.margin" = "lightgreen", "Internal.margin.CoI" = "darkgreen",
                "External.margin" = "lightblue", "External.margin.CoI" = "blue")

comp <- all_pos %>%
  filter(!is.na(Structure), !is.na(final_annotation), !is.na(category_bin)) %>%
  count(slide, category_bin, Structure, final_annotation) %>%
  group_by(slide, category_bin, Structure) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


fa_pal  <- c(group_pal, anno_pal)

# ── Plot 1: per-sample barplot, aggregated in one plot, faceted by category_bin ──
p1 <- ggplot(comp, aes(category_bin, prop, fill = final_annotation)) +
  geom_col(position = "fill") +
  facet_grid(Structure ~ slide) +
  scale_fill_manual(values = fa_pal, name = "final_annotation") +
  scale_y_continuous(labels = percent) +
  labs(y = "proportion", title = "Cell-type composition by tumour-border structure") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(barplt_dir, "persample_barplot.png"), p1,
       width = max(10, length(samples) * 1.5), height = 12, dpi = 200, limitsize = FALSE)
write.csv(comp, file.path(barplt_dir, "persample_barplot.csv"), row.names = FALSE)

# ── Plot 2: across-sample average barplot, faceted by category_bin ───────────────
avg <- comp %>%
  group_by(category_bin, Structure, final_annotation) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

p2 <- ggplot(avg, aes(category_bin, mean_prop, fill = final_annotation)) +
  geom_col(position = "stack") +
  facet_wrap(~ Structure) +
  scale_fill_manual(values = fa_pal, name = "final_annotation") +
  scale_y_continuous(labels = percent) +
  labs(y = "mean proportion across samples",
       title = "Mean cell-type composition by tumour-border structure") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(barplt_dir, "average_barplot.png"), p2,
       width = 10, height = 6, dpi = 150, limitsize = FALSE)
write.csv(avg, file.path(barplt_dir, "average_barplot.csv"), row.names = FALSE)

# ── Plot 3: aggregated spatial plot, all samples faceted, coloured by Structure ──
p3 <- ggplot(all_pos[!is.na(all_pos$Structure), ], aes(x, y, colour = Structure)) +
  geom_point(size = 0.15) +
  scale_y_reverse() +
  scale_colour_manual(values = struct_pal, drop = FALSE, name = "Structure") +
  facet_wrap(~ slide, scales = "free", ncol = 4) +
  theme_classic(base_size = 18) + theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(spatial_dir, "all_samples_structure_spatial.png"), p3,
       width = 20, height = 10, dpi = 300, limitsize = FALSE)

cat("\n==================== 10.4.aggregate_tumour_border done ====================\n")
