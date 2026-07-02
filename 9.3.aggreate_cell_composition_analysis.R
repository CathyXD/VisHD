#!/usr/bin/env Rscript
# 9.3.aggreate_cell_composition_analysis.R   (run-once, all 8 samples)
# Aggregate the per-sample composition bar-plot tables written by
# 9.1.per_sample_tumour_normal.R (LUT-245-XX/final_png/barplots/
# barplot_<group_var>_by_<fill_var>.csv, columns: <group_var>,<fill_var>,Freq)
# into cross-sample visualizations. For every (group_var x fill_var) combination:
#
#   1. combined_tables/  all slides concatenated + within-(slide,group) proportion
#   2. barplots/         composition (stacked proportion) faceted by slide; the
#                        facet label shows that slide's total cell number
#   3. heatmaps/         x = "<slide>-<group level>", y = fill category,
#                        fill = proportion
#   4. boxplots/         per-sample composition proportions, x = group level,
#                        one panel per fill category, each point = a sample,
#                        Wilcoxon between group levels (annotated on 2-level vars)
#   5. wilcoxon/         Wilcoxon rank-sum results per combo + a master table
#
# Proportion = Freq / sum(Freq) within each (slide, group level) — the same
# composition the 10.0 bar plots show with position = "fill".
#
# NOTE: cluster group-variables (banksy_clusters / pearson_clusters /
# SpaNorm_snn_res.1 / seurat_clusters) use per-slide cluster IDs — cluster "1"
# in one slide is unrelated to "1" in another, so their shared-axis heatmaps and
# cross-sample boxplots are NOT biologically aligned (included by request).
#
#   Rscript 9.3.aggreate_cell_composition_analysis.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(pals)
})

setwd("~/VisHD")
set.seed(1)                                     # reproducible jitter

out_dir <- "~/VisHD/10.0.1.aggregate_cell_composition"
tbl_dir <- file.path(out_dir, "combined_tables")
bar_dir <- file.path(out_dir, "barplots")
hm_dir  <- file.path(out_dir, "heatmaps")
box_dir <- file.path(out_dir, "boxplots")
wx_dir  <- file.path(out_dir, "wilcoxon")
for (d in c(tbl_dir, bar_dir, hm_dir, box_dir, wx_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ───────────────────────────────────────────────────────────────────
# Natural ordering: numeric-like levels sorted numerically, else alphabetically.
order_levels <- function(x) {
  ux <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", ux))) ux[order(as.integer(ux))] else sort(ux)
}

# Qualitative palette sized to n categories (polychrome up to 36, then ramped).
make_pal <- function(levs) {
  n <- length(levs)
  base <- as.vector(pals::polychrome(36))
  cols <- if (n <= 36) base[seq_len(max(n, 1))] else
          grDevices::colorRampPalette(base)(n)
  setNames(cols[seq_len(n)], levs)
}

# ── Fill palettes reused verbatim from 9.1.per_sample_tumour_normal.R ─────────
# module-group colours (module_anno fill).
group_pal <- c("Neg" = "lightblue", "G1" = "red", "G2" = "gold", "G3" = "royalblue",
               "G1/G2" = "orange", "G1/G3" = "purple", "G2/G3" = "green",
               "G1/G2/G3" = "grey")
canon <- function(x) vapply(strsplit(x, "/"),
                            function(p) paste(sort(p), collapse = "/"), character(1))

# cell_type fill: Set3 by sorted level, Tumour -> red/gold/orange (10.0 lines 253-255).
build_ct_pal <- function(levs) {
  levs <- sort(unique(as.character(levs)))
  pal  <- setNames(as.vector(pals::brewer.set3(max(3, length(levs))))[seq_along(levs)], levs)
  tum  <- grep("Tumour", names(pal))
  if (length(tum)) pal[tum] <- c("red", "gold", "orange")[seq_along(tum)]
  pal
}
# module_anno fill: group_pal by canonical name + Normal = lightpink (10.0 lines 229-231).
build_mg_pal <- function(levs) {
  levs <- unique(as.character(levs))
  pal  <- setNames(group_pal[canon(levs)], levs)
  pal["Normal"] <- "lightpink"
  pal
}

# All pairwise Wilcoxon rank-sum tests of `prop` across group levels, for one
# fill category. Returns a tidy data.frame (NA stats on failure / low n).
wilcox_pairs <- function(sub, combo, group_var, fill_var, fill_level) {
  gl <- order_levels(sub$group_level[!is.na(sub$prop)])
  if (length(gl) < 2) return(NULL)
  rows <- lapply(combn(gl, 2, simplify = FALSE), function(pr) {
    x  <- sub$prop[as.character(sub$group_level) == pr[1] & !is.na(sub$prop)]
    y  <- sub$prop[as.character(sub$group_level) == pr[2] & !is.na(sub$prop)]
    tt <- tryCatch(wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
    data.frame(combo = combo, group_var = group_var, fill_var = fill_var,
               fill_level = fill_level, group1 = pr[1], group2 = pr[2],
               n1 = length(x), n2 = length(y),
               statistic = if (is.null(tt)) NA_real_ else unname(tt$statistic),
               p_value   = if (is.null(tt)) NA_real_ else tt$p.value,
               method    = "wilcoxon rank-sum",
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

# ── 1. Discover per-sample composition CSVs and group them by combo ────────────
paths  <- system("realpath ~/VisHD/LUT-245-*/", intern = TRUE)
slides <- basename(paths)
cat("Slides:", paste(slides, collapse = ", "), "\n")

bar_files <- unlist(lapply(paths, function(p)
  list.files(file.path(p, "final_png", "barplots"),
             pattern = "^barplot_.*_by_(cell_type|module_anno)\\.csv$",
             full.names = TRUE)))
if (length(bar_files) == 0) stop("No barplot_*_by_*.csv files found under final_png/barplots/")

combos <- split(bar_files, basename(bar_files))
cat("Found", length(bar_files), "CSVs across", length(combos), "combos\n")

slide_pal  <- make_pal(slides)
all_long   <- list()   # for a master combined table
all_wilcox <- list()   # for a master wilcoxon table

# Build each fill palette ONCE from the global union of levels, so a given cell
# type / module group keeps one colour across every combo (10.0's palettes).
fill_levels <- list()
for (f in bar_files) {
  d  <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  fv <- names(d)[2]
  fill_levels[[fv]] <- union(fill_levels[[fv]], as.character(d[[2]]))
}
fill_pals <- list()
if (!is.null(fill_levels[["cell_type"]]))
  fill_pals[["cell_type"]]   <- build_ct_pal(fill_levels[["cell_type"]])
if (!is.null(fill_levels[["module_anno"]]))
  fill_pals[["module_anno"]] <- build_mg_pal(fill_levels[["module_anno"]])

# ── 2-6. Process each (group_var x fill_var) combo ─────────────────────────────
for (combo in names(combos)) {
  files <- combos[[combo]]
  name  <- sub("\\.csv$", "", combo)
  hdr       <- names(read.csv(files[1], check.names = FALSE, nrows = 1))
  group_var <- hdr[1]; fill_var <- hdr[2]
  cat("\n== ", combo, " (", group_var, " x ", fill_var, ", ",
      length(files), " slides) ==\n", sep = "")

  # 2. Load + combine + within-(slide,group) proportion
  df <- do.call(rbind, lapply(files, function(f) {
    d     <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    slide <- basename(dirname(dirname(dirname(f))))   # LUT-XX/final_png/barplots/f
    data.frame(slide       = slide,
               group_level = as.character(d[[1]]),
               fill_level  = as.character(d[[2]]),
               Freq        = as.numeric(d[[3]]),
               stringsAsFactors = FALSE)
  }))
  df <- df %>%
    group_by(slide, group_level) %>%
    mutate(sg_total = sum(Freq), prop = ifelse(sg_total > 0, Freq / sg_total, NA_real_)) %>%
    ungroup() %>%
    group_by(slide) %>%
    mutate(slide_total = sum(Freq)) %>%
    ungroup()
  df$group_level <- factor(df$group_level, levels = order_levels(df$group_level))
  df$fill_level  <- factor(df$fill_level,  levels = sort(unique(df$fill_level)))

  out_tbl <- df %>% mutate(group_level = as.character(group_level),
                           fill_level  = as.character(fill_level),
                           group_var   = group_var, fill_var = fill_var)
  write.csv(out_tbl, file.path(tbl_dir, paste0(name, ".csv")), row.names = FALSE)
  all_long[[combo]] <- out_tbl

  fill_pal <- if (!is.null(fill_pals[[fill_var]])) fill_pals[[fill_var]] else
              make_pal(levels(df$fill_level))
  n_g <- nlevels(df$group_level); n_f <- nlevels(df$fill_level); n_s <- length(unique(df$slide))

  # 3. Composition bar plots faceted by slide (facet label carries slide total)
  slide_tot <- df %>% distinct(slide, slide_total)
  lab_map   <- setNames(sprintf("%s\n(n=%s)", slide_tot$slide,
                                formatC(slide_tot$slide_total, format = "d", big.mark = ",")),
                        slide_tot$slide)
  ncf_bar <- min(3, n_s)
  p_bar <- ggplot(df, aes(group_level, prop, fill = fill_level)) +
    geom_col() +
    facet_wrap(~ slide, ncol = ncf_bar, labeller = labeller(slide = lab_map)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = fill_pal, na.value = "grey80", name = fill_var) +
    labs(title = sprintf("%s composition per %s, by slide", fill_var, group_var),
         x = group_var, y = "proportion") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(bar_dir, paste0(name, ".png")), p_bar,
         width  = max(12, ncf_bar * (n_g * 0.25 + 2.5)),
         height = max(8,  ceiling(n_s / ncf_bar) * 3.5),
         dpi = 150, limitsize = FALSE)

  # 4. Heatmap: x = "<slide>-<group level>", y = fill category, fill = proportion
  hd <- df %>% arrange(slide, group_level) %>%
    mutate(xkey = paste(slide, group_level, sep = "-"))
  hd$xkey <- factor(hd$xkey, levels = unique(hd$xkey))
  nx <- nlevels(hd$xkey)
  p_hm <- ggplot(hd, aes(xkey, fill_level, fill = prop)) +
    geom_tile(colour = "grey90") +
    scale_fill_viridis_c(labels = scales::percent, name = "proportion", limits = c(0, 1)) +
    labs(title = sprintf("%s composition: %s (proportion within each slide x %s)",
                         fill_var, group_var, group_var),
         x = paste0("slide - ", group_var), y = fill_var) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(hm_dir, paste0(name, ".png")), p_hm,
         width  = max(8, nx * 0.18 + 3),
         height = max(5, n_f * 0.24 + 2),
         dpi = 150, limitsize = FALSE)

  # 5-6. Boxplots (point = sample) + Wilcoxon between group levels per fill
  wx <- do.call(rbind, lapply(levels(df$fill_level), function(fl)
    wilcox_pairs(df[df$fill_level == fl, , drop = FALSE], combo, group_var, fill_var, fl)))
  if (!is.null(wx)) {
    write.csv(wx, file.path(wx_dir, paste0(name, "_wilcox.csv")), row.names = FALSE)
    all_wilcox[[combo]] <- wx
  }

  ncf_box <- min(4, n_f)
  p_box <- ggplot(df, aes(group_level, prop)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92") +
    geom_jitter(aes(colour = slide), width = 0.15, height = 0, size = 1.4) +
    facet_wrap(~ fill_level, ncol = ncf_box, scales = "free_y") +
    scale_colour_manual(values = slide_pal, name = "slide") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = sprintf("%s composition proportions across samples, by %s",
                         fill_var, group_var),
         subtitle = "each point = one slide; Wilcoxon between group levels",
         x = group_var, y = "proportion (per sample)") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # On-plot p-value only when the group variable has exactly 2 levels (readable);
  # multi-level combos report every pair in the saved wilcoxon CSV instead.
  if (n_g == 2 && !is.null(wx)) {
    ymax  <- df %>% group_by(fill_level) %>%
      summarise(y = max(prop, na.rm = TRUE), .groups = "drop")
    annot <- wx %>% group_by(fill_level) %>% slice(1) %>% ungroup() %>%
      mutate(label = ifelse(is.na(p_value), "p=NA", paste0("p=", signif(p_value, 2))))
    annot$fill_level <- factor(annot$fill_level, levels = levels(df$fill_level))
    annot <- left_join(annot, ymax, by = "fill_level")
    p_box <- p_box +
      geom_text(data = annot, aes(x = 1.5, y = y * 1.08, label = label),
                inherit.aes = FALSE, size = 2.6)
  }
  ggsave(file.path(box_dir, paste0(name, ".png")), p_box,
         width  = max(8, ncf_box * 3),
         height = max(5, ceiling(n_f / ncf_box) * 3),
         dpi = 150, limitsize = FALSE)
}

# ── Master combined tables ─────────────────────────────────────────────────────
write.csv(do.call(rbind, all_long),
          file.path(tbl_dir, "all_composition_long.csv"), row.names = FALSE)
if (length(all_wilcox) > 0)
  write.csv(do.call(rbind, all_wilcox),
            file.path(wx_dir, "all_wilcoxon_results.csv"), row.names = FALSE)

cat("\n==================== 10.0.1 aggregate composition done ====================\n")
cat("Outputs under", out_dir, "\n")
