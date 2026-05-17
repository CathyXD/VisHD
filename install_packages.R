#!/usr/bin/env Rscript
# Install all broken/missing packages into the shared R library.
#
# Must be run INSIDE the container so packages compile against the correct R ABI.
#
# On Setonix (interactive):
#   singularity exec -c \
#     -B /tmp/tmp_$USER:/tmp \
#     -B ${MYSCRATCH}:$HOME \
#     spatial_4.5v6.sif \
#     Rscript ~/VisHD/install_packages.R
#
# On Peter Mac (interactive):
#   apptainer exec -c \
#     -B /tmp/tmp_$USER:/tmp \
#     -B /home/sweng:$HOME \
#     spatial_4.5v6.sif \
#     Rscript ~/VisHD/install_packages.R

lib     <- "/scratch/pawsey1172/sweng/R_Library/4.5"
cran    <- "https://cran.r-project.org"

# ==============================================================================
# 1. Remove broken stub directories (dirs missing DESCRIPTION)
# ==============================================================================
cat("\n=== Cleaning broken stubs ===\n")
pkgs  <- list.dirs(lib, full.names = TRUE, recursive = FALSE)
stubs <- pkgs[!file.exists(file.path(pkgs, "DESCRIPTION"))]
cat("Found", length(stubs), "broken stubs — removing...\n")
unlink(stubs, recursive = TRUE)

# Helper: only install packages not already present
install_new <- function(pkgs, lib, ...) {
  have     <- rownames(installed.packages(lib.loc = lib))
  to_get   <- setdiff(pkgs, have)
  if (length(to_get) == 0) { cat("  All already installed.\n"); return(invisible()) }
  cat("  Installing:", paste(to_get, collapse = ", "), "\n")
  install.packages(to_get, lib = lib, ...)
}

# ==============================================================================
# 2. CRAN packages
# ==============================================================================
cat("\n=== CRAN packages ===\n")

cran_pkgs <- c(
  # data wrangling / utilities
  "assertthat", "backports", "brio", "callr", "checkmate", "cli", "clipr",
  "collections", "covr", "credentials", "dbscan", "dendextend", "desc",
  "devtools", "diffobj", "downlit", "ellipsis", "fansi", "fs", "gert",
  "gh", "gitcreds", "ini", "janitor", "jsonlite", "languageserver",
  "lintr", "mclust", "mnormt", "oompaBase", "pbmcapply", "phangorn",
  "pkgbuild", "pkgdown", "pkgload", "praise", "processx", "profvis",
  "proxy", "ps", "psych", "R6", "R.cache", "rcmdcheck", "readxl",
  "renv", "rex", "rlang", "R.methodsS3", "R.oo", "roxygen2", "rstudioapi",
  "R.utils", "rversions", "sessioninfo", "snakecase", "stringi",
  "styler", "testthat", "urlchecker", "usethis", "vdiffr", "waldo",
  "whisker", "xml2", "xmlparsedata", "xopen", "zip",
  # math / stats
  "aricode", "e1071", "GPArotation", "interp", "pROC",
  # graphics / fonts
  "fontBitstreamVera", "fontLiberation", "fontquiver", "gdtools",
  "ggforce", "ggfun", "ggiraph", "ggnewscale", "ggplotify", "ggprism",
  "ggthemes", "gridGraphics", "Hmisc", "htmlTable", "jpeg",
  "latticeExtra", "paletteer", "palr", "pals", "prismatic", "ragg",
  "scatterpie", "scico", "systemfonts", "textshaping", "tweenr", "xgboost",
  # spatial
  "classInt", "s2", "sf", "units", "wk",
  # scRNA / spatial biology
  "aplot", "babelgene", "brew", "ggtangle", "gson", "leidenAlg",
  "msigdbr", "MuDataSeurat", "RApiSerialize", "RcppHungarian",
  "sccore", "scCustomize", "stringfish", "tidytree", "yulab.utils",
  # qs2 (successor to qs, available on CRAN)
  "qs2"
)

install_new(cran_pkgs, lib = lib, repos = cran)

# ==============================================================================
# 3. qs — removed from CRAN; install from source archive
#    Last archived version: 0.27.2  (verify at cran.r-project.org/src/contrib/Archive/qs/)
# ==============================================================================
cat("\n=== qs from CRAN archive ===\n")
if (!"qs" %in% rownames(installed.packages(lib.loc = lib))) {
  qs_url <- "https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.2.tar.gz"
  cat("  Installing qs from:", qs_url, "\n")
  install.packages(qs_url, lib = lib, repos = NULL, type = "source")
} else {
  cat("  qs already installed.\n")
}

# ==============================================================================
# 4. Bioconductor packages
# ==============================================================================
cat("\n=== Bioconductor packages ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = lib, repos = cran)

bioc_pkgs <- c(
  "aroma.light", "clusterProfiler", "DOSE", "EDASeq", "enrichplot",
  "fgsea", "genefilter", "GO.db", "GOSemSim", "ggtree", "hwriter",
  "qvalue", "RUVSeq", "ShortRead", "treeio"
)

have <- rownames(installed.packages(lib.loc = lib))
bioc_to_get <- setdiff(bioc_pkgs, have)
if (length(bioc_to_get) > 0) {
  cat("  Installing:", paste(bioc_to_get, collapse = ", "), "\n")
  BiocManager::install(bioc_to_get, lib = lib, update = FALSE, ask = FALSE)
} else {
  cat("  All already installed.\n")
}

# ==============================================================================
# 5. GitHub packages
# ==============================================================================
cat("\n=== GitHub packages ===\n")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", lib = lib, repos = cran)

github_pkgs <- list(
  fastCNV    = "must-bioinfo/fastCNV",
  harrypotter = "aljrico/harrypotter"
)

have <- rownames(installed.packages(lib.loc = lib))
for (pkg in names(github_pkgs)) {
  if (!pkg %in% have) {
    cat("  Installing", pkg, "from", github_pkgs[[pkg]], "\n")
    remotes::install_github(github_pkgs[[pkg]], lib = lib, upgrade = "never")
  } else {
    cat(" ", pkg, "already installed.\n")
  }
}

# ==============================================================================
# 6. Unknown source — mcprogress, scPearsonPCA
#    These were not on CRAN/Bioc at time of writing; check GitHub manually if needed.
# ==============================================================================
cat("\n=== Attempting mcprogress / scPearsonPCA from CRAN (may fail) ===\n")
tryCatch(install_new(c("mcprogress", "scPearsonPCA"), lib = lib, repos = cran),
         error = function(e) cat("  Skipped (not on CRAN):", conditionMessage(e), "\n"))

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Installation complete ===\n")
have <- rownames(installed.packages(lib.loc = lib))
cat("Packages now in lib:", length(have), "\n")
