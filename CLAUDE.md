# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Spatial transcriptomics analysis pipeline for Visium HD (single-cell resolution) data, focused on prostate cancer. Processes 10X Genomics Visium HD output to identify tumour vs normal cells, infer copy-number variations (CNVs), and characterize tumour subclones across 8 samples (`LUT-245-07` through `LUT-245-20`).

## Running Scripts

Scripts are designed to run as **Slurm job arrays** on HPC. The array index (`$SLURM_ARRAY_TASK_ID`) selects the sample to process (1–8). There are two HPC environments:

**Setonix (Pawsey)** — uses Singularity container:
```bash
sbatch runarray.slurm       # or runarray2.slurm
```
Key settings: `--account=pawsey1172`, `--partition=work`, `dir=$MYSCRATCH`, image `docker://swweng/spatial:4.5v6`

**Peter Mac** — uses Apptainer or native R:
```bash
sbatch runarrayPM.slurm          # Apptainer container
sbatch runarrayPM_simple.slurm   # Native R via module load ood_r/4.5.0
```

To run a single script manually for one sample (e.g., sample index 1):
```bash
Rscript ~/VisHD/0.generate_srt.R 1
```

Slurm logs go to `slurmout/`.

## Pipeline Steps

Scripts are numbered by their order of execution:

| Script | Purpose |
|--------|---------|
| `0.generate_srt.R` | Load raw 10X HDF5 + GeoJSON, build Seurat object with spatial coordinates and cell/nucleus areas; also generates 16µm binned object |
| `0.1.category_annotation.R` | Initial category annotation |
| `1.banksy.R` | BANKSY spatial clustering at λ=0.2 and λ=0.8 on SpaNorm-normalized data |
| `2.ATAClone.R` | ATAClone CNV analysis on 16µm binned data using 10 Mbp genomic windows |
| `2.infercnv_test.R` / `2.infercnv_ATAClone.R` | InferCNV CNV inference using tumour/normal cell labels |
| `3.2.tumour_normal_expression.R` | Expression visualizations for tumour vs normal |
| `4.1refine_subclones.R` | **Manual step** — hardcoded cluster→tumour/normal/subclone mappings per sample; run interactively |
| `4.2infercnv_subclone.R` | InferCNV on refined subclones |
| `4.convert_h5ad.R` | Export Seurat to `.h5ad` (AnnData) via MuDataSeurat |
| `5.1.tumour_normal_split.R` | Subset tumour cells, re-run SpaNorm, BANKSY, module scoring |
| `5.DT_deg.Rmd` | Differential expression and pathway enrichment (interactive Rmd) |
| `6.celltype_annotation.R` | Cell type annotation |
| `7.1.integration_cells.R` | Merge raw counts from all 8 slides, run scPearsonPCA (with and without batch correction by `slide`), compute archetype + tumour/normal module scores, export to `.h5ad`; outputs to `integration/` |
| `LUT-245-XX/4.refine_cnv.Rmd` | Per-sample interactive CNV refinement notebook |

## Data Layout

```
~/VisHD/
├── raw/LUT-245-XX/outs/          # 10X pipeline outputs (raw data, gitignored)
│   └── segmented_outputs/
│       ├── filtered_feature_cell_matrix.h5
│       ├── graphclust_annotated_cell_segmentations.geojson
│       └── nucleus_segmentations.geojson
├── LUT-245-XX/                   # Per-sample analysis outputs
│   ├── raw_srt.qs                # Raw Seurat object
│   ├── spanorm_srt.qs            # After SpaNorm normalization
│   ├── banksy_srt.qs/qs2         # After BANKSY clustering
│   ├── tumour_anno_srt.qs2       # After tumour/normal annotation
│   ├── bined_ouput/              # 16µm binned Seurat objects
│   │   ├── srt.qs
│   │   └── srt_infercnv.qs2
│   ├── tumour/                   # Tumour-only re-analysis
│   │   ├── tumour_srt.qs2
│   │   └── vishd_counts.h5ad
│   └── category.csv              # CB (CellBiome) binary cell classification
└── integration/                  # Cross-sample integration outputs (7.1.integration_cells.R)
    ├── integrated_pearson_srt.qs2  # Merged Seurat with pearsonpca/umap reductions
    └── integrated_cells.h5ad       # AnnData export with spatial obsm
```

`$HOME` inside the container maps to `$MYSCRATCH` on Setonix.

## Key Concepts

**Cell classification metadata fields:**
- `SR_Cluster`: Initial graph cluster from 10X pipeline
- `category`: CB (CellBiome) classification (`CB 0`, `CB 1`) or `DT` (unclassified cells)
- `tumour_anno`: `Tumour` / `Normal` / `Removed` (from `4.1refine_subclones.R`)
- `ATAClone_cluster`: CNV-based clonal cluster labels
- `subclone`: Manually curated integer subclone assignment
- `slide`: Sample name (e.g. `LUT-245-07`); added during cross-sample merge in `7.1.integration_cells.R`
- `x_centroid` / `y_centroid`: Spatial centroid coordinates (pixels) extracted via `GetTissueCoordinates(srt, which = "centroids")` and stored in `@meta.data` so they survive merge/export
- `pearson_clusters`: Leiden clusters from unbatch-corrected scPearsonPCA graph
- `pearson_clusters_batch`: Leiden clusters from batch-corrected scPearsonPCA graph

**Two data resolutions run in parallel:**
- **Cell-level** (`raw_srt.qs`): Single-cell segmented Visium HD data
- **Binned** (`bined_ouput/`): 16µm bins (used for CNV analysis — more counts per unit)

Annotations are transferred from bins to cells using `transfer_visiumhd_to_cells()` (k-NN in spatial coordinates).

**Tumour markers used throughout:** `AR`, `FOLH1`, `KLK2`, `KLK3`, `KLK4`, `TMPRSS2`, `NKX3-1`, `HOXB13`, `TRPM8`

**Tumour archetypes** (from `clean_module.Rds`): `AR`, `Inflammation`, `NE1`, `NE2`, `Cycling`, `Glycolysis`

## Key Functions (`functions.R`)

| Function | Description |
|----------|-------------|
| `do.spanorm(srt)` | Runs SpaNorm normalization → SVG detection → PCA → UMAP → Leiden clustering (res=1); returns modified `srt` |
| `do.banksy(srt)` | Wraps BANKSY spatial clustering (uses SpaNorm assay, λ=0.2 and 0.8) |
| `spatial_plot(srt, outdir, name)` | Saves standard diagnostic PNGs (ImageDimPlot, FeaturePlot, VlnPlot) to `outdir` |
| `transfer_visiumhd_to_cells(srt_cell, srt_bin, ...)` | Transfers metadata from binned to cell-level object via FNN k-NN |
| `binarise_expression(expr, ...)` | Gaussian mixture model binarization; returns 0/1 integer vector |
| `srt2anndata(srt, ..., save_name)` | Full Seurat → AnnData export including SpaNorm data, PCA, UMAP, spatial coords |
| `filter_artefacts_knn(srt, ...)` | Removes isolated cells by k-NN spatial density |
| `plot_cnv_heatmap(mat, labels, ...)` | ComplexHeatmap CNV matrix with cluster row annotations |
| `pathwayenrich_plot(top_n, gsea_result, ...)` | GSEA visualization using Hallmark/C5/C6 gene sets |

## Custom R Library

Several packages are installed to `~/R_Library/4.5` (not the system library). Always load these with:
```r
library(SpaNorm,      lib.loc = "~/R_Library/4.5")
library(leidenbase,   lib.loc = "~/R_Library/4.5")
library(UCell,        lib.loc = "~/R_Library/4.5")
library(anndataR,     lib.loc = "~/R_Library/4.5")
library(scPearsonPCA, lib.loc = "~/R_Library/4.5")
library(qs,           lib.loc = "~/R_Library/4.5")
```

Key `scPearsonPCA` functions: `sparse_quasipoisson_pca_seurat()` (no batch), `sparse_quasipoisson_pca_seurat_batch()` (batch by slide), `make_umap()`, `gene_frequency()`.

## File Formats

- `.qs` / `.qs2`: Fast R object serialization (use `qread()`/`qsave()` for qs, `qs_read()`/`qs_save()` for qs2)
- `.h5ad`: AnnData format for downstream Python/scanpy analysis
- `.Rds`: Standard R serialization for smaller objects (gene signatures, annotations)
- `.Rmd` / `.nb.html`: Per-sample interactive analysis notebooks (CNV refinement)

## Reference Data Files

| File | Contents |
|------|----------|
| `clean_module.Rds` | Tumour archetype gene sets (AR, Inflammation, NE1, NE2, Cycling, Glycolysis) |
| `proberef.Rds` | Visium probe genomic coordinates (chr, start, end per gene) |
| `gene_ord2.Rds` | Gene order file for InferCNV |
| `gene_ord.txt` | Gene ordering reference |
| `infercnv_10Mbp_genomeref.Rds` | 10 Mbp genomic window reference for ATAClone |
| `stable_counts_ref.Rds` | Stable count reference |
| `Hall.Rds`, `C5.Rds`, `C6.Rds` | MSigDB gene sets (Hallmark, C5 ontology, C6 oncogenic) |
| `normal_markers.R` | Loads `all_marker` — named list of normal cell type gene signatures |

## Coding Guidelines (Karpathy)

**1. Think Before Coding** — State assumptions explicitly. If multiple interpretations exist, present them. If something is unclear, ask rather than guess.

**2. Simplicity First** — Minimum code that solves the problem. No speculative features, abstractions, or error handling for impossible scenarios.

**3. Surgical Changes** — Touch only what the request requires. Don't improve adjacent code. Remove only imports/variables YOUR changes made unused; leave pre-existing dead code unless asked.

**4. Goal-Driven Execution** — Define verifiable success criteria before starting. For multi-step tasks, state a brief plan with a check per step.
