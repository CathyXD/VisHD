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

## tumour_anno Rerun DAG (stages 4–10)

The downstream pipeline is re-run whenever the per-sample tumour/normal split is
refined. The corrected split lives in the `tumour_anno` field of the per-sample
object **`LUT-245-XX/tumour_subclone_srt.qs2`** — the canonical entry input for
`4.1`/`4.2`/`4.3` and `9.1`. Cells labelled `"Removed"` are dropped during the
tumour/normal splits before processing.

> **Numbering note:** scripts were renumbered to a contiguous sequence (`0.1 → 10.1`)
> and their output **directory** names were renumbered to match (e.g.
> `4.4.integrate_tumour_anno_srt.R` → `4.4.integrate_tumour_anno/`). The legacy binned
> track and the old `7.1.integration_cells.R` (hardcoded-cluster split) were retired to
> `bin/`, and stages 8–11 were cascaded down one to fill the gap. Paths below are the
> real on-disk dir names. Exceptions that keep their name: `6.2archetype_downstream_tumour/`
> and `9.normalcell_annotation/` (external SCimilarity-hint step, no numbered writer).

**Script roles + I/O (read → write):**

| Script | Mode | Reads | Writes |
|--------|------|-------|--------|
| `4.1.tumour_split.R` | per-sample (1-8) | `tumour_subclone_srt.qs2` (keeps `Tumour`) | `tumour/tumour_srt.qs2` + `tumour/tumour.h5ad` |
| `4.2.normal_split.R` | per-sample (1-8) | `tumour_subclone_srt.qs2` (keeps `Normal`) | `normal/normal_srt.qs2` |
| `4.3.public_signature_exp.R` | per-sample (1-8) | `tumour_subclone_srt.qs2` (drops `Removed`) | `tumour_anno_srt_with_public_signatures.qs2` |
| `4.4.integrate_tumour_anno_srt.R` | run-once | all 8 `tumour/tumour_srt.qs2` | `4.4.integrate_tumour_anno/integrated_pearson_srt.qs2` (+`.h5ad`) |
| `4.5.infercnv_tumour_anno.R` | run-once | `4.4.../integrated_pearson_srt.qs2` | InferCNV (**excluded** from rerun — needs a Normal reference the tumour-only integration lacks) |
| `5.1.DT_deg.R` | per-sample (1-8) | `tumour/tumour_srt.qs2` | `tumour/deg_DTvsCB_*.Rds`, `enrich_*.Rds` |
| `5.2.aggregate_DT_deg.R` | run-once | per-sample 5.1 `.Rds` | `5.2.aggregate_results/` |
| `5.3.integrate_DESeq2_deg.R` | run-once | per-sample `tumour/tumour_srt.qs2` | `5.3.DESeq2_results/` (pseudobulk DESeq2) |
| `5.4.tumour_expression_proportion.R` | run-once | `tumour/tumour_srt.qs2` + 5.3 results | `5.4.expression_proportion/` |
| `6.1.archetypal_analysis_tumour.ipynb` | per-sample (1-8, `SAMPLE_IDX`) | `tumour/tumour.h5ad` | `tumour/archetype_result/` |
| `6.2.archetype_downstream.R` | run-once | per-sample `archetype_result/` | `6.2archetype_downstream_tumour/` |
| `6.3.archetype_module.r` | run-once | 6.2 downstream outputs | archetype/`groupdeg` module gene sets |
| `6.4.signature_analysis.R` | run-once | `tumour/tumour_srt.qs2` + 6.3 groupdeg | `6.4.signature_analysis/` |
| `6.5.integrate_tumour_archetype.ipynb` | run-once | `4.4.../integrated_tumour_anno.h5ad` | `6.5.integrate_tumour_archetype/` |
| `6.6.archetype_visualisation.R` | run-once | `4.4.../integrated_pearson_srt.qs2` + 6.5 csv | `6.5.integrate_tumour_archetype/viz/` |
| `7.1.normalcell_integration_pearson.r` | run-once | all 8 `normal/normal_srt.qs2` | `7.1.normal_cell_integration/integrated_pearson_srt2.qs2` (+`.h5ad`) |
| `7.2.infercnv_check.R` | run-once | `7.1.../integrated_pearson_srt2.qs2` | InferCNV check (read-only) |
| `7.3.normalcell_integration_scvi.r` | run-once | all 8 `normal/normal_srt.qs2` | scVI normal integration (alternative to 7.1) |
| `8.1.scimilarity_check.R` | run-once | `7.1.../integrated_pearson_srt2.qs2` + SCimilarity hints | `8.1.scimilarity_check/normal_srt_annotated.qs2` |
| `8.2.DEG_cluster_annotation.R` | run-once | `8.1.../normal_srt_annotated.qs2` | `8.2.DEG_cluster_annotation/normal_srt_DEcluster.qs2` |
| `8.3.additional_visual.R` | run-once | `8.2.../normal_srt_DEcluster.qs2` | extra cluster visuals (chained after 8.2) |
| `8.4.final_clear_normal_integration.R` | run-once | `8.2.../normal_srt_DEcluster.qs2` | `8.4.final_clear_normal_integration/normal_srt_final_anno.qs2` (+`final_annotation.csv`) |
| `8.5.additional_analysis.R` | run-once | `8.4.../normal_srt_final_anno.qs2` | final-annotation analysis / viz |
| `8.6.final_clear_normal_persample.R` | per-sample (1-8) | `normal/normal_srt.qs2` + 8.4 `final_annotation.csv` | `normal/normal_anno_srt.qs2` |
| `9.1.per_sample_tumour_normal.R` | per-sample (1-8) | `tumour/tumour_srt.qs2` + `normal/normal_anno_srt.qs2` + `tumour_subclone_srt.qs2` | `tumour_normal_anno_srt.qs2` |
| `9.2.per_sample_visual.R` | per-sample (1-8) | `tumour_normal_anno_srt.qs2` | per-sample visuals |
| `9.3.aggreate_cell_composition_analysis.R` | run-once | all 8 `tumour_normal_anno_srt.qs2` | `9.3.aggregate_cell_composition/` |
| `9.4.tumour_normal_integration.R` | run-once | all 8 `tumour_normal_anno_srt.qs2` | `9.4.tumour_normal_integration/integrated_pearson_srt.qs2` |
| `9.5.integrative_infercnv.r` | run-once | `9.4.tumour_normal_integration/...` | InferCNV on the integrated object |
| `10.1.Statial.R` | run-once | all 8 `tumour_normal_anno_srt.qs2` | `combined_sce.qs2` (Statial Kontextual + SpatioMark) |

Stage groups: **7** = normal-cell integration, **8** = normal annotation, **9** =
tumour/normal combined, **10** = spatial statistics.

**Wave structure** (submit with `sbatch --dependency=afterok:<jid>`; jobs within a
wave are mutually independent and run in parallel):

```
Wave 1:  4.1 · 4.2 · 4.3                                   (entry; array 1-8)
Wave 2:  4.4 · 5.1 · 5.3 · 6.1 ← 4.1 ;  7.1 · 7.3 ← 4.2
Wave 3:  5.2 ← 5.1 ; 5.4 ← 5.3 ; 6.2 ← 6.1 ; 6.5 ← 4.4 ; 7.2 · 8.1 ← 7.1
Wave 4:  6.3 ← 6.2 ; 6.6 ← 6.5 ; 8.2 ← 8.1
Wave 5:  6.4 ← 6.3 ; 8.3 · 8.4 ← 8.2
Wave 6:  8.5 · 8.6 ← 8.4
Wave 7:  9.1 ← 4.1 + 8.6                                   (per-sample; array 1-8)
Wave 8:  9.2 · 9.3 · 9.4 · 10.1 ← 9.1
Wave 9:  9.5 ← 9.4
```

Critical paths: normal `4.2 → 7.1 → 8.1 → 8.2 → 8.4 → 8.6 → 9.1`; tumour
`4.1 → 4.4 → 6.5 → 6.6`; combined `9.1 → 9.4 → 9.5`.

**Slurm wrapper resource tiers** (`run_<step>.slurm`, `--partition=work --account=pawsey1172`):
- **Single object / per-sample** (4.1, 4.2, 4.3, 5.1, 6.1, 8.6, 9.1, 9.2): 80G / 1 cpu / 3 h, `--array=1-8`
- **Integrative** (4.4, 5.3, 7.1, 8.1, 9.4): 150G / 4 cpu / 5 h  (6.5 = 200G; 8.1 SCimilarity is the most OOM-prone)
- **Visualization / lightweight** (5.2, 5.4, 6.2, 6.4, 6.6, 8.3, 8.5, 9.3): 50G / 1 cpu / 1 h

To force a clean rerun, delete the existence-guarded integrated caches first
(`4.4.../integrated_pearson_srt.qs2` + `..._tumour_anno.h5ad`, `7.1.../integrated_pearson_srt2.qs2`,
`9.4.tumour_normal_integration/integrated_pearson_srt.qs2`); per-sample `.qs2` are
overwritten unconditionally.

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
