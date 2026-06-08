"""SCimilarity cell-type annotation for the normal-cell Pearson integration.

Reads ~/VisHD/8.1.normal_cell_integration/integrated_normal_cells.h5ad,
runs Genentech's SCimilarity foundation model against the prepackaged
23 M-cell reference, relabels low-confidence calls as "Unknown", and writes
the annotated AnnData + per-cell CSV + UMAP PNG into ~/VisHD/9.normalcell_annotation.

Tutorial:
  https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
"""

from pathlib import Path
import os
import sys

import matplotlib
matplotlib.use("Agg")

import scanpy as sc
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts


# ── Paths ──────────────────────────────────────────────────────────────────
# Use $MYSCRATCH as the project root — works both natively on Pawsey and
# under the singularity bind-mount that maps $MYSCRATCH to $HOME.
ROOT    = Path(os.environ.get("MYSCRATCH", str(Path.home()))) / "VisHD"
IN_H5AD = ROOT / "1.integrate_raw_cell/normal.h5ad"
MODEL   = Path("/scratch/pawsey1172/sweng/VisHD/SCimilarity/model_v1.1")
OUT     = ROOT / "9.normalcell_annotation"
OUT.mkdir(parents=True, exist_ok=True)


# ── Fail-fast model guard ─────────────────────────────────────────────────
# SCimilarity v1.1 layout expected by CellAnnotation():
#   model_v1.1/encoder.ckpt
#   model_v1.1/gene_order.tsv
#   model_v1.1/annotation/labelled_kNN.bin
#   model_v1.1/annotation/reference_labels.tsv
# If any are missing the download is incomplete — bail with a clear error
# rather than crashing inside CellAnnotation() with an opaque traceback.
required = [
    MODEL / "encoder.ckpt",
    MODEL / "gene_order.tsv",
    MODEL / "annotation" / "labelled_kNN.bin",
    MODEL / "annotation" / "reference_labels.tsv",
]
for p in required:
    if not p.exists() or p.stat().st_size == 0:
        sys.exit(f"[FATAL] {p} missing/empty — model download incomplete")


# ── Load + layer swap ─────────────────────────────────────────────────────
# srt2anndata() in functions.R writes X = log-norm data, layers['counts'] = raw.
# lognorm_counts() needs raw integers in .X.
print(f"Loading {IN_H5AD}", flush=True)
adata = sc.read_h5ad(IN_H5AD)
assert "counts" in adata.layers, "layers['counts'] missing — re-export with srt2anndata"
adata.X = adata.layers["counts"].copy()
print(f"  {adata.n_obs} cells × {adata.n_vars} genes", flush=True)


# ── Preprocess & annotate ─────────────────────────────────────────────────
ca = CellAnnotation(model_path=str(MODEL))
adata = align_dataset(adata, ca.gene_order)
lognorm_counts(adata)
ca.annotate_dataset(adata)


# ── Low-confidence relabel ────────────────────────────────────────────────
DIST_THRESHOLD = 0.05
adata.obs["celltype_hint_raw"] = adata.obs["celltype_hint"]
adata.obs.loc[adata.obs["min_dist"] > DIST_THRESHOLD, "celltype_hint"] = "Unknown"
print(f"\ncelltype_hint counts (min_dist > {DIST_THRESHOLD} → Unknown):", flush=True)
print(adata.obs["celltype_hint"].value_counts(), flush=True)


# ── Save annotated h5ad + per-cell CSV ────────────────────────────────────
out_h5ad = OUT / "annotated_normal_cells.h5ad"
adata.write_h5ad(out_h5ad, compression="gzip")
print(f"\nWrote {out_h5ad}", flush=True)

out_csv = OUT / "celltype_hint_per_cell.csv"
adata.obs[["celltype_hint", "celltype_hint_raw", "min_dist"]].to_csv(out_csv)
print(f"Wrote {out_csv}", flush=True)


# ── UMAP plot ─────────────────────────────────────────────────────────────
# Re-use the batch-corrected pearson UMAP if present; otherwise any X_*umap.
umap_key = next(
    (k for k in adata.obsm if "batch" in k.lower() and "umap" in k.lower()),
    None,
)
if umap_key is None:
    umap_key = next((k for k in adata.obsm if "umap" in k.lower()), None)

if umap_key is None:
    print("[warn] no UMAP embedding found in adata.obsm — skipping UMAP plot",
          flush=True)
else:
    print(f"Plotting UMAP from obsm['{umap_key}']", flush=True)
    adata.obsm["X_umap"] = adata.obsm[umap_key]
    sc.pl.umap(adata, color="celltype_hint", show=False, frameon=False,
               save="_celltype_hint.png")
    src = Path("figures") / "umap_celltype_hint.png"
    if src.exists():
        os.replace(src, OUT / "umap_celltype_hint.png")
        print(f"Wrote {OUT / 'umap_celltype_hint.png'}", flush=True)

print("\nDone.", flush=True)
