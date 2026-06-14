"""SCimilarity cell-type annotation + cell search for the normal-cell integration.

Reads ~/VisHD/8.1.normal_cell_integration/integrated_normal_cells2.h5ad and:
  1. Runs Genentech's SCimilarity foundation model (CellAnnotation) against the
     prepackaged 23 M-cell reference, relabels low-confidence calls as "Unknown".
  2. Builds a SCimilarity embedding UMAP + Leiden clustering ("new embedding").
  3. Signature-based cell search (cell_search_tutorial_2): for each meta-program
     signature (same gene sets FeaturePlotted in 9.2.scimilarity_check.R, plus
     SVEC), scores cells, marks the top 5%, builds a centroid and finds which
     reference cell types it matches.
  4. Cluster-based cell type search (cell_search_tutorial_3): predicts each
     SCimilarity-Leiden cluster's cell type from the reference.

Writes the annotated AnnData + per-cell CSV + UMAP PNG, the signature search
results, and the cluster search results into ~/VisHD/9.normalcell_annotation.

Tutorials:
  https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
  https://genentech.github.io/scimilarity/notebooks/cell_search_tutorial_2.html
  https://genentech.github.io/scimilarity/notebooks/cell_search_tutorial_3.html
"""

from pathlib import Path
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scanpy as sc
from scimilarity import CellAnnotation, CellQuery
from scimilarity.utils import align_dataset, lognorm_counts


# ── Paths ──────────────────────────────────────────────────────────────────
# Use $MYSCRATCH as the project root — works both natively on Pawsey and
# under the singularity bind-mount that maps $MYSCRATCH to $HOME.
ROOT      = Path(os.environ.get("MYSCRATCH", str(Path.home()))) / "VisHD"
IN_H5AD   = ROOT / "8.1.normal_cell_integration/integrated_normal_cells2.h5ad"
MODEL     = Path("/scratch/pawsey1172/sweng/VisHD/SCimilarity/model_v1.1")
META_XLSX = ROOT / "public_signature/meta_programs_2025-01-29.xlsx"
OUT       = ROOT / "9.normalcell_annotation"
SIG_DIR   = OUT / "cellsearch_signature"
CLU_DIR   = OUT / "cellsearch_cluster"
for d in (OUT, SIG_DIR, SIG_DIR / "figures", CLU_DIR, CLU_DIR / "figures"):
    d.mkdir(parents=True, exist_ok=True)


# ── Fail-fast model guard ─────────────────────────────────────────────────
# SCimilarity v1.1 layout expected by CellAnnotation() + CellQuery():
#   model_v1.1/encoder.ckpt
#   model_v1.1/gene_order.tsv
#   model_v1.1/annotation/labelled_kNN.bin
#   model_v1.1/annotation/reference_labels.tsv
#   model_v1.1/cellsearch/full_kNN.bin
#   model_v1.1/cellsearch/cell_metadata/   (TileDB)
#   model_v1.1/cellsearch/cell_embedding/  (TileDB)
# If any are missing the download is incomplete — bail with a clear error
# rather than crashing inside CellAnnotation()/CellQuery() with an opaque traceback.
required_files = [
    MODEL / "encoder.ckpt",
    MODEL / "gene_order.tsv",
    MODEL / "annotation" / "labelled_kNN.bin",
    MODEL / "annotation" / "reference_labels.tsv",
    MODEL / "cellsearch" / "full_kNN.bin",
]
for p in required_files:
    if not p.exists() or p.stat().st_size == 0:
        sys.exit(f"[FATAL] {p} missing/empty — model download incomplete")
required_dirs = [
    MODEL / "cellsearch" / "cell_metadata",
    MODEL / "cellsearch" / "cell_embedding",
]
for p in required_dirs:
    if not p.is_dir():
        sys.exit(f"[FATAL] {p} missing — cell search index incomplete")


# ── Load + layer swap ─────────────────────────────────────────────────────
# srt2anndata() in functions.R writes X = log-norm data, layers['counts'] = raw.
# lognorm_counts() needs raw integers in .X.
print(f"Loading {IN_H5AD}", flush=True)
adata = sc.read_h5ad(IN_H5AD)
assert "counts" in adata.layers, "layers['counts'] missing — re-export with srt2anndata"
adata.X = adata.layers["counts"].copy()
print(f"  {adata.n_obs} cells × {adata.n_vars} genes", flush=True)


# ── Preprocess & annotate ─────────────────────────────────────────────────
# After align_dataset + lognorm_counts: X = lognorm, layers['counts'] = aligned
# raw counts (both required by the cell search below). annotate_dataset() also
# populates obsm['X_scimilarity'].
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

# Per-cell hint CSV (read downstream by 9.2.scimilarity_check.R)
out_csv = OUT / "celltype_hint_per_cell.csv"
adata.obs[["celltype_hint", "celltype_hint_raw", "min_dist"]].to_csv(out_csv)
print(f"Wrote {out_csv}", flush=True)


# ── celltype_hint UMAP (re-use the pearson batch UMAP if present) ──────────
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
    print(f"Plotting celltype_hint UMAP from obsm['{umap_key}']", flush=True)
    adata.obsm["X_umap"] = adata.obsm[umap_key]
    sc.pl.umap(adata, color="celltype_hint", show=False, frameon=False,
               save="_celltype_hint.png")
    src = Path("figures") / "umap_celltype_hint.png"
    if src.exists():
        os.replace(src, OUT / "umap_celltype_hint.png")
        print(f"Wrote {OUT / 'umap_celltype_hint.png'}", flush=True)

# Free the 8.7 G annotation kNN before loading the 21 G cell-search kNN.
del ca


# ── SCimilarity embedding UMAP + Leiden (the "new embedding") ──────────────
print("\nBuilding SCimilarity-embedding neighbours / UMAP / Leiden ...", flush=True)
sc.pp.neighbors(adata, use_rep="X_scimilarity", key_added="scim", random_state=0)
sc.tl.umap(adata, neighbors_key="scim", random_state=0)
adata.obsm["X_scimilarity_umap"] = adata.obsm["X_umap"].copy()
sc.tl.leiden(adata, neighbors_key="scim", key_added="scimilarity_leiden",
             resolution=1.0, flavor="igraph", n_iterations=2, directed=False,
             random_state=0)
n_clu = adata.obs["scimilarity_leiden"].nunique()
print(f"Leiden on X_scimilarity: {n_clu} clusters", flush=True)


# ── Build signature dict (SVEC + meta-programs, mirrors 9.2) ───────────────
SVEC = ["SEMG1", "SEMG2", "MUC6", "PGC", "CYP4F8", "CLU", "PDK4",
        "SLPI", "AKR1B1", "KRT7", "SLC26A3", "PATE1", "PAX8"]

signatures = {"SVEC": SVEC}
sheets = pd.read_excel(META_XLSX, sheet_name=None)
for sheet, df in sheets.items():
    if sheet == "Malignant":
        continue
    for col in df.columns:
        genes = [g.strip() for g in df[col].dropna().astype(str) if g.strip()]
        if genes:
            signatures[f"{sheet}__{col}"] = genes
print(f"{len(signatures)} candidate signatures (SVEC + meta-programs)", flush=True)


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def save_embedding_plot(color, path, title=None, categorical=False):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    kwargs = dict(basis="scimilarity_umap", color=color, ax=ax, show=False,
                  frameon=False, title=title)
    if categorical:
        kwargs.update(legend_fontsize=5, legend_loc="right margin")
    else:
        kwargs.update(cmap="viridis")
    sc.pl.embedding(adata, **kwargs)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ── CellQuery (loads the 21 G cellsearch kNN + reference metadata) ─────────
cq = CellQuery(model_path=str(MODEL))


# ── Signature-based cell search (cell_search_tutorial_2) ───────────────────
print("\nSignature-based cell search (top 5% query cells per signature) ...",
      flush=True)
sig_summary = []
score_cols = []
for name, genes in signatures.items():
    used = [g for g in genes if g in adata.var_names]
    if len(used) < 3:
        print(f"  [skip] {name}: only {len(used)} genes in dataset", flush=True)
        continue

    sc.tl.score_genes(adata, used, score_name="_sig", random_state=0)
    score_col = f"sigscore__{safe_name(name)}"
    adata.obs[score_col] = adata.obs["_sig"].astype(float)
    score_cols.append(score_col)

    cut = np.quantile(adata.obs["_sig"].values, 0.95)
    adata.obs["_q"] = (adata.obs["_sig"].values >= cut).astype(int)
    n_q = int(adata.obs["_q"].sum())

    _, _, _, meta, qc = cq.search_centroid_nearest(adata, "_q", k=10000, qc=True)
    vc = meta["prediction"].value_counts(normalize=True)
    coherence = float(qc.get("query_coherence", np.nan))

    meta["prediction"].value_counts().head(50).to_csv(
        SIG_DIR / f"{safe_name(name)}_celltype_counts.csv",
        header=["n_neighbors"])
    sig_summary.append({
        "signature":      name,
        "n_genes_total":  len(genes),
        "n_genes_used":   len(used),
        "n_query_cells":  n_q,
        "query_coherence": coherence,
        "top1_celltype":  vc.index[0],
        "top1_frac":      float(vc.iloc[0]),
        "top5_celltypes": "; ".join(f"{c}:{v:.3f}" for c, v in vc.head(5).items()),
    })
    save_embedding_plot(score_col, SIG_DIR / "figures" / f"{safe_name(name)}.png",
                        title=f"{name} ({len(used)} genes)")
    print(f"  {name}: query_coherence={coherence:.1f}, top={vc.index[0]} "
          f"({vc.iloc[0]:.2f})", flush=True)

adata.obs.drop(columns=["_sig", "_q"], errors="ignore", inplace=True)

sig_summary_df = pd.DataFrame(sig_summary)
sig_summary_df.to_csv(SIG_DIR / "signature_celltype_summary.csv", index=False)
if score_cols:
    adata.obs[score_cols].to_csv(SIG_DIR / "signature_scores_per_cell.csv")
print(f"Wrote signature search results to {SIG_DIR}", flush=True)


# ── Cluster-based cell type search (cell_search_tutorial_3) ────────────────
print("\nCluster-based cell type search on scimilarity_leiden ...", flush=True)
_, cluster_idx, _, _, cmeta = cq.search_cluster_centroids_nearest(
    adata, "scimilarity_leiden")

cluster_rows = []
full_counts = []
for clu, sub in cmeta.groupby("centroid", observed=True):
    vc = sub["prediction"].value_counts(normalize=True)
    n_cells = int((adata.obs["scimilarity_leiden"].astype(str) == str(clu)).sum())
    cluster_rows.append({
        "cluster":            clu,
        "n_cells_in_cluster": n_cells,
        "n_neighbors":        len(sub),
        "predicted_celltype": vc.index[0],
        "predicted_frac":     float(vc.iloc[0]),
        "top5_celltypes":     "; ".join(f"{c}:{v:.3f}" for c, v in vc.head(5).items()),
    })
    counts = sub["prediction"].value_counts().rename("n_neighbors").reset_index()
    counts.columns = ["prediction", "n_neighbors"]
    counts.insert(0, "cluster", clu)
    full_counts.append(counts)

pred_df = pd.DataFrame(cluster_rows).sort_values("cluster")
pred_df.to_csv(CLU_DIR / "cluster_celltype_predictions.csv", index=False)
pd.concat(full_counts, ignore_index=True).to_csv(
    CLU_DIR / "cluster_celltype_full_counts.csv", index=False)

mapping = dict(zip(pred_df["cluster"].astype(str), pred_df["predicted_celltype"]))
adata.obs["scimilarity_cluster_celltype"] = (
    adata.obs["scimilarity_leiden"].astype(str).map(mapping))
print(pred_df[["cluster", "n_cells_in_cluster", "predicted_celltype",
               "predicted_frac"]].to_string(index=False), flush=True)

save_embedding_plot("scimilarity_leiden",
                    CLU_DIR / "figures" / "leiden_clusters.png",
                    title="SCimilarity Leiden clusters", categorical=True)
save_embedding_plot("scimilarity_cluster_celltype",
                    CLU_DIR / "figures" / "cluster_celltype.png",
                    title="Cluster cell-type search", categorical=True)
print(f"Wrote cluster search results to {CLU_DIR}", flush=True)


# ── Save annotated h5ad (carries X_scimilarity, X_scimilarity_umap, ────────
#    scimilarity_leiden, scimilarity_cluster_celltype, sigscore__* obs) ─────
if "X_umap" in adata.obsm:          # scratch key (== X_scimilarity_umap now)
    del adata.obsm["X_umap"]
out_h5ad = OUT / "annotated_normal_cells.h5ad"
adata.write_h5ad(out_h5ad, compression="gzip")
print(f"\nWrote {out_h5ad}", flush=True)

print("\nDone.", flush=True)
