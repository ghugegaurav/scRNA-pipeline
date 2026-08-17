"""
Stage 10 - Cell type annotation.

Corresponds to slides "What is Cell Annotation?" / "Different Annotation
Methods":
  1. Manual annotation  -> score canonical marker genes per cluster and
     assign the label whose markers score highest (slide example: high CD3E
     -> T cells, high CD19 -> B cells).
  3. Machine-learning based -> CellTypist, a logistic-regression classifier
     pretrained on large labeled reference atlases. The slide's own example
     loads `celltypist.models.Model.load("Cells_Adult_Breast.pkl")` for a
     breast dataset; here we default to "Immune_All_Low.pkl" because pbmc3k
     is immune/PBMC data -- swap the model filename in config.yaml to match
     whatever tissue your own data comes from (see celltypist.org for the
     full model list).

CellTypist needs real human gene symbols, so it only runs when
`annotation.run_celltypist: true` in config.yaml AND the dataset has gene
symbol var_names (pbmc3k does; the GitHub sample's probe IDs do not).
"""
import pandas as pd
import scanpy as sc

from utils import load_config, step, REPO_ROOT

# Slide example: "High CD3E -> T cells, High CD19 -> B cells".
# Extended with the standard pbmc3k tutorial marker set so every one of the
# dataset's major PBMC populations gets a label.
MARKER_GENES = {
    "T cells": ["CD3D", "CD3E", "CD3G", "IL7R"],
    "CD14+ Monocytes": ["CD14", "LYZ"],
    "B cells": ["CD19", "MS4A1", "CD79A"],
    "NK cells": ["GNLY", "NKG7"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
    "Dendritic cells": ["FCER1A", "CST3"],
    "Megakaryocytes": ["PPBP"],
}


def manual_annotation(adata, tables_dir):
    step("Manual marker-based annotation")
    available = {
        cell_type: [g for g in genes if g in adata.var_names]
        for cell_type, genes in MARKER_GENES.items()
    }
    available = {k: v for k, v in available.items() if v}
    if not available:
        print("  SKIPPED: none of the canonical PBMC marker genes "
              f"({sorted({g for gs in MARKER_GENES.values() for g in gs})}) "
              "are present in this dataset's var_names, so manual marker "
              "scoring does not apply here (expected for the GitHub sample's "
              "Affymetrix probe IDs).")
        return adata

    for cell_type, genes in available.items():
        score_name = f"score_{cell_type.replace(' ', '_')}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name)

    score_cols = [f"score_{k.replace(' ', '_')}" for k in available]
    cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    best_type = cluster_scores.idxmax(axis=1).str.replace("score_", "").str.replace("_", " ")
    mapping = best_type.to_dict()
    adata.obs["manual_annotation"] = adata.obs["leiden"].map(mapping).astype("category")

    print(cluster_scores.round(3))
    print("\nCluster -> manual annotation:")
    print(best_type)
    cluster_scores.to_csv(tables_dir / "manual_annotation_scores.csv")
    return adata


def celltypist_annotation(adata, cfg):
    step("Automated annotation (CellTypist)")
    if not cfg["annotation"]["run_celltypist"]:
        print("  SKIPPED: set annotation.run_celltypist: true in config.yaml to enable.")
        return adata
    if not adata.uns.get("has_gene_symbols", False):
        print("  SKIPPED: CellTypist matches genes by human gene symbol; this "
              "dataset's var_names are not gene symbols (expected for the "
              "GitHub sample).")
        return adata

    import celltypist
    from celltypist import models

    model_name = cfg["annotation"]["celltypist_model"]
    models.download_models(model=model_name)  # no-op if already cached locally
    model = models.Model.load(model_name)

    # CellTypist expects log1p-normalized-to-10k data, matching what stage 05
    # already produced before scaling overwrote adata.X -- use the saved layer.
    adata_for_ct = adata.copy()
    if "lognorm" in adata_for_ct.layers:
        adata_for_ct.X = adata_for_ct.layers["lognorm"]

    predictions = celltypist.annotate(adata_for_ct, model=model, majority_voting=True)
    result = predictions.to_adata()
    adata.obs["celltypist_label"] = result.obs["predicted_labels"]
    adata.obs["celltypist_majority_voting"] = result.obs["majority_voting"]
    print(adata.obs[["leiden", "celltypist_majority_voting"]].value_counts().head(20))
    return adata


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    tables_dir = REPO_ROOT / cfg["paths"]["tables_dir"]
    tables_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(processed_dir / "adata_07_batch.h5ad")

    adata = manual_annotation(adata, tables_dir)
    adata = celltypist_annotation(adata, cfg)

    out_path = processed_dir / "adata_final.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved final annotated object -> {out_path}")
    print("Pipeline complete. Open notebooks/scRNA_teaching_walkthrough.ipynb to explore the results.")


if __name__ == "__main__":
    main()
