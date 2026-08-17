"""
Stage 02 - Build a standardized AnnData object.

Corresponds to slides: "4) Different File Formats in Single Cell" and
"5) What is AnnData (adata)?". Whatever the source format (10x .h5/.mtx or
a plain .csv), this script's job is to end up with one consistent object:

    adata.X      expression matrix (cells x genes)
    adata.obs    per-cell metadata
    adata.var    per-gene metadata
    adata.uns    free-form info, incl. which dataset/track this run is using

saved to data/processed/adata_00_raw.h5ad, exactly like the slide's
`adata.write_h5ad("BC_data_combined.h5ad")` example.
"""
import anndata as ad
import pandas as pd
import scanpy as sc

from utils import load_config, step, has_gene_symbol_annotation, looks_like_raw_counts, REPO_ROOT


def build_pbmc3k(raw_dir):
    adata = sc.read_h5ad(raw_dir / "pbmc3k_raw.h5ad")
    adata.var_names_make_unique()
    adata.obs["sample"] = "pbmc3k"
    adata.uns["data_source"] = "pbmc3k (10x Genomics, real raw UMI counts)"
    return adata


def build_github_sample(raw_dir):
    csv_path = raw_dir / "RNA-seq.csv"
    df = pd.read_csv(csv_path, index_col=0)  # rows = cells, columns = probes
    adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    adata.var_names_make_unique()
    adata.obs["sample"] = "github_sample"

    labels_path = raw_dir / "labels.csv"
    if labels_path.exists():
        labels = pd.read_csv(labels_path, index_col=0)
        labels = labels.reindex(adata.obs_names)
        adata.obs["reference_label"] = labels["labels"].astype("category")
        print(f"Attached ground-truth 'reference_label' for {labels['labels'].notna().sum()} / {adata.n_obs} cells "
              f"(from the original assignment's data/labels.csv).")

    adata.uns["data_source"] = (
        "GiatrasKon/scRNAseq-Analysis-Pipeline sample (137 cells x 54,675 Affymetrix "
        "microarray probes, already log2-normalized -- NOT raw UMI counts)"
    )
    return adata


def main():
    cfg = load_config()
    raw_dir = REPO_ROOT / cfg["paths"]["raw_dir"]
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    processed_dir.mkdir(parents=True, exist_ok=True)

    step("Building AnnData object")
    dataset = cfg["dataset"]
    if dataset == "pbmc3k":
        adata = build_pbmc3k(raw_dir)
    elif dataset == "github_sample":
        adata = build_github_sample(raw_dir)
    else:
        raise ValueError(f"Unknown dataset '{dataset}' in config.yaml")

    # Record data characteristics that later stages use to decide which QC
    # steps are scientifically applicable (see utils.looks_like_raw_counts /
    # utils.has_gene_symbol_annotation).
    adata.uns["is_raw_counts"] = looks_like_raw_counts(adata.X)
    adata.uns["has_gene_symbols"] = has_gene_symbol_annotation(adata.var_names)

    print(f"\nadata: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  data_source        : {adata.uns['data_source']}")
    print(f"  looks like raw counts?      {adata.uns['is_raw_counts']}")
    print(f"  has gene-symbol var_names?  {adata.uns['has_gene_symbols']}")

    out_path = processed_dir / "adata_00_raw.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/03_qc_metrics.py")


if __name__ == "__main__":
    main()
