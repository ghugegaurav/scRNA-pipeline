"""
Stage 01 - Fetch data.

Downloads whichever dataset config.yaml points at:

  dataset: "pbmc3k"
      Real 10x Genomics droplet scRNA-seq data (2,700 PBMCs), raw UMI counts.
      Fetched via scanpy's own `sc.datasets.pbmc3k()`, which downloads from
      10x Genomics' official file server. This is the dataset used by the
      Scanpy project's own tutorials, so every QC/doublet/annotation step
      taught in the slides applies to it unmodified.

  dataset: "github_sample"
      The sample dataset bundled with the course GitHub repository
      (GiatrasKon/scRNAseq-Analysis-Pipeline). Downloaded directly from the
      repo's raw file URLs defined in config.yaml.

Corresponds to slide: "3) From FASTQ to Gene-Cell Matrix" -- this script is
the "Step 2: Alignment & Counting" output already sitting on a server;
we start from the finished gene-cell matrix, same as the slide does.
"""
import gzip
import shutil
from pathlib import Path

import requests
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def fetch_pbmc3k(raw_dir: Path) -> Path:
    step("Fetching pbmc3k (10x Genomics PBMC, raw UMI counts)")
    # sc.datasets.pbmc3k() downloads and caches
    # https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
    # under ./data by default; we point its cache at our raw_dir.
    sc.settings.datasetdir = raw_dir
    adata = sc.datasets.pbmc3k()
    out_path = raw_dir / "pbmc3k_raw.h5ad"
    adata.write_h5ad(out_path)
    print(f"Saved raw AnnData -> {out_path}  ({adata.n_obs} cells x {adata.n_vars} genes)")
    return out_path


def fetch_github_sample(raw_dir: Path, cfg: dict) -> Path:
    step("Fetching the sample dataset from the course GitHub repo")
    csv_gz_url = cfg["github_sample"]["raw_csv_gz_url"]
    labels_url = cfg["github_sample"]["raw_labels_url"]

    csv_gz_path = raw_dir / "RNA-seq.csv.gz"
    csv_path = raw_dir / "RNA-seq.csv"
    labels_path = raw_dir / "labels.csv"

    for url, dest in [(csv_gz_url, csv_gz_path), (labels_url, labels_path)]:
        print(f"Downloading {url}")
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        print(f"  -> {dest}  ({dest.stat().st_size:,} bytes)")

    with gzip.open(csv_gz_path, "rb") as f_in, open(csv_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed -> {csv_path}")
    return csv_path


def main():
    cfg = load_config()
    raw_dir = REPO_ROOT / cfg["paths"]["raw_dir"]
    raw_dir.mkdir(parents=True, exist_ok=True)

    dataset = cfg["dataset"]
    if dataset == "pbmc3k":
        fetch_pbmc3k(raw_dir)
    elif dataset == "github_sample":
        fetch_github_sample(raw_dir, cfg)
    else:
        raise ValueError(f"Unknown dataset '{dataset}' in config.yaml (expected 'pbmc3k' or 'github_sample')")

    print("\nDone. Next: python scripts/02_build_anndata.py")


if __name__ == "__main__":
    main()
