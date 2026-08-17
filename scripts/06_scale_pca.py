"""
Stage 06 - Scaling and PCA.

Corresponds to slides:
  "What is Scaling?" -> subtract the mean, divide by the standard deviation
      per gene, so no single high-magnitude gene dominates PCA.
  "What is PCA?"      -> reduces ~2000 gene dimensions down to a handful of
      principal components that capture the major axes of variation, making
      downstream neighbor-graph/clustering steps both faster and less noisy.
"""
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_03_hvg.h5ad")

    step("Scaling (zero mean, unit variance per gene)")
    adata.layers["lognorm"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)
    print("  scaled (clipped at max_value=10 to limit the influence of extreme outliers)")

    step("PCA")
    n_pcs = min(cfg["dimred"]["n_pcs"], adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    var_explained = adata.uns["pca"]["variance_ratio"].sum()
    print(f"  computed {n_pcs} PCs, explaining {var_explained:.1%} of total variance")

    out_path = processed_dir / "adata_04_pca.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/07_neighbors_umap.py")


if __name__ == "__main__":
    main()
