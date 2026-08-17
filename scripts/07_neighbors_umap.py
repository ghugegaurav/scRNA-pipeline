"""
Stage 07 - Neighbor graph and UMAP.

Corresponds to slides:
  "Why Neighbors?" -> using the first N PCs, find each cell's k nearest
      neighbors and connect them into a graph; clustering and UMAP both
      operate on this graph, not on the raw expression matrix.
  "What is UMAP?"   -> Uniform Manifold Approximation and Projection lays the
      neighbor graph out in 2D while preserving local (and some global)
      structure, so biologists can visually spot cell populations, batch
      effects, and heterogeneity.
"""
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_04_pca.h5ad")

    n_pcs = adata.obsm["X_pca"].shape[1]
    n_neighbors = min(cfg["dimred"]["n_neighbors"], adata.n_obs - 1)

    step(f"Neighbor graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    print("  k-nearest-neighbor graph built on PCA space.")

    step("UMAP")
    sc.tl.umap(adata)
    print("  UMAP embedding stored in adata.obsm['X_umap'].")

    out_path = processed_dir / "adata_05_umap.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/08_clustering_de.py")


if __name__ == "__main__":
    main()
