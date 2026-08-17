"""
Stage 08 - Clustering and differential expression.

Corresponds to slides:
  "Graph Theory Logic" (Leiden) -> starts with every cell as its own cluster,
      then merges cells to maximize modularity (are connections stronger
      within a cluster than between clusters?).
  "What is Resolution?" -> low resolution = few big clusters, high resolution
      = many small clusters. We run the configured resolution plus a low/high
      pair so students can see the effect directly.
  "What is rank_genes_groups?" -> for each cluster, Wilcoxon rank-sum test of
      that cluster's cells vs. all other cells -> marker genes per cluster.

When running on the GitHub sample (which ships a `reference_label` column
from the original assignment's data/labels.csv), this stage also compares
Leiden's clusters against that ground truth, and reproduces the original
repo's own GMM/DBSCAN-on-PCA approach for a side-by-side comparison between
"classical ML clustering" (the assignment's own method) and the
graph-based Leiden workflow taught in this course.
"""
import numpy as np
import pandas as pd
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def compare_to_classical_ml(adata, tables_dir):
    """Reproduce the original repo's GMM + DBSCAN approach on the PCA
    embedding, and compare all clusterings to the reference_label ground
    truth using Adjusted Rand Index / Normalized Mutual Information."""
    from sklearn.cluster import DBSCAN
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
    from sklearn.mixture import GaussianMixture

    step("Comparing Leiden vs. the original repo's classical ML clustering (GMM / DBSCAN)")
    X_pca = adata.obsm["X_pca"]
    rows = []

    # GMM with BIC-selected number of components (same idea as the original
    # repo's codebase.py, using a smaller search range for speed here).
    best_bic, best_gmm = np.inf, None
    for k in range(2, 11):
        gmm = GaussianMixture(n_components=k, covariance_type="full", random_state=42).fit(X_pca)
        bic = gmm.bic(X_pca)
        if bic < best_bic:
            best_bic, best_gmm = bic, gmm
    gmm_labels = best_gmm.predict(X_pca)
    adata.obs["gmm_cluster"] = pd.Categorical(gmm_labels.astype(str))

    dbscan_labels = DBSCAN(eps=1.0, min_samples=5).fit_predict(X_pca)
    adata.obs["dbscan_cluster"] = pd.Categorical(dbscan_labels.astype(str))

    for name, labels in [("leiden", adata.obs["leiden"]), ("gmm", gmm_labels), ("dbscan", dbscan_labels)]:
        n_clusters = len(set(labels)) - (1 if -1 in set(labels) else 0)
        row = {"method": name, "n_clusters": n_clusters}
        if n_clusters > 1:
            row["silhouette"] = silhouette_score(X_pca, labels)
        if "reference_label" in adata.obs and adata.obs["reference_label"].notna().any():
            ref = adata.obs["reference_label"]
            row["ARI_vs_reference"] = adjusted_rand_score(ref, labels)
            row["NMI_vs_reference"] = normalized_mutual_info_score(ref, labels)
        rows.append(row)

    comparison = pd.DataFrame(rows).set_index("method")
    print(comparison)
    comparison.to_csv(tables_dir / "clustering_comparison_vs_classical_ml.csv")
    print(f"Saved -> {tables_dir / 'clustering_comparison_vs_classical_ml.csv'}")


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    tables_dir = REPO_ROOT / cfg["paths"]["tables_dir"]
    tables_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(processed_dir / "adata_05_umap.h5ad")

    resolution = cfg["clustering"]["resolution"]
    step(f"Leiden clustering (resolution={resolution})")
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", flavor="igraph", n_iterations=2)
    print(f"  {adata.obs['leiden'].nunique()} clusters at resolution={resolution}")

    step("Resolution sweep (slide: 'What is Resolution?')")
    for tag, res in [("low", 0.2), ("high", 2.0)]:
        key = f"leiden_{tag}"
        sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2)
        print(f"  resolution={res:<4} ({tag:<4}) -> {adata.obs[key].nunique()} clusters")

    step("Differential expression: rank_genes_groups (Wilcoxon rank-sum test)")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(tables_dir / "cluster_markers.csv", index=False)
    print(f"  top marker per cluster:")
    print(markers.sort_values("scores", ascending=False).groupby("group").head(1)[["group", "names", "scores"]])
    print(f"Saved full marker table -> {tables_dir / 'cluster_markers.csv'}")

    if "reference_label" in adata.obs.columns:
        compare_to_classical_ml(adata, tables_dir)

    out_path = processed_dir / "adata_06_clustered.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/09_batch_correction.py")


if __name__ == "__main__":
    main()
