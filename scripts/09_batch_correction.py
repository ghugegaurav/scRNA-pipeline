"""
Stage 09 - Batch correction demo.

Corresponds to slide "What is Batch Effect?": technical variation (different
reagents/day/depth) that is unrelated to biology, which can make a UMAP
separate cells by batch instead of by cell type.
  BBKNN  -> forces each cell's neighbors to be drawn evenly from every batch.
  Harmony -> keeps PCA-space, but iteratively shifts each batch's embedding
             to align with the others while preserving biological structure.

Neither bundled dataset (pbmc3k, or the GitHub sample) has a real
experimental batch variable, so this stage creates a SYNTHETIC 2-way split
purely so students can see the mechanics of BBKNN/Harmony run end-to-end.
This is clearly labeled -- it is a mechanics demo, not a real batch-effect
correction, and is skipped by setting `batch_correction.method: "none"` in
config.yaml.
"""
import numpy as np
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_06_clustered.h5ad")

    method = cfg["batch_correction"]["method"]
    step(f"Batch correction demo (method={method})")

    if method == "none":
        print("  SKIPPED (config.yaml: batch_correction.method = 'none')")
        adata.write_h5ad(processed_dir / "adata_07_batch.h5ad")
        return

    if "batch" not in adata.obs.columns:
        rng = np.random.default_rng(42)
        adata.obs["batch"] = rng.choice(["batch_1", "batch_2"], size=adata.n_obs)
        print("  NOTE: this dataset has no real experimental batch variable, so a "
              "SYNTHETIC random 2-way split was created purely to demonstrate how "
              f"{method} is invoked. Do not interpret the 'corrected' UMAP as "
              "removing a real technical effect -- there is nothing real to remove "
              "here. On your own multi-sample data, adata.obs['batch'] should come "
              "from real sample/lane/run metadata.")

    if method == "harmony":
        sc.external.pp.harmony_integrate(adata, key="batch")
        rep = "X_pca_harmony"
        print(f"  Harmony-adjusted PCA embedding stored in adata.obsm['{rep}']")
    elif method == "bbknn":
        import scanpy.external as sce

        sce.pp.bbknn(adata, batch_key="batch")
        rep = None
        print("  BBKNN-balanced neighbor graph replaces the standard neighbor graph.")
    else:
        raise ValueError(f"Unknown batch_correction.method '{method}' in config.yaml")

    step("Recomputing neighbors/UMAP on the batch-corrected representation")
    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_pre_batch_correction"] = adata.obsm["X_umap"].copy()
    if rep:
        sc.pp.neighbors(adata, use_rep=rep)
    sc.tl.umap(adata)
    print("  batch-corrected UMAP stored in adata.obsm['X_umap'] "
          "(previous UMAP kept in adata.obsm['X_umap_pre_batch_correction']).")

    out_path = processed_dir / "adata_07_batch.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/10_annotation.py")


if __name__ == "__main__":
    main()
