"""
Stage 05 - Normalization, log transform, highly variable genes.

Corresponds to slides:
  "What is sc.pp and sc.tl in Scanpy?" / "Normalization"
      -> sc.pp.normalize_total(adata, target_sum=1e4): each cell is scaled so
         sequencing-depth differences between cells stop looking biological.
  "Log Transformation"
      -> sc.pp.log1p(adata): compresses the long right tail of expression
         values so a few very high genes don't dominate downstream PCA.
  "Highly Variable Genes (HVG)"
      -> sc.pp.highly_variable_genes(..., n_top_genes=2000, flavor='seurat'):
         keep genes that vary across cells instead of flat housekeeping genes.

If the dataset arrived already normalized (no raw counts), normalize_total /
log1p are skipped -- applying them again would double-transform the values.
"""
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_02_doublets.h5ad")

    adata.layers["counts"] = adata.X.copy()  # slide: adata.layers["counts"]

    step("Normalization + log transform")
    if adata.uns.get("is_raw_counts", False):
        target_sum = cfg["normalize"]["target_sum"]
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        print(f"  normalize_total(target_sum={target_sum:g}) + log1p applied.")
    else:
        print("  SKIPPED normalize_total/log1p: this dataset's values are already "
              "normalized (e.g. log2 RMA-normalized microarray intensities for the "
              "GitHub sample), so re-normalizing would distort them.")

    step("Highly variable genes (HVG)")
    n_top_genes = cfg["hvg"]["n_top_genes"]
    n_top_genes = min(n_top_genes, adata.n_vars - 1)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, subset=True, flavor=cfg["hvg"]["flavor"]
    )
    print(f"  kept {adata.n_vars} highly variable genes (flavor='{cfg['hvg']['flavor']}')")

    out_path = processed_dir / "adata_03_hvg.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/06_scale_pca.py")


if __name__ == "__main__":
    main()
