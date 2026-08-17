"""
Stage 03 - Quality control.

Corresponds to slides:
  "Why Quality Control is Needed?"        -> dead/broken cells, empty droplets,
                                              doublets, ambient RNA contamination.
  "1) Mitochondrial Genes (MT-)"          -> dying cells leak cytoplasmic RNA,
                                              so mitochondrial RNA % spikes.
  "2) Ribosomal Genes (RPS, RPL)"         -> housekeeping genes; dominance can
                                              signal technical bias.
  "3) Hemoglobin Genes (HB)"              -> flags red-blood-cell contamination.
  "What is MAD (Median Absolute Deviation)?" -> robust outlier detection.

If the loaded dataset has no human gene-symbol var_names (e.g. the GitHub
sample's Affymetrix probe IDs), the mito/ribo/hb steps are skipped with an
explanation instead of silently producing zeros.
"""
import scanpy as sc

from utils import load_config, step, is_outlier, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_00_raw.h5ad")

    step("Step 1: Flagging mitochondrial / ribosomal / hemoglobin genes")
    if adata.uns.get("has_gene_symbols", False):
        adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
        adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB[^(P)]")
        print(f"  mitochondrial genes found: {adata.var['mt'].sum()}")
        print(f"  ribosomal genes found:     {adata.var['ribo'].sum()}")
        print(f"  hemoglobin genes found:    {adata.var['hb'].sum()}")
    else:
        adata.var["mt"] = False
        adata.var["ribo"] = False
        adata.var["hb"] = False
        print("  SKIPPED: var_names are not human gene symbols "
              f"(e.g. '{adata.var_names[0]}'), so MT-/RPS-/RPL-/HB prefix "
              "matching does not apply to this dataset. This is expected for "
              "the GitHub sample (Affymetrix probe IDs).")

    step("Step 2: Calculating QC metrics (sc.pp.calculate_qc_metrics)")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True
    )
    print(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

    step("Step 3: MAD-based outlier detection (nmads=%d)" % cfg["qc"]["mad_nmads"])
    nmads = cfg["qc"]["mad_nmads"]
    adata.obs["outlier"] = (
        is_outlier(adata.obs["log1p_total_counts"], nmads)
        | is_outlier(adata.obs["log1p_n_genes_by_counts"], nmads)
        | is_outlier(adata.obs["pct_counts_in_top_20_genes"], nmads)
    )
    print(f"  flagged as statistical outliers: {adata.obs['outlier'].sum()} / {adata.n_obs} cells")

    if adata.uns.get("has_gene_symbols", False):
        max_pct_mt = cfg["qc"]["max_pct_mt"]
        adata.obs["mt_outlier"] = is_outlier(adata.obs["pct_counts_mt"], nmads) | (
            adata.obs["pct_counts_mt"] > max_pct_mt
        )
        print(f"  flagged for high mitochondrial %% (> {max_pct_mt}%% or MAD outlier): "
              f"{adata.obs['mt_outlier'].sum()} / {adata.n_obs} cells")
    else:
        adata.obs["mt_outlier"] = False

    step("Step 4: Filtering low-quality cells and rarely-detected genes")
    n_before = adata.n_obs
    keep = ~(adata.obs["outlier"] | adata.obs["mt_outlier"])
    adata = adata[keep].copy()
    print(f"  cells: {n_before} -> {adata.n_obs} (removed {n_before - adata.n_obs})")

    sc.pp.filter_cells(adata, min_genes=cfg["qc"]["min_genes_per_cell"])
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=cfg["qc"]["min_cells_per_gene"])
    print(f"  genes: {n_genes_before} -> {adata.n_vars} (removed genes detected in "
          f"< {cfg['qc']['min_cells_per_gene']} cells)")
    print(f"  cells after min_genes={cfg['qc']['min_genes_per_cell']} filter: {adata.n_obs}")

    out_path = processed_dir / "adata_01_qc.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/04_doublet_detection.py")


if __name__ == "__main__":
    main()
