"""
Stage 04 - Doublet detection.

Corresponds to slide "What is a Doublet?": two cells captured in one droplet
produce an artificial hybrid transcriptome that can masquerade as a fake
cluster. Scrublet simulates artificial doublets by combining pairs of real
cells and scores every real cell by how similar it looks to a simulated one.

Scrublet's simulation only makes sense on raw (integer) UMI counts -- if the
loaded dataset is already normalized (no raw counts available), this stage
is skipped with an explanation rather than run on the wrong kind of values.
"""
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    adata = sc.read_h5ad(processed_dir / "adata_01_qc.h5ad")

    step("Doublet detection (Scrublet)")
    if adata.uns.get("is_raw_counts", False):
        sc.pp.scrublet(adata, expected_doublet_rate=cfg["doublets"]["expected_doublet_rate"])
        n_doublets = int(adata.obs["predicted_doublet"].sum())
        print(f"  predicted doublets: {n_doublets} / {adata.n_obs} cells "
              f"(mean doublet_score={adata.obs['doublet_score'].mean():.3f})")
        n_before = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        print(f"  removed predicted doublets: {n_before} -> {adata.n_obs} cells")
    else:
        adata.obs["doublet_score"] = float("nan")
        adata.obs["predicted_doublet"] = False
        print("  SKIPPED: Scrublet simulates doublets by adding together raw UMI "
              "counts from real cell pairs, which is only meaningful on raw counts. "
              "This dataset is already normalized, so this step is skipped "
              "(expected for the GitHub sample).")

    out_path = processed_dir / "adata_02_doublets.h5ad"
    adata.write_h5ad(out_path)
    print(f"\nSaved -> {out_path}")
    print("Next: python scripts/05_normalize_hvg.py")


if __name__ == "__main__":
    main()
