"""
Stage 11 - Generate teaching figures.

Saves a handful of PNGs to results/figures/ so students have something to
look at without needing to run Jupyter: QC violin plots, UMAP colored by
cluster/annotation, and (dense vs sparse) a marker-gene dotplot -- tying
back to slides "6) Dense vs Sparse Matrix" (memory) and the annotation
slides (biological interpretation).
"""
import matplotlib

matplotlib.use("Agg")  # headless (no display needed inside WSL)
import matplotlib.pyplot as plt
import scanpy as sc

from utils import load_config, step, REPO_ROOT


def savefig(fig_dir, name):
    path = fig_dir / name
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved -> {path}")


def main():
    cfg = load_config()
    processed_dir = REPO_ROOT / cfg["paths"]["processed_dir"]
    fig_dir = REPO_ROOT / cfg["paths"]["figures_dir"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(processed_dir / "adata_final.h5ad")

    sc.settings.figdir = fig_dir

    step("QC violin plots")
    qc_vars = [v for v in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if v in adata.obs]
    if qc_vars:
        sc.pl.violin(adata, qc_vars, jitter=0.4, multi_panel=True, show=False)
        savefig(fig_dir, "01_qc_violin.png")

    step("UMAP colored by Leiden cluster")
    sc.pl.umap(adata, color="leiden", show=False)
    savefig(fig_dir, "02_umap_leiden.png")

    if "manual_annotation" in adata.obs:
        step("UMAP colored by manual annotation")
        sc.pl.umap(adata, color="manual_annotation", show=False)
        savefig(fig_dir, "03_umap_manual_annotation.png")

    if "celltypist_majority_voting" in adata.obs:
        step("UMAP colored by CellTypist annotation")
        sc.pl.umap(adata, color="celltypist_majority_voting", show=False)
        savefig(fig_dir, "04_umap_celltypist.png")

    if "batch" in adata.obs:
        step("UMAP colored by (synthetic) batch")
        sc.pl.umap(adata, color="batch", show=False)
        savefig(fig_dir, "05_umap_batch.png")

    step("Top marker genes per cluster (dotplot)")
    if "rank_genes_groups" in adata.uns:
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, show=False)
        savefig(fig_dir, "06_marker_dotplot.png")

    print(f"\nAll figures saved under {fig_dir}")


if __name__ == "__main__":
    main()
