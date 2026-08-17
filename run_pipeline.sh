#!/usr/bin/env bash
# Master runner -- executes every pipeline stage in order.
#
# Run this from inside WSL, from the scRNA_Teaching_Pipeline/ directory:
#   conda activate scrna-teaching
#   bash run_pipeline.sh
#
# Each stage reads the previous stage's saved .h5ad from data/processed/ and
# writes its own, so you can also re-run a single stage on its own, e.g.:
#   python scripts/08_clustering_de.py
set -euo pipefail

if [[ "$(uname -r)" != *microsoft* && "$(uname -r)" != *Microsoft* ]]; then
    echo "WARNING: this does not look like a WSL environment (uname -r: $(uname -r))." >&2
    echo "The pipeline will still run on plain Linux/macOS, but the guide book assumes WSL." >&2
fi

cd "$(dirname "${BASH_SOURCE[0]}")"

if ! python -c "import scanpy" &>/dev/null; then
    echo "ERROR: the 'scrna-teaching' conda environment is not active (or not installed)." >&2
    echo "Run: conda activate scrna-teaching" >&2
    echo "If it doesn't exist yet: bash setup/install_env.sh" >&2
    exit 1
fi

cd scripts
STAGES=(
    01_fetch_data.py
    02_build_anndata.py
    03_qc_metrics.py
    04_doublet_detection.py
    05_normalize_hvg.py
    06_scale_pca.py
    07_neighbors_umap.py
    08_clustering_de.py
    09_batch_correction.py
    10_annotation.py
    11_generate_figures.py
)

for stage in "${STAGES[@]}"; do
    echo ""
    echo "############################################################"
    echo "# Running $stage"
    echo "############################################################"
    python "$stage"
done

echo ""
echo "Pipeline finished. Explore results in ../results/ or open"
echo "../notebooks/scRNA_teaching_walkthrough.ipynb"
