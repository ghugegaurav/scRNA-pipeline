#!/usr/bin/env bash
# One-time setup script -- run ONCE inside your WSL Ubuntu terminal.
#
#   bash setup/install_env.sh
#
# What it does:
#   1. Installs Miniforge (a lightweight conda distribution defaulting to the
#      free conda-forge channel) into $HOME/miniforge3, if not already present.
#   2. Creates the "scrna-teaching" conda environment from environment.yml.
#   3. Registers it as a Jupyter kernel so it shows up in JupyterLab.
#
# Safe to re-run: every step is skipped if already done.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."   # repo root (scRNA_Teaching_Pipeline/)

MINIFORGE_DIR="$HOME/miniforge3"
# Always-current "latest release" redirect maintained by the conda-forge project:
# https://github.com/conda-forge/miniforge
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

echo "== Step 1/3: conda =="
if command -v conda &>/dev/null; then
    echo "conda already installed: $(command -v conda)"
elif [ -x "$MINIFORGE_DIR/bin/conda" ]; then
    echo "Miniforge already installed at $MINIFORGE_DIR"
    source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
else
    echo "Downloading Miniforge from $MINIFORGE_URL ..."
    curl -fsSL "$MINIFORGE_URL" -o /tmp/miniforge_installer.sh
    bash /tmp/miniforge_installer.sh -b -p "$MINIFORGE_DIR"
    rm -f /tmp/miniforge_installer.sh
    "$MINIFORGE_DIR/bin/conda" init bash
    source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
    echo "Miniforge installed. (A new terminal, or 'source ~/.bashrc', will be needed"
    echo "for the 'conda' command to be available in future sessions.)"
fi

# Make sure `conda` resolves in this script even on a fresh install.
if ! command -v conda &>/dev/null; then
    source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
fi

echo ""
echo "== Step 2/3: scrna-teaching environment =="
if conda env list | grep -qE '^\s*scrna-teaching\s'; then
    echo "Environment 'scrna-teaching' already exists -- updating it from environment.yml"
    conda env update -n scrna-teaching -f environment.yml --prune
else
    echo "Creating environment 'scrna-teaching' from environment.yml (this can take several minutes)..."
    conda env create -f environment.yml
fi

echo ""
echo "== Step 3/3: Jupyter kernel =="
conda run -n scrna-teaching python -m ipykernel install --user --name scrna-teaching --display-name "Python (scrna-teaching)"

echo ""
echo "Setup complete. Next steps:"
echo "  conda activate scrna-teaching"
echo "  bash run_pipeline.sh"
echo "or open notebooks/scRNA_teaching_walkthrough.ipynb in JupyterLab:"
echo "  jupyter lab"
