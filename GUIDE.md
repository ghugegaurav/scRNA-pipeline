# scRNA-seq Teaching Pipeline -- Student Guide Book

This guide walks you, step by step, through running the single-cell RNA-seq
(scRNA-seq) teaching pipeline on your Windows PC using **WSL (Windows
Subsystem for Linux)**. No prior Linux experience is assumed.

The pipeline itself teaches the exact workflow covered in the course slides
(`SC_analysis.pptx`): quality control, doublet detection, normalization,
dimensionality reduction, clustering, batch correction, and cell type
annotation, using **Scanpy** (the standard Python single-cell toolkit) --
and it can run either on a real 10x Genomics dataset or on the sample
dataset from the course's GitHub repository
([GiatrasKon/scRNAseq-Analysis-Pipeline](https://github.com/GiatrasKon/scRNAseq-Analysis-Pipeline)).

---

## Table of Contents

1. [What you need](#1-what-you-need)
2. [Part A -- Set up WSL (one-time)](#2-part-a----set-up-wsl-one-time)
3. [Part B -- Get the pipeline files into WSL](#3-part-b----get-the-pipeline-files-into-wsl)
4. [Part C -- Install the Python environment (one-time)](#4-part-c----install-the-python-environment-one-time)
5. [Part D -- Run the pipeline](#5-part-d----run-the-pipeline)
6. [Part E -- Explore results / use the notebook](#6-part-e----explore-results--use-the-notebook)
7. [Switching datasets](#7-switching-datasets)
8. [What each stage does (concept -> code map)](#8-what-each-stage-does-concept----code-map)
9. [Troubleshooting](#9-troubleshooting)
10. [Sources](#10-sources)

---

## 1) What you need

- A Windows 10 (build 19041+) or Windows 11 PC.
- Administrator access (needed once, to enable WSL).
- An internet connection (to download WSL, Ubuntu, conda packages, and the datasets).
- About 6 GB of free disk space (conda environment + downloaded datasets).

---

## 2) Part A -- Set up WSL (one-time)

WSL lets you run a real Ubuntu Linux terminal directly inside Windows. The
bioinformatics tools this pipeline uses (Scanpy and its dependencies) are
built and tested for Linux, so WSL is the smoothest way to run them on a
Windows PC.

### 2.1 Check whether WSL is already installed

Open **PowerShell** (Start menu -> type `PowerShell` -> Enter) and run:

```powershell
wsl --status
```

- If it prints a WSL version and a default distribution name, WSL is already
  installed -- skip to [2.3](#23-install-ubuntu).
- If it says something like `Class not registered` or
  `Wsl/CallMsi/Install/REGDB_E_CLASSNOTREG`, WSL is present as a Windows
  feature but not correctly registered -- go to
  [Troubleshooting: "Class not registered"](#class-not-registered-error-on-wsl---status)
  before continuing.
- If the command is not recognized at all, install WSL as below.

### 2.2 Install WSL

In an **elevated (Run as Administrator)** PowerShell window:

```powershell
wsl --install
```

This single command enables the required Windows features (Windows
Subsystem for Linux + Virtual Machine Platform), downloads the WSL2 Linux
kernel, and installs **Ubuntu** as the default distribution. Reboot when
prompted.

Official reference: <https://learn.microsoft.com/en-us/windows/wsl/install>

### 2.3 Install Ubuntu

If `wsl --install` didn't already set up a distribution (e.g. WSL was
present but empty), install Ubuntu explicitly:

```powershell
wsl --install -d Ubuntu
```

The first launch asks you to create a UNIX username and password for
*inside* WSL -- this is separate from your Windows login and can be
anything you like. Remember the password; you'll need it for `sudo`
commands.

### 2.4 Confirm it works

```powershell
wsl -l -v
```

You should see `Ubuntu` listed with `VERSION 2` and `STATE Running` (or
`Stopped`, which is fine -- it starts on demand). Open it with:

```powershell
wsl
```

You're now at a Linux `$` prompt inside Windows.

---

## 3) Part B -- Get the pipeline files into WSL

This pipeline lives in your Windows filesystem under:

```
C:\Users\<you>\OneDrive\Documents\GauravAcademics\NGS\scRNA_Teaching_Pipeline
```

WSL can see your whole Windows C: drive under `/mnt/c/...`. From your WSL
terminal:

```bash
cd "/mnt/c/Users/<you>/OneDrive/Documents/GauravAcademics/NGS/scRNA_Teaching_Pipeline"
ls
```

Replace `<you>` with your actual Windows username. You should see
`environment.yml`, `config.yaml`, `run_pipeline.sh`, `scripts/`,
`notebooks/`, `setup/`, etc.

> **Tip:** running from a path under OneDrive works fine, but OneDrive's
> background syncing can occasionally slow down large file writes. If you
> notice it's sluggish, you can instead copy the whole
> `scRNA_Teaching_Pipeline` folder to a plain Linux path such as
> `~/scRNA_Teaching_Pipeline` (run `cp -r "/mnt/c/Users/<you>/OneDrive/Documents/GauravAcademics/NGS/scRNA_Teaching_Pipeline" ~/`)
> and work from there instead -- it will run noticeably faster either way,
> since native Linux filesystem I/O is faster than the `/mnt/c/...` bridge.

---

## 4) Part C -- Install the Python environment (one-time)

From inside the `scRNA_Teaching_Pipeline` folder in your WSL terminal:

```bash
bash setup/install_env.sh
```

This script (safe to re-run):
1. Installs **Miniforge** (a lightweight conda distribution) into
   `~/miniforge3` if conda isn't already available.
2. Creates a conda environment called **`scrna-teaching`** from
   `environment.yml` (Scanpy, leiden clustering, Scrublet, Harmony, BBKNN,
   CellTypist, JupyterLab, etc.). This step downloads several hundred MB and
   can take 5-15 minutes depending on your connection.
3. Registers `scrna-teaching` as a Jupyter kernel.

If this is the very first time conda has been installed, close and reopen
your WSL terminal (or run `source ~/.bashrc`) once the script finishes, so
the `conda` command is available in new terminals.

---

## 5) Part D -- Run the pipeline

Activate the environment, then run everything:

```bash
conda activate scrna-teaching
bash run_pipeline.sh
```

This runs all 11 stages in order (see
[section 8](#8-what-each-stage-does-concept----code-map)), printing progress
and explanations as it goes -- by default on **pbmc3k**, a real 10x
Genomics PBMC dataset that scanpy downloads automatically on stage 1.

Runtime: a few minutes on a normal laptop (pbmc3k has ~2,700 cells).

You can also run any single stage on its own, once earlier stages have run
at least once:

```bash
python scripts/08_clustering_de.py
```

### Outputs

- `data/processed/*.h5ad` -- the AnnData object after each stage (numbered,
  so you can load any intermediate stage in a notebook to inspect it).
- `results/figures/*.png` -- QC violin plots, UMAP plots, marker dotplot.
- `results/tables/*.csv` -- marker genes per cluster, manual annotation
  scores, and (for the GitHub sample dataset) a clustering comparison table.

---

## 6) Part E -- Explore results / use the notebook

To open the interactive walkthrough notebook (mirrors the slide deck
section by section):

```bash
conda activate scrna-teaching
jupyter lab
```

This prints a URL like `http://localhost:8888/lab?token=...`. WSL2
automatically forwards `localhost` ports to Windows, so you can simply
**copy that URL into your normal Windows web browser** (Edge, Chrome,
etc.) -- no extra networking setup required. In JupyterLab, open
`notebooks/scRNA_teaching_walkthrough.ipynb`, select the
**Python (scrna-teaching)** kernel if it isn't already selected, and run
the cells top to bottom (Shift+Enter).

---

## 7) Switching datasets

Open `config.yaml` and change:

```yaml
dataset: "pbmc3k"          # real 10x Genomics data -- every teaching step applies
# or
dataset: "github_sample"   # the sample dataset from the course GitHub repo
```

Then re-run `bash run_pipeline.sh` (or delete `data/processed/*.h5ad` first
if you want a completely clean run).

The **GitHub sample** (`GiatrasKon/scRNAseq-Analysis-Pipeline`,
`data/RNA-seq.csv.gz`) is 137 cells x 54,675 Affymetrix microarray probes,
**already log2-normalized** -- not raw 10x UMI counts, and not labeled with
human gene symbols. The pipeline detects this automatically
(`scripts/utils.py`) and skips mitochondrial/ribosomal/hemoglobin QC,
Scrublet doublet detection, and re-normalization for this dataset, printing
an explanation instead of producing meaningless numbers. Everything else
(HVG, PCA, neighbors/UMAP, Leiden clustering, differential expression,
manual annotation) still runs. This dataset also ships ground-truth cluster
labels from the original course assignment
(`data/labels.csv`), which stage 8 automatically uses to score Leiden's
clusters (Adjusted Rand Index / Normalized Mutual Information) against both
the ground truth and the original assignment's own GMM/DBSCAN clustering
approach -- a direct comparison between the classical-ML method taught in
that assignment and the modern graph-based Leiden method taught in this
course.

---

## 8) What each stage does (concept -> code map)

| Stage | Script | Slide concept |
|---|---|---|
| 1 | `01_fetch_data.py` | "From FASTQ to Gene-Cell Matrix" -- downloads the count matrix |
| 2 | `02_build_anndata.py` | "Different File Formats" / "What is AnnData?" |
| 3 | `03_qc_metrics.py` | "Why QC is Needed", mitochondrial/ribo/hb genes, MAD outliers |
| 4 | `04_doublet_detection.py` | "What is a Doublet?" -- Scrublet |
| 5 | `05_normalize_hvg.py` | "Normalization", "Log Transformation", "Highly Variable Genes" |
| 6 | `06_scale_pca.py` | "What is Scaling?", "What is PCA?" |
| 7 | `07_neighbors_umap.py` | "Why Neighbors?", "What is UMAP?" |
| 8 | `08_clustering_de.py` | Leiden clustering, "What is Resolution?", `rank_genes_groups` |
| 9 | `09_batch_correction.py` | "What is Batch Effect?" -- BBKNN / Harmony |
| 10 | `10_annotation.py` | "What is Cell Annotation?" -- manual markers + CellTypist |
| 11 | `11_generate_figures.py` | Saves the plots used throughout the notebook |

---

## 9) Troubleshooting

### "Class not registered" error on `wsl --status`

Full error text looks like:

```
Class not registered
Error code: Wsl/CallMsi/Install/REGDB_E_CLASSNOTREG
```

This means the WSL Windows feature is present but its installer component
isn't correctly registered with Windows -- a known WSL issue, not something
specific to this pipeline. Fix, in order (try each, moving on only if the
previous one doesn't resolve it):

1. **Update WSL** (elevated PowerShell):
   ```powershell
   wsl --update
   wsl --shutdown
   ```
   Then retry `wsl --status`.
2. **Re-enable the Windows features** (elevated PowerShell), then reboot:
   ```powershell
   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
   ```
   Reboot, then run `wsl --install` again.
3. **Reinstall the "Windows Subsystem for Linux" app** from the Microsoft
   Store (search "Windows Subsystem for Linux"), which carries its own
   updated installer/registration logic independent of the OS image.
4. **Run the System File Checker**, in case a system DLL/COM registration
   is corrupted (elevated PowerShell):
   ```powershell
   sfc /scannow
   DISM /Online /Cleanup-Image /RestoreHealth
   ```
   Reboot and retry.

Full official troubleshooting reference:
<https://learn.microsoft.com/en-us/windows/wsl/troubleshooting>

### `conda: command not found` in a new terminal after setup

Run `source ~/.bashrc`, or simply close and reopen the WSL terminal. If it
still isn't found, `~/miniforge3/bin/conda init bash` and reopen the
terminal again.

### `bash run_pipeline.sh` says the environment isn't active

Run `conda activate scrna-teaching` first. If the environment doesn't
exist yet, run `bash setup/install_env.sh`.

### `environment.yml` install is very slow or hangs

Conda's own dependency solver can be slow. Miniforge (which this guide
uses) ships **mamba** in newer releases -- if `conda env create` is
taking a very long time, cancel it and instead run:
```bash
conda install -n base -c conda-forge mamba -y
mamba env create -f environment.yml
```
`mamba` is a drop-in, much faster reimplementation of the same command.

### Out of memory / the process gets killed

pbmc3k (~2,700 cells) comfortably fits in a few hundred MB of RAM; the
GitHub sample dataset (137 cells) is even smaller. If you still hit memory
issues (e.g. on a machine with very little RAM allocated to WSL), create
or edit `C:\Users\<you>\.wslconfig` with:
```ini
[wsl2]
memory=6GB
```
then run `wsl --shutdown` in PowerShell and reopen WSL.

### JupyterLab opens but the page doesn't load in the Windows browser

Copy the *exact* URL (including the `?token=...` part) printed in the
terminal. If `localhost` doesn't work, try replacing it with `127.0.0.1` in
the URL.

### `celltypist` / model download fails (no internet in WSL, or blocked)

CellTypist annotation is optional (`annotation.run_celltypist: false` by
default in `config.yaml`). It downloads its model file on first use from
the Sanger Institute's `celltypist.cog.sanger.ac.uk` file server (model
catalog: <https://www.celltypist.org/models>). If your network blocks that
domain, leave `run_celltypist: false` -- the rest of the pipeline does not
depend on it.

---

## 10) Sources

- Course slides: `SC_analysis.pptx`
- Course GitHub repository:
  <https://github.com/GiatrasKon/scRNAseq-Analysis-Pipeline>
- Scanpy documentation: <https://scanpy.readthedocs.io/>
- WSL installation docs: <https://learn.microsoft.com/en-us/windows/wsl/install>
- WSL troubleshooting docs: <https://learn.microsoft.com/en-us/windows/wsl/troubleshooting>
- Miniforge (conda distribution): <https://github.com/conda-forge/miniforge>
- 10x Genomics pbmc3k dataset (downloaded automatically by scanpy):
  <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>
- Scrublet: <https://pypi.org/project/scrublet/>
- Harmony (harmonypy): <https://pypi.org/project/harmonypy/>
- BBKNN: <https://pypi.org/project/bbknn/>
- CellTypist: <https://www.celltypist.org/> / <https://pypi.org/project/celltypist/>
