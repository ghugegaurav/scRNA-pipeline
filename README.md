# scRNA-seq Teaching Pipeline

A Scanpy-based single-cell RNA-seq analysis pipeline built to teach the
workflow covered in `SC_analysis.pptx`: quality control, doublet detection,
normalization, dimensionality reduction, clustering, batch correction, and
cell type annotation. Designed to run inside **WSL (Windows Subsystem for
Linux)**.

**New here? Start with [`GUIDE.md`](GUIDE.md)** -- it walks through setting
up WSL, installing the environment, and running the pipeline, step by step.

## Layout

```
scRNA_Teaching_Pipeline/
├── GUIDE.md                 <- full setup + execution guide (start here)
├── config.yaml               <- all tunable parameters live here
├── environment.yml           <- conda environment definition
├── run_pipeline.sh           <- runs every stage in order
├── setup/install_env.sh      <- one-time WSL/conda setup script
├── scripts/                  <- the 11 pipeline stages (01_ ... 11_)
├── notebooks/                <- interactive walkthrough notebook
├── data/                     <- raw + processed data (created at runtime)
└── results/                  <- figures + tables (created at runtime)
```

## Two datasets, one pipeline

Set `dataset:` in `config.yaml` to:
- `"pbmc3k"` (default) -- a real 10x Genomics PBMC dataset, downloaded
  automatically. Every step taught in the slides applies directly.
- `"github_sample"` -- the sample dataset from
  [GiatrasKon/scRNAseq-Analysis-Pipeline](https://github.com/GiatrasKon/scRNAseq-Analysis-Pipeline),
  downloaded automatically from that repository. Raw-count-only steps
  (mitochondrial QC, Scrublet) are automatically skipped for this dataset
  with an explanation, since it is pre-normalized microarray data rather
  than raw droplet counts -- see `GUIDE.md` section 7 for details.

## Quick start (after WSL is set up -- see GUIDE.md)

```bash
bash setup/install_env.sh
conda activate scrna-teaching
bash run_pipeline.sh
```
