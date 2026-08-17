"""
Shared helper functions used by every numbered stage script.

Kept deliberately small and readable -- students should be able to open this
file and understand every line, since several of these functions are the
direct code equivalents of concepts explained in the course slides.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_config(config_path: str | Path | None = None) -> dict:
    """Load config.yaml from the repo root (or a custom path)."""
    path = Path(config_path) if config_path else REPO_ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def step(title: str):
    """Print a visible section header so console output mirrors the slide deck's structure."""
    line = "=" * 70
    print(f"\n{line}\n{title}\n{line}")


def is_outlier(series: pd.Series, nmads: int = 5) -> pd.Series:
    """
    MAD-based outlier detection (course slide: 'What is MAD (Median Absolute
    Deviation)?'). Flags values further than `nmads` Median Absolute
    Deviations from the median -- robust to extreme values, unlike a
    standard-deviation cutoff.
    """
    from scipy.stats import median_abs_deviation

    med = np.median(series)
    mad = median_abs_deviation(series)
    if mad == 0:
        return pd.Series(np.zeros(len(series), dtype=bool), index=series.index)
    return (series < med - nmads * mad) | (series > med + nmads * mad)


def looks_like_raw_counts(matrix, sample_size: int = 5000) -> bool:
    """
    Heuristic used to decide whether count-based QC steps (mitochondrial %,
    Scrublet doublet detection) are scientifically meaningful for the loaded
    dataset. Raw UMI/read counts are non-negative integers; normalized or
    log-transformed data is not.
    """
    import scipy.sparse as sp

    if sp.issparse(matrix):
        sample = matrix[: min(sample_size, matrix.shape[0])].toarray()
    else:
        sample = np.asarray(matrix[: min(sample_size, matrix.shape[0])])
    if sample.size == 0:
        return False
    finite = sample[np.isfinite(sample)]
    if finite.size == 0:
        return False
    non_negative = (finite >= 0).all()
    near_integer = np.allclose(finite, np.round(finite), atol=1e-6)
    return bool(non_negative and near_integer)


def has_gene_symbol_annotation(var_names, prefixes=("MT-", "RPS", "RPL", "HB")) -> bool:
    """
    Checks whether `var_names` look like human gene symbols (as opposed to
    e.g. Affymetrix microarray probe IDs like '1007_s_at'), which is what
    the mitochondrial/ribosomal/hemoglobin prefix-matching QC step needs.
    """
    upper = pd.Index(var_names).astype(str).str.upper()
    return bool(upper.str.startswith(prefixes).any())
