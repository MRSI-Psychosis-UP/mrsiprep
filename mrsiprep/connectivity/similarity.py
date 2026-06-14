"""Similarity matrix computation."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr


def compute_similarity(matrix: pd.DataFrame, method: str = "spearman") -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    labels = matrix.index
    n = values.shape[0]
    out = np.eye(n, dtype=float)
    if method == "cosine":
        norms = np.linalg.norm(values, axis=1, keepdims=True)
        normed = np.divide(values, norms, out=np.zeros_like(values), where=norms > 0)
        out = normed @ normed.T
    elif method == "euclidean_distance":
        out = squareform(pdist(np.nan_to_num(values), metric="euclidean"))
    elif method in {"pearson", "spearman"}:
        for i in range(n):
            for j in range(i + 1, n):
                x = values[i]
                y = values[j]
                valid = np.isfinite(x) & np.isfinite(y)
                if valid.sum() < 2:
                    val = np.nan
                elif method == "pearson":
                    val = pearsonr(x[valid], y[valid]).statistic
                else:
                    val = spearmanr(x[valid], y[valid]).statistic
                out[i, j] = out[j, i] = val
    else:
        raise ValueError(f"Unsupported connectivity method: {method}")
    return pd.DataFrame(out, index=labels, columns=labels)
