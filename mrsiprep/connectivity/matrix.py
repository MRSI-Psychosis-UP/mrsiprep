"""Build regional metabolite matrices."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_regional_matrix(regional_table: str | Path, value_col: str = "weighted_mean", min_coverage: float | None = None) -> pd.DataFrame:
    df = pd.read_csv(regional_table, sep="\t")
    if min_coverage is not None:
        df = df[df["coverage"] >= min_coverage]
    matrix = df.pivot_table(index="parcel_id", columns="metabolite", values=value_col, aggfunc="mean")
    return matrix.sort_index()


def zscore_columns(matrix: pd.DataFrame) -> pd.DataFrame:
    return (matrix - matrix.mean(axis=0)) / matrix.std(axis=0).replace(0, 1)
