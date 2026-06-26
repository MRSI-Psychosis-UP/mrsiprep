"""Minimal HTML reports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mrsiprep.io.naming import figure_derivative


def generate_subject_report(config, subject: str, session: str | None, outputs: dict) -> Path:
    out = figure_derivative(config.derivative_dir, subject, session, extension="html", desc="report")
    out.parent.mkdir(parents=True, exist_ok=True)
    qc_html = ""
    qc_path = outputs.get("qc_summary")
    if qc_path and Path(qc_path).exists():
        qc_html = pd.read_csv(qc_path, sep="\t").to_html(index=False, border=0)
    regional_html = ""
    regional = outputs.get("regional_table")
    if regional and Path(regional).exists():
        regional_html = pd.read_csv(regional, sep="\t").head(50).to_html(index=False, border=0)
    parcel_qc_html = ""
    parcel_qc_summary = ""
    parcel_qc = outputs.get("parcel_qc")
    if parcel_qc and Path(parcel_qc).exists():
        parcel_df = pd.read_csv(parcel_qc, sep="\t")
        overview = (
            parcel_df.groupby(["parcel_id", "parcel_name", "hemisphere"], dropna=False)
            .agg(
                anatomical_coverage_percent=("anatomical_coverage_percent", "first"),
                mean_crlb=("mean_crlb", "mean"),
                qc_valid_fraction=("qc_valid_fraction", "mean"),
            )
            .reset_index()
            .sort_values("anatomical_coverage_percent")
        )
        parcel_qc_html = overview.to_html(index=False, border=0, float_format=lambda value: f"{value:.3f}")
        parcel_qc_summary = (
            f"<p>Mean anatomical MRSI coverage: <strong>{overview['anatomical_coverage_percent'].mean():.1f}%</strong>; "
            f"mean parcel CRLB: <strong>{parcel_df['mean_crlb'].mean():.2f}</strong>.</p>"
        )
    lines = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'><title>MRSIPrep report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:2rem;line-height:1.4}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:4px 8px}code{background:#f3f3f3;padding:2px 4px}</style>",
        "</head><body>",
        f"<h1>MRSIPrep report: sub-{subject}" + (f" ses-{session}" if session else "") + "</h1>",
        "<h2>Inputs</h2>",
        f"<p>BIDS directory: <code>{config.bids_dir}</code></p>",
        f"<p>Output directory: <code>{config.derivative_dir}</code></p>",
        f"<p>Processing mode: <code>{config.processing_mode}</code></p>",
        "<h2>MRSI QC</h2>",
        qc_html or "<p>No QC table available.</p>",
        "<h2>Parcelwise Coverage and CRLB</h2>",
        parcel_qc_summary,
        parcel_qc_html or "<p>No parcelwise QC table available.</p>",
        "<h2>Regional Metabolites</h2>",
        regional_html or "<p>No regional table available.</p>",
        "<h2>Outputs</h2>",
        "<ul>",
    ]
    for key, value in sorted(outputs.items()):
        lines.append(f"<li><strong>{key}</strong>: <code>{value}</code></li>")
    lines.extend(["</ul>", "</body></html>"])
    out.write_text("\n".join(lines), encoding="utf-8")
    return out
