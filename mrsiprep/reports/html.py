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
    lines = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'><title>MRSIPrep report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:2rem;line-height:1.4}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:4px 8px}code{background:#f3f3f3;padding:2px 4px}</style>",
        "</head><body>",
        f"<h1>MRSIPrep report: sub-{subject}" + (f" ses-{session}" if session else "") + "</h1>",
        "<h2>Inputs</h2>",
        f"<p>BIDS directory: <code>{config.bids_dir}</code></p>",
        f"<p>Output directory: <code>{config.derivative_dir}</code></p>",
        "<h2>MRSI QC</h2>",
        qc_html or "<p>No QC table available.</p>",
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
