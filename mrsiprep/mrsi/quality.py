"""Voxelwise MRSI QC masks."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.io.naming import mrsi_derivative
from mrsiprep.utils.images import save_nifti
from mrsiprep.utils.tables import write_tsv


def make_quality_masks(config, subject: str, session: str | None, metabolite_maps: dict[str, Path], crlb_maps: dict[str, Path], snr_map: Path | None, linewidth_map: Path | None, brainmask: Path) -> tuple[dict[str, Path], Path]:
    brain = nib.load(str(brainmask)).get_fdata().astype(bool)
    snr = nib.load(str(snr_map)).get_fdata() if snr_map and snr_map.exists() else None
    linewidth = nib.load(str(linewidth_map)).get_fdata() if linewidth_map and linewidth_map.exists() else None
    qcmasks: dict[str, Path] = {}
    rows = []

    for met, path in metabolite_maps.items():
        img = nib.load(str(path))
        data = img.get_fdata()
        mask = np.isfinite(data) & brain
        crlb = nib.load(str(crlb_maps[met])).get_fdata() if met in crlb_maps and crlb_maps[met].exists() else None
        if snr is not None:
            mask &= snr >= config.snr_min
        if linewidth is not None:
            mask &= linewidth <= config.linewidth_max
        if crlb is not None:
            mask &= crlb <= config.crlb_max
        out = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", met=met, desc="qcmask", suffix_override="mask")
        qcmasks[met] = save_nifti(mask.astype(np.uint8), img, out, dtype=np.uint8)
        valid = mask & np.isfinite(data)
        rows.append(
            {
                "metabolite": met,
                "n_total_voxels": int(brain.sum()),
                "n_valid_voxels": int(valid.sum()),
                "valid_fraction": float(valid.sum() / max(brain.sum(), 1)),
                "mean_snr": _safe_mean(snr, valid),
                "median_snr": _safe_median(snr, valid),
                "mean_linewidth": _safe_mean(linewidth, valid),
                "median_linewidth": _safe_median(linewidth, valid),
                "mean_crlb": _safe_mean(crlb, valid),
                "median_crlb": _safe_median(crlb, valid),
            }
        )
    summary = mrsi_derivative(config.derivative_dir, subject, session, desc="mrsiqc", suffix_override="tsv")
    write_tsv(rows, summary)
    return qcmasks, summary


def _safe_mean(data, mask):
    if data is None or not np.any(mask):
        return np.nan
    return float(np.nanmean(data[mask]))


def _safe_median(data, mask):
    if data is None or not np.any(mask):
        return np.nan
    return float(np.nanmedian(data[mask]))
