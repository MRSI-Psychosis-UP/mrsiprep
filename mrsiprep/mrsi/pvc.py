"""Partial-volume correction via PETPVC."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.io.naming import mrsi_derivative
from mrsiprep.utils.images import save_nifti


class PVCError(RuntimeError):
    """Raised when partial-volume correction fails."""


def create_tissue_4d(config, subject: str, session: str | None, tissue_mrsi: dict[str, Path], reference: Path) -> Path:
    out = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", desc="4Dtissue", suffix_override="mrsi")
    if out.exists() and not (config.overwrite_pve or config.overwrite):
        return out
    ref_img = nib.load(str(reference))
    data = np.stack([nib.load(str(tissue_mrsi[label])).get_fdata() for label in ("GM", "WM", "CSF")], axis=-1)
    return save_nifti(data.astype(np.float32), ref_img, out, dtype=np.float32)


def run_pvc(config, subject: str, session: str | None, metabolite_maps: dict[str, Path], tissue_4d: Path, brainmask: Path, psf_width: float = 5.0) -> dict[str, Path]:
    if shutil.which("petpvc") is None:
        raise PVCError("petpvc command not found on PATH. Use --no-pvc to skip partial-volume correction.")
    brain = nib.load(str(brainmask)).get_fdata().astype(bool)
    out_maps: dict[str, Path] = {}
    for met, path in metabolite_maps.items():
        out = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", met=met, desc="pvc", suffix_override="mrsi")
        if out.exists() and not (config.overwrite_pve or config.overwrite):
            out_maps[met] = out
            continue
        tmp_out = out.with_name(out.name.replace("_desc-pvc_", "_desc-petpvcraw_"))
        cmd = ["petpvc", "-i", str(path), "-m", str(tissue_4d), "-p", "RBV", "-x", str(psf_width), "-y", str(psf_width), "-z", str(psf_width), "-o", str(tmp_out)]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        img = nib.load(str(tmp_out))
        raw = nib.load(str(path)).get_fdata()
        data = img.get_fdata()
        data[data > 2 * raw] = 0
        data[data < 0] = 0
        data[~brain] = 0
        out_maps[met] = save_nifti(data.astype(np.float32), img, out, dtype=np.float32)
    return out_maps
