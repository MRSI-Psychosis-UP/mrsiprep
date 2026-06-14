"""Biharmonic-style MRSI spike filtering."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn import image as nil_image
from skimage.restoration import inpaint_biharmonic

from mrsiprep.io.naming import mrsi_derivative
from mrsiprep.utils.images import save_nifti


def filter_metabolite_maps(config, subject: str, session: str | None, metabolite_maps: dict[str, Path], brainmask: Path) -> dict[str, Path]:
    if not config.filter_biharmonic:
        return metabolite_maps
    brain = nib.load(str(brainmask)).get_fdata().astype(bool)
    filtered: dict[str, Path] = {}
    for met, path in metabolite_maps.items():
        out = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", met=met, desc="preproc", suffix_override="mrsi")
        if out.exists() and not (config.overwrite_filt or config.overwrite):
            filtered[met] = out
            continue
        img = nib.load(str(path))
        data = np.nan_to_num(img.get_fdata(dtype=np.float32), nan=0.0).squeeze()
        spike_mask = get_spike_mask(data, percentile=config.spike_percentile)
        missing = ((data == 0) & brain) | spike_mask
        repaired = data.copy()
        if np.any(missing):
            try:
                repaired = inpaint_biharmonic(repaired, missing)
            except Exception:
                repaired[missing] = np.nanmedian(data[brain & ~missing])
        voxel_dims = np.array(img.header.get_zooms()[:3])
        fwhm = float(np.round(voxel_dims.mean() * np.sqrt(2)))
        smooth = nil_image.smooth_img(nib.Nifti1Image(repaired.astype(np.float32), img.affine, img.header), fwhm=fwhm).get_fdata()
        repaired[missing] = smooth[missing]
        repaired[~brain] = 0
        filtered[met] = save_nifti(repaired.astype(np.float32), img, out, dtype=np.float32)
        spike_out = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", met=met, desc="spikemask", suffix_override="mask")
        save_nifti(spike_mask.astype(np.uint8), img, spike_out, dtype=np.uint8)
    return filtered


def get_spike_mask(data: np.ndarray, percentile: float = 99.0) -> np.ndarray:
    inside = data > 0
    if not np.any(inside):
        return np.zeros_like(data, dtype=bool)
    threshold = np.percentile(data[inside], percentile)
    return data > threshold
