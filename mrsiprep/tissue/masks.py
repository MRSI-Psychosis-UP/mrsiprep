"""Tissue segmentation masks."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti


def build_brain_csf_seed_mask(
    config,
    subject: str,
    session: str | None,
    skull_t1: Path,
    raw_t1: Path,
    brain_mask: Path | None = None,
) -> Path:
    """Build an intracranial mask seed for Atropos.

    The mask is derived from skull-stripping geometry only. It is not derived
    from CAT12 p1/p2/p3 maps, so Atropos can use it to create those maps.
    """
    out_mask = anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="brainCSFseed", suffix_override="mask")
    if out_mask.exists() and not config.overwrite:
        return out_mask

    raw_img = nib.load(str(raw_t1))
    seed_img = nib.load(str(brain_mask)) if brain_mask and Path(brain_mask).exists() else nib.load(str(skull_t1))
    if raw_img.shape[:3] != seed_img.shape[:3]:
        raise ValueError("Cannot build Atropos brain+CSF mask: raw T1 and skull-strip seed have different shapes.")
    if not np.allclose(raw_img.affine, seed_img.affine, atol=1e-3):
        raise ValueError("Cannot build Atropos brain+CSF mask: raw T1 and skull-strip seed do not share the same affine.")

    seed = np.nan_to_num(seed_img.get_fdata(dtype=np.float32).squeeze(), copy=False) > 0
    if not np.any(seed):
        raise ValueError("Cannot build Atropos brain+CSF mask: skull-strip seed is empty.")

    mask = ndi.binary_fill_holes(seed)
    mask = ndi.binary_closing(mask, structure=ndi.generate_binary_structure(3, 2), iterations=1)
    mask = _largest_component(mask)
    iterations = _dilation_iterations(raw_img, config.atropos_mask_dilation_mm)
    if iterations > 0:
        mask = ndi.binary_dilation(mask, structure=ndi.generate_binary_structure(3, 2), iterations=iterations)
    mask = ndi.binary_closing(mask, structure=ndi.generate_binary_structure(3, 2), iterations=1)
    mask = ndi.binary_fill_holes(mask)
    mask = _largest_component(mask)
    return save_nifti(mask.astype(np.uint8), raw_img, out_mask, dtype=np.uint8)


def _dilation_iterations(img: nib.Nifti1Image, dilation_mm: float) -> int:
    if dilation_mm <= 0:
        return 0
    min_zoom = min(float(z) for z in img.header.get_zooms()[:3])
    if min_zoom <= 0:
        return int(round(dilation_mm))
    return max(1, int(round(dilation_mm / min_zoom)))


def _largest_component(mask: np.ndarray) -> np.ndarray:
    labels, n_labels = ndi.label(mask)
    if n_labels == 0:
        return mask.astype(bool)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return labels == int(np.argmax(counts))
