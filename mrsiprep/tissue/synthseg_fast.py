"""SynthSeg-constrained FSL FAST tissue backend."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy import ndimage as ndi

from mrsiprep.interfaces.fsl import run_fast
from mrsiprep.interfaces.hdbet import run_hd_bet
from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti


CSF_VENTRICLE_LABELS = {4, 5, 14, 15, 24, 43, 44}


def synthseg_fast_pve_path(config, subject: str, session: str | None, index: int) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, desc=f"p{index}")


def synthseg_fast_brain_path(config, subject: str, session: str | None) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="hdbetSynthsegFastBrain")


def synthseg_fast_brain_mask_path(config, subject: str, session: str | None) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="hdbetSynthsegFastBrain", suffix_override="mask")


def segment_t1_synthseg_fast(config, subject: str, session: str | None, t1_path: Path) -> dict[str, Path]:
    """Segment GM/WM/CSF with HD-BET + SynthSeg-constrained FAST.

    SynthSeg contributes a CSF/ventricle prior and HD-BET contributes the
    brain mask. The combined mask is applied to the raw T1 before FAST, and
    FAST supplies the partial-volume estimates.
    """

    t1_path = Path(t1_path)
    work_dir = config.work_dir / f"sub-{subject}" / (f"ses-{session}" if session else "ses-none") / "synthseg_fast"
    work_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "GM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="GM", suffix_override="probseg"),
        "WM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="WM", suffix_override="probseg"),
        "CSF": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="CSF", suffix_override="probseg"),
    }
    pve_outputs = {
        "GM": synthseg_fast_pve_path(config, subject, session, 1),
        "WM": synthseg_fast_pve_path(config, subject, session, 2),
        "CSF": synthseg_fast_pve_path(config, subject, session, 3),
    }
    brain_out = synthseg_fast_brain_path(config, subject, session)
    brain_mask_out = synthseg_fast_brain_mask_path(config, subject, session)
    if all(path.exists() for path in [*outputs.values(), *pve_outputs.values(), brain_out, brain_mask_out]) and not config.overwrite:
        return outputs

    t1_img = nib.load(str(t1_path))
    labels_native = _run_or_load_synthseg(config, t1_path, work_dir)
    csf_vent_mask = _synthseg_csf_ventricle_mask(labels_native)
    hdbet_mask = _run_hdbet_mask(config, t1_path, brain_out, brain_mask_out)
    fast_mask = hdbet_mask | csf_vent_mask
    masked_t1 = _write_masked_t1(t1_img, fast_mask, work_dir / "synthsegFast_maskedT1w.nii.gz")

    fast_sources = run_fast(masked_t1, work_dir / "fast" / "fast", verbose=config.verbose)
    fast_maps = {}
    for label, source in fast_sources.items():
        img = nib.load(str(source))
        fast_maps[label] = np.clip(np.nan_to_num(img.get_fdata(dtype=np.float32).squeeze(), copy=False), 0.0, 1.0)
    fast_maps = _apply_synthseg_csf_tissue_correction(fast_maps, labels_native)
    for label, data in fast_maps.items():
        outputs[label] = save_nifti(data, t1_img, outputs[label], dtype=np.float32)
        save_nifti(data, t1_img, pve_outputs[label], dtype=np.float32)
    return outputs


def _run_or_load_synthseg(config, t1_path: Path, work_dir: Path) -> np.ndarray:
    native_labels = work_dir / "synthseg_labels_space-T1w.nii.gz"
    if native_labels.exists() and not config.overwrite:
        return _load_labels(native_labels)

    synthseg_labels = work_dir / "synthseg_labels.nii.gz"
    if not synthseg_labels.exists() or config.overwrite:
        synthseg_cmd = _find_mri_synthseg()
        subprocess.run(
            [
                synthseg_cmd,
                "--i",
                str(t1_path),
                "--o",
                str(synthseg_labels),
                "--robust",
            ],
            check=True,
            stdout=None if config.verbose else subprocess.PIPE,
            stderr=None if config.verbose else subprocess.PIPE,
            text=True,
            env=_synthseg_env(synthseg_cmd),
        )

    labels_img = nib.load(str(synthseg_labels))
    t1_img = nib.load(str(t1_path))
    if labels_img.shape[:3] != t1_img.shape[:3] or not np.allclose(labels_img.affine, t1_img.affine, atol=1e-3):
        labels_img = resample_from_to(labels_img, (t1_img.shape[:3], t1_img.affine), order=0)
    labels = np.rint(np.nan_to_num(labels_img.get_fdata(dtype=np.float32).squeeze(), copy=False)).astype(np.uint8)
    save_nifti(labels, t1_img, native_labels, dtype=np.uint8)
    return labels


def _run_hdbet_mask(config, t1_path: Path, brain_out: Path, mask_out: Path) -> np.ndarray:
    if not brain_out.exists() or config.overwrite:
        hdbet_input = t1_path
        temp_dir = None
        tmp_input = None
        if t1_path.name.endswith(".nii.gz"):
            temp_dir = tempfile.TemporaryDirectory(prefix="mrsiprep_synthseg_fast_hdbet_")
            tmp_input = Path(temp_dir.name) / t1_path.name[:-7]
            nib.save(nib.load(str(t1_path)), str(tmp_input))
            hdbet_input = tmp_input
        try:
            run_hd_bet(hdbet_input, brain_out, device="cuda", verbose=config.verbose)
        finally:
            if tmp_input is not None and tmp_input.exists():
                tmp_input.unlink()
            if temp_dir is not None:
                temp_dir.cleanup()

    brain_img = nib.load(str(brain_out))
    brain = np.nan_to_num(brain_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    mask = ndi.binary_fill_holes(brain > 0)
    save_nifti(mask.astype(np.uint8), brain_img, mask_out, dtype=np.uint8)
    return mask.astype(bool, copy=False)


def _synthseg_csf_ventricle_mask(labels: np.ndarray) -> np.ndarray:
    return np.isin(np.asarray(labels), list(CSF_VENTRICLE_LABELS))


def _apply_synthseg_csf_tissue_correction(fast_maps: dict[str, np.ndarray], labels: np.ndarray) -> dict[str, np.ndarray]:
    corrected = {label: data.copy() for label, data in fast_maps.items()}
    synthseg_labels = np.asarray(labels)
    synthseg_background = synthseg_labels == 0
    if np.any(synthseg_background):
        for label in ("GM", "WM", "CSF"):
            corrected[label][synthseg_background] = 0.0

    synthseg_csf = synthseg_labels == 24
    tissue_in_synthseg_csf = synthseg_csf & ((corrected["GM"] > 0.01) | (corrected["WM"] > 0.01))
    if np.any(tissue_in_synthseg_csf):
        corrected["GM"][tissue_in_synthseg_csf] = 0.0
        corrected["WM"][tissue_in_synthseg_csf] = 0.0
        corrected["CSF"][tissue_in_synthseg_csf] = 1.0
    return corrected


def _write_masked_t1(t1_img: nib.Nifti1Image, mask: np.ndarray, out_path: Path) -> Path:
    t1 = np.nan_to_num(t1_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    masked = np.where(mask, t1, 0.0).astype(np.float32)
    return save_nifti(masked, t1_img, out_path, dtype=np.float32)


def _load_labels(path: Path) -> np.ndarray:
    return np.rint(np.nan_to_num(nib.load(str(path)).get_fdata(dtype=np.float32).squeeze(), copy=False)).astype(np.uint8)


def _find_mri_synthseg() -> str:
    found = shutil.which("mri_synthseg")
    if found:
        return found

    candidates: list[Path] = []
    local_apps = Path.home() / "Apps"
    candidates.extend(sorted(local_apps.glob("freesurfer-*/*/bin/mri_synthseg"), reverse=True))
    candidates.extend(sorted(local_apps.glob("freesurfer-*/*/python/scripts/mri_synthseg"), reverse=True))
    fs_home = os.environ.get("FREESURFER_HOME")
    if fs_home:
        root = Path(fs_home)
        candidates.extend(
            [
                root / "bin" / "mri_synthseg",
                root / "python" / "scripts" / "mri_synthseg",
            ]
        )
        candidates.extend(root.glob("*/bin/mri_synthseg"))
        candidates.extend(root.glob("*/python/scripts/mri_synthseg"))
    candidates.extend(
        [
            Path("/usr/local/freesurfer/8.2.0/bin/mri_synthseg"),
            Path("/usr/local/freesurfer/8.2.0/python/scripts/mri_synthseg"),
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("mri_synthseg was not found on PATH or in FREESURFER_HOME. Source FreeSurfer before running synthseg-fast.")


def _synthseg_env(command: str) -> dict[str, str]:
    env = os.environ.copy()
    path = Path(command)
    if path.name == "mri_synthseg" and path.parent.name == "bin":
        fs_home = path.parent.parent
        if (fs_home / "python" / "scripts" / "mri_synthseg").exists():
            env["FREESURFER_HOME"] = str(fs_home)
    return env
