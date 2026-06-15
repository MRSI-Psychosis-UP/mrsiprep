"""FreeSurfer tissue backend."""

from __future__ import annotations

import shutil
from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.interfaces.freesurfer import convert_to_t1_space, freesurfer_subject_id, run_recon_all
from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti


GM_LABELS = {
    3,
    8,
    10,
    11,
    12,
    13,
    17,
    18,
    26,
    28,
    42,
    47,
    49,
    50,
    51,
    52,
    53,
    54,
    58,
    60,
}
WM_LABELS = {2, 7, 41, 46, 77, 251, 252, 253, 254, 255}
CSF_LABELS = {4, 5, 14, 15, 24, 31, 43, 44, 63, 72}


def freesurfer_pve_path(config, subject: str, session: str | None, index: int) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, desc=f"p{index}")


def freesurfer_brain_path(config, subject: str, session: str | None) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="brain")


def freesurfer_brain_mask_path(config, subject: str, session: str | None) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="brain", suffix_override="mask")


def freesurfer_aseg_path(config, subject: str, session: str | None) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, space="T1w", desc="freesurferAseg", suffix_override="dseg")


def segment_t1_freesurfer(config, subject: str, session: str | None, t1_path: Path, raw_t1_path: Path | None = None) -> dict[str, Path]:
    raw_t1 = Path(raw_t1_path or t1_path)
    fs_subject = freesurfer_subject_id(raw_t1)
    fs_dir = run_recon_all(
        raw_t1,
        config.freesurfer_dir,
        fs_subject,
        force=config.overwrite,
        nthreads=config.nthreads,
    )
    _write_native_brain_derivatives(config, subject, session, raw_t1, fs_dir)
    outputs = _write_native_tissue_maps(config, subject, session, raw_t1, fs_dir)
    return outputs


def _write_native_brain_derivatives(config, subject: str, session: str | None, raw_t1: Path, fs_dir: Path) -> None:
    brain_out = freesurfer_brain_path(config, subject, session)
    mask_out = freesurfer_brain_mask_path(config, subject, session)
    if brain_out.exists() and mask_out.exists() and not config.overwrite:
        return

    brain_mgz = fs_dir / "mri" / "brain.mgz"
    if not brain_mgz.exists():
        raise FileNotFoundError(f"FreeSurfer brain.mgz not found: {brain_mgz}")
    convert_to_t1_space(brain_mgz, raw_t1, brain_out, interpolation="trilin", verbose=config.verbose)
    brain_img = nib.load(str(brain_out))
    brain = np.nan_to_num(brain_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    save_nifti((brain > 0).astype(np.uint8), brain_img, mask_out, dtype=np.uint8)


def _write_native_tissue_maps(config, subject: str, session: str | None, raw_t1: Path, fs_dir: Path) -> dict[str, Path]:
    outputs = {
        "GM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="GM", suffix_override="probseg"),
        "WM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="WM", suffix_override="probseg"),
        "CSF": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="CSF", suffix_override="probseg"),
    }
    pve_outputs = {
        "GM": freesurfer_pve_path(config, subject, session, 1),
        "WM": freesurfer_pve_path(config, subject, session, 2),
        "CSF": freesurfer_pve_path(config, subject, session, 3),
    }
    aseg_out = freesurfer_aseg_path(config, subject, session)
    if all(path.exists() for path in [*outputs.values(), *pve_outputs.values(), aseg_out]) and not config.overwrite:
        return outputs

    aseg_mgz = fs_dir / "mri" / "aseg.mgz"
    if not aseg_mgz.exists():
        raise FileNotFoundError(f"FreeSurfer aseg.mgz not found: {aseg_mgz}")
    convert_to_t1_space(aseg_mgz, raw_t1, aseg_out, interpolation="nearest", verbose=config.verbose)
    fs_masks = _write_freesurfer_space_class_masks(fs_dir, aseg_mgz)

    for label, fs_mask in fs_masks.items():
        convert_to_t1_space(fs_mask, raw_t1, outputs[label], interpolation="trilin", verbose=config.verbose)

    ref = nib.load(str(raw_t1))
    probs = {label: _load_probability(path) for label, path in outputs.items()}
    total = probs["GM"] + probs["WM"] + probs["CSF"]
    overfull = total > 1.0
    for label in probs:
        data = probs[label]
        data[overfull] = data[overfull] / total[overfull]
        save_nifti(data, ref, outputs[label], dtype=np.float32)
        save_nifti(data, ref, pve_outputs[label], dtype=np.float32)
    return outputs


def _write_freesurfer_space_class_masks(fs_dir: Path, aseg_mgz: Path) -> dict[str, Path]:
    aseg_img = nib.load(str(aseg_mgz))
    aseg = np.rint(aseg_img.get_fdata(dtype=np.float32)).astype(np.int16)
    specs = {
        "GM": (GM_LABELS, fs_dir / "mri" / "mrsiprep_p1_gm.mgz"),
        "WM": (WM_LABELS, fs_dir / "mri" / "mrsiprep_p2_wm.mgz"),
        "CSF": (CSF_LABELS, fs_dir / "mri" / "mrsiprep_p3_csf.mgz"),
    }
    outputs = {}
    for label, (values, path) in specs.items():
        if not path.exists():
            data = np.isin(aseg, list(values)).astype(np.float32)
            header = aseg_img.header.copy()
            header.set_data_dtype(np.float32)
            nib.save(nib.MGHImage(data, aseg_img.affine, header), str(path))
        outputs[label] = path
    return outputs


def _load_probability(path: Path) -> np.ndarray:
    data = np.nan_to_num(nib.load(str(path)).get_fdata(dtype=np.float32).squeeze(), copy=False)
    return np.clip(data.astype(np.float32), 0.0, 1.0)


def copy_freesurfer_subject_for_chimera(config, raw_t1: Path, target_subjects_dir: Path | None = None) -> Path:
    source = config.freesurfer_dir / freesurfer_subject_id(raw_t1)
    if not source.exists():
        raise FileNotFoundError(f"FreeSurfer subject directory not found: {source}")
    if target_subjects_dir is None or Path(target_subjects_dir).resolve() == config.freesurfer_dir.resolve():
        return source
    target = Path(target_subjects_dir) / source.name
    if target.exists():
        return target
    shutil.copytree(source, target)
    return target
