"""ANTs Atropos tissue segmentation."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.interfaces.ants import ANTsError, _import_ants
from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti


def atropos_pve_path(config, subject: str, session: str | None, index: int) -> Path:
    return anat_derivative(config.derivative_dir, subject, session, desc=f"p{index}")


def segment_t1_atropos(config, subject: str, session: str | None, t1_path: Path, brain_mask: Path | None = None) -> dict[str, Path]:
    ants = _import_ants()
    outputs = {
        "GM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="GM", suffix_override="probseg"),
        "WM": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="WM", suffix_override="probseg"),
        "CSF": anat_derivative(config.derivative_dir, subject, session, space="T1w", label="CSF", suffix_override="probseg"),
    }
    pve_outputs = {
        "GM": atropos_pve_path(config, subject, session, 1),
        "WM": atropos_pve_path(config, subject, session, 2),
        "CSF": atropos_pve_path(config, subject, session, 3),
    }
    if all(path.exists() for path in [*outputs.values(), *pve_outputs.values()]) and not config.overwrite:
        return outputs

    image = ants.image_read(str(t1_path))
    mask = ants.image_read(str(brain_mask)) if brain_mask and Path(brain_mask).exists() else ants.get_mask(image)
    corrected = ants.n4_bias_field_correction(image, mask=mask)
    result = ants.atropos(a=corrected, x=mask, i="KMeans[3]", m="[0.2,1x1x1]", c="[5,0]", priorweight=0.0)
    seg = result["segmentation"].numpy()
    probs = result["probabilityimages"]
    corrected_np = corrected.numpy()
    class_means = []
    for idx in range(1, 4):
        vals = corrected_np[seg == idx]
        class_means.append((float(np.nanmean(vals)) if vals.size else np.inf, idx))
    class_means.sort(key=lambda item: item[0])
    # T1 intensities are usually CSF < GM < WM in a skull-stripped image.
    label_for_idx = {
        class_means[0][1]: "CSF",
        class_means[1][1]: "GM",
        class_means[2][1]: "WM",
    }
    ref = nib.load(str(t1_path))
    for idx, prob_img in enumerate(probs, start=1):
        label = label_for_idx.get(idx)
        if label is None:
            continue
        data = prob_img.numpy().astype(np.float32)
        save_nifti(data, ref, outputs[label], dtype=np.float32)
        save_nifti(data, ref, pve_outputs[label], dtype=np.float32)
    if not all(path.exists() for path in [*outputs.values(), *pve_outputs.values()]):
        raise ANTsError("Atropos did not produce all GM/WM/CSF and p1/p2/p3 probability maps.")
    return outputs
