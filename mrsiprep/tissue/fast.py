"""FSL FAST tissue backend."""

from __future__ import annotations

from pathlib import Path

from mrsiprep.interfaces.fsl import run_fast
from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti


def segment_t1_fast(config, subject: str, session: str | None, t1_path: Path) -> dict[str, Path]:
    import nibabel as nib

    tmp_prefix = config.work_dir / f"sub-{subject}" / (f"ses-{session}" if session else "ses-none") / "fast" / "fast"
    fast_maps = run_fast(t1_path, tmp_prefix, verbose=config.verbose)
    out = {}
    for label, path in fast_maps.items():
        target = anat_derivative(config.derivative_dir, subject, session, space="T1w", label=label, suffix_override="probseg")
        img = nib.load(str(path))
        out[label] = save_nifti(img.get_fdata().astype("float32"), img, target, dtype="float32")
    return out
