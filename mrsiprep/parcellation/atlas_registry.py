"""MNI atlas loading."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn import datasets, image

from mrsiprep.parcellation.labels import write_labels


def load_mni_atlas(config, work_dir: str | Path) -> tuple[Path, Path, str]:
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    atlas = config.atlas.lower()
    if atlas == "custom":
        if not config.custom_atlas or not config.custom_atlas_lut:
            raise ValueError("--custom-atlas and --custom-atlas-lut are required for atlas=custom.")
        return Path(config.custom_atlas), Path(config.custom_atlas_lut), "custom"
    if atlas.startswith("schaefer"):
        n_rois = int(atlas.replace("schaefer", ""))
        fetched = datasets.fetch_atlas_schaefer_2018(n_rois=n_rois, yeo_networks=7, resolution_mm=1)
        atlas_img = image.resample_to_img(fetched.maps, datasets.load_mni152_template(), interpolation="nearest", force_resample=True)
        atlas_path = work_dir / f"atlas-{atlas}_space-MNI152NLin2009cAsym_dseg.nii.gz"
        nib.save(atlas_img, atlas_path)
        data = atlas_img.get_fdata().astype(int)
        indices = np.unique(data)
        indices = indices[indices != 0]
        labels_path = work_dir / f"atlas-{atlas}_labels.tsv"
        labels = [label.decode() if isinstance(label, bytes) else str(label) for label in fetched.labels]
        write_labels(indices, labels[: len(indices)], labels_path)
        return atlas_path, labels_path, atlas
    if atlas in {"mist197", "mist-197"}:
        fetched = datasets.fetch_atlas_basc_multiscale_2015()
        atlas_img = image.resample_to_img(fetched.scale197, datasets.load_mni152_template(), interpolation="nearest", force_resample=True)
        atlas_path = work_dir / "atlas-mist197_space-MNI152NLin2009cAsym_dseg.nii.gz"
        nib.save(atlas_img, atlas_path)
        indices = np.unique(atlas_img.get_fdata().astype(int))
        indices = indices[indices != 0]
        write_labels(indices, [f"MIST-{i}" for i in indices], work_dir / "atlas-mist197_labels.tsv")
        return atlas_path, work_dir / "atlas-mist197_labels.tsv", "mist197"
    raise ValueError(f"Unsupported MNI atlas: {config.atlas}")
