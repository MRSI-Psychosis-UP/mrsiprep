"""Chimera native-space parcellation workflow."""

from __future__ import annotations

from pathlib import Path

from mrsiprep.interfaces.ants import apply_transforms
from mrsiprep.io.bids import BIDSLayout
from mrsiprep.io.naming import parcellation_derivative
from mrsiprep.parcellation.base import ParcellationResult
from mrsiprep.parcellation.labels import copy_labels


def run_chimera_parcellation(config, subject: str, session: str | None, mrsi_reference: Path, t1_to_mrsi_transforms: list[Path]) -> ParcellationResult:
    layout = BIDSLayout(config.bids_dir)
    source_atlas = layout.chimera_atlas(subject, session, config.chimera_scheme, config.chimera_scale, config.chimera_grow, space="orig")
    if source_atlas is None:
        raise FileNotFoundError(
            f"Chimera atlas not found for sub-{subject} ses-{session} scheme={config.chimera_scheme} scale={config.chimera_scale}."
        )
    scale = f"scale{config.chimera_scale}"
    atlas_name = f"chimera{config.chimera_scheme}"
    t1_out = parcellation_derivative(config.derivative_dir, subject, session, space="T1w", atlas=atlas_name, scale=scale)
    mrsi_out = parcellation_derivative(config.derivative_dir, subject, session, space="MRSI", atlas=atlas_name, scale=scale)
    labels_out = parcellation_derivative(config.derivative_dir, subject, session, atlas=atlas_name, scale=scale, suffix_override="tsv")
    if not t1_out.exists() or config.overwrite:
        import shutil

        t1_out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_atlas, t1_out)
    if not mrsi_out.exists() or config.overwrite:
        apply_transforms(mrsi_reference, t1_out, t1_to_mrsi_transforms, mrsi_out, interpolation="genericLabel")
    source_labels = source_atlas.with_suffix("").with_suffix(".tsv") if source_atlas.name.endswith(".nii.gz") else source_atlas.with_suffix(".tsv")
    if source_labels.exists():
        copy_labels(source_labels, labels_out)
    else:
        _labels_from_image(mrsi_out, labels_out)
    return ParcellationResult(atlas_t1=t1_out, atlas_mrsi=mrsi_out, labels=labels_out, mode="chimera", atlas_name=atlas_name, scale=scale)


def _labels_from_image(image_path: Path, labels_path: Path) -> None:
    import nibabel as nib
    import numpy as np

    from mrsiprep.parcellation.labels import write_labels

    data = nib.load(str(image_path)).get_fdata().astype(int)
    indices = np.unique(data)
    indices = indices[indices != 0]
    write_labels(indices, [str(i) for i in indices], labels_path)
