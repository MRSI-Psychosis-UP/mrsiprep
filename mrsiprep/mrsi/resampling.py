"""MRSI resampling helpers."""

from __future__ import annotations

from pathlib import Path

from nilearn import datasets

from mrsiprep.io.naming import mrsi_derivative
from mrsiprep.registration.transforms import apply_image_transform


def transform_mrsi_maps(config, subject: str, session: str | None, maps: dict[str, Path], mrsi_to_t1: list[Path], t1_to_mni: list[Path] | None, t1_reference: Path) -> dict[str, dict[str, Path]]:
    outputs: dict[str, dict[str, Path]] = {}
    if "T1w" in config.output_spaces or "t1w" in config.transform:
        outputs["T1w"] = {}
        for met, path in maps.items():
            out = mrsi_derivative(config.derivative_dir, subject, session, space="T1w", met=met, desc="preproc", suffix_override="mrsi")
            outputs["T1w"][met] = apply_image_transform(t1_reference, path, mrsi_to_t1, out, interpolation="linear") if not out.exists() or config.overwrite_transform else out
    if ("MNI152NLin2009cAsym" in config.output_spaces or "mni" in config.transform) and t1_to_mni:
        outputs["MNI152NLin2009cAsym"] = {}
        template = datasets.load_mni152_template()
        transforms = list(t1_to_mni) + list(mrsi_to_t1)
        for met, path in maps.items():
            out = mrsi_derivative(config.derivative_dir, subject, session, space="MNI152NLin2009cAsym", met=met, desc="preproc", suffix_override="mrsi")
            outputs["MNI152NLin2009cAsym"][met] = apply_image_transform(template, path, transforms, out, interpolation="linear") if not out.exists() or config.overwrite_transform else out
    return outputs
