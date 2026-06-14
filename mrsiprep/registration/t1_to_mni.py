"""T1-to-MNI registration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from nilearn import datasets

from mrsiprep.interfaces.ants import register
from mrsiprep.registration.transforms import all_exist, ants_transform_prefix, transform_paths
from mrsiprep.utils.images import mean_resolution


@dataclass
class T1ToMNIResult:
    forward: list[Path]
    inverse: list[Path]
    prefix: Path
    template: object


def run_t1_to_mni(config, subject: str, session: str | None, t1_path: Path) -> T1ToMNIResult:
    prefix = ants_transform_prefix(config.derivative_dir / "transforms" / "ants", subject, session, "anat")
    forward = transform_paths(prefix, "forward")
    inverse = transform_paths(prefix, "inverse")
    resolution = mean_resolution(t1_path)
    template = datasets.load_mni152_template(resolution)
    if all_exist(forward) and all_exist(inverse) and not (config.overwrite_mni_reg or config.overwrite):
        return T1ToMNIResult(forward, inverse, prefix, template)
    transform = "s" if config.normalization in {"simple", "ants-syn"} else "r"
    register(template, t1_path, prefix, transform=transform, verbose=config.verbose)
    return T1ToMNIResult(transform_paths(prefix, "forward"), transform_paths(prefix, "inverse"), prefix, template)
