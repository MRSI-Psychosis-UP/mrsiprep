"""Tissue workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mrsiprep.tissue.ants_atropos import segment_t1_atropos
from mrsiprep.tissue.fast import segment_t1_fast
from mrsiprep.tissue.fractions import copy_tissue_to_derivatives, load_existing_cat12, resample_tissue_to_mrsi


@dataclass
class TissueResult:
    t1: dict[str, Path]
    mrsi: dict[str, Path]


def run_tissue_workflow(config, subject: str, session: str | None, t1_path: Path, brain_mask: Path | None, mrsi_reference: Path, t1_to_mrsi_transforms: list[Path]) -> TissueResult:
    backend = config.tissue_backend
    if backend == "existing":
        tissue_t1 = copy_tissue_to_derivatives(config, subject, session, load_existing_cat12(config, subject, session))
    elif backend == "ants-atropos":
        tissue_t1 = segment_t1_atropos(config, subject, session, t1_path, brain_mask)
    elif backend == "fast":
        tissue_t1 = segment_t1_fast(config, subject, session, t1_path)
    elif backend == "freesurfer":
        from mrsiprep.tissue.freesurfer import segment_t1_freesurfer

        tissue_t1 = segment_t1_freesurfer(config, subject, session, t1_path)
    else:
        raise ValueError(f"Unsupported tissue backend: {backend}")
    tissue_mrsi = resample_tissue_to_mrsi(config, subject, session, tissue_t1, mrsi_reference, t1_to_mrsi_transforms)
    return TissueResult(t1=tissue_t1, mrsi=tissue_mrsi)
