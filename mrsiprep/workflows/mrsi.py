"""MRSI preprocessing workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mrsiprep.io.loaders import MRSIInputs
from mrsiprep.mrsi.filtering import filter_metabolite_maps
from mrsiprep.mrsi.masks import ensure_brainmask
from mrsiprep.mrsi.quality import make_quality_masks
from mrsiprep.mrsi.reference import generate_reference


@dataclass
class MRSIResult:
    raw_maps: dict[str, Path]
    preproc_maps: dict[str, Path]
    corrected_maps: dict[str, Path]
    crlb_maps: dict[str, Path]
    snr_map: Path | None
    linewidth_map: Path | None
    brainmask: Path
    reference: Path
    qcmasks: dict[str, Path]
    qc_summary: Path


def run_mrsi_workflow(config, subject: str, session: str | None, inputs: MRSIInputs) -> MRSIResult:
    brainmask = ensure_brainmask(config, subject, session, inputs.brainmask, inputs.water_map, inputs.metabolite_maps)
    preproc = filter_metabolite_maps(config, subject, session, inputs.metabolite_maps, brainmask)
    reference = generate_reference(config, subject, session, preproc, preferred_met=config.ref_met)
    qcmasks, qc_summary = make_quality_masks(
        config,
        subject,
        session,
        preproc,
        inputs.crlb_maps,
        inputs.snr_map,
        inputs.linewidth_map,
        brainmask,
    )
    return MRSIResult(
        raw_maps=inputs.metabolite_maps,
        preproc_maps=preproc,
        corrected_maps=preproc,
        crlb_maps=inputs.crlb_maps,
        snr_map=inputs.snr_map,
        linewidth_map=inputs.linewidth_map,
        brainmask=brainmask,
        reference=reference,
        qcmasks=qcmasks,
        qc_summary=qc_summary,
    )
