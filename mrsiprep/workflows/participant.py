"""Participant workflow orchestration."""

from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from pathlib import Path

from mrsiprep.io.bids import BIDSLayout, Recording
from mrsiprep.io.derivatives import init_derivative
from mrsiprep.io.validators import ValidationError, validate_recording
from mrsiprep.mrsi.pvc import create_tissue_4d, run_pvc
from mrsiprep.mrsi.resampling import transform_mrsi_maps
from mrsiprep.parcellation.extraction import extract_regional_metabolites
from mrsiprep.tissue.ants_atropos import atropos_pve_path, segment_t1_atropos
from mrsiprep.tissue.freesurfer import (
    freesurfer_brain_mask_path,
    freesurfer_brain_path,
    freesurfer_pve_path,
    segment_t1_freesurfer,
)
from mrsiprep.tissue.masks import build_brain_csf_seed_mask
from mrsiprep.utils.logging import LOGGER
from mrsiprep.utils.misc import normalize_session, normalize_subject, read_participant_pairs
from mrsiprep.utils.provenance import write_provenance
from mrsiprep.workflows.anatomical import prepare_anatomical
from mrsiprep.workflows.base import ensure_work_dirs
from mrsiprep.workflows.connectivity import run_connectivity_workflow
from mrsiprep.workflows.mrsi import run_mrsi_workflow
from mrsiprep.workflows.parcellation import run_parcellation_workflow
from mrsiprep.workflows.registration import run_registration_workflow
from mrsiprep.workflows.reports import run_reports_workflow
from mrsiprep.workflows.tissue import run_tissue_workflow


@dataclass
class RecordingStatus:
    subject: str
    session: str | None
    status: str
    outputs: dict = field(default_factory=dict)
    error: str | None = None


def collect_recordings(config) -> list[Recording]:
    if config.participants_file:
        return [Recording(sub, ses) for sub, ses in read_participant_pairs(config.participants_file)]
    subjects = config.participant_label or []
    sessions = config.session_label or []
    if subjects:
        if sessions:
            return [Recording(normalize_subject(sub), normalize_session(ses)) for sub in subjects for ses in sessions]
        return [Recording(normalize_subject(sub), None) for sub in subjects]
    return BIDSLayout(config.bids_dir).discover_recordings()


def run_participant_workflow(config) -> list[RecordingStatus]:
    ensure_work_dirs(config)
    init_derivative(config.derivative_dir)
    statuses: list[RecordingStatus] = []
    for recording in collect_recordings(config):
        try:
            outputs = run_single_recording(config, recording.subject, recording.session)
            statuses.append(RecordingStatus(recording.subject, recording.session, "success", outputs=outputs))
        except Exception as exc:  # batch-safe failure
            msg = f"sub-{recording.subject}" + (f" ses-{recording.session}" if recording.session else "")
            LOGGER.error("%s failed: %s", msg, exc)
            if config.verbose:
                LOGGER.error(traceback.format_exc())
            statuses.append(RecordingStatus(recording.subject, recording.session, "failed", error=str(exc)))
            continue
    return statuses


def run_single_recording(config, subject: str, session: str | None) -> dict:
    subject = normalize_subject(subject)
    session = normalize_session(session)
    LOGGER.info("Processing sub-%s%s", subject, f" ses-{session}" if session else "")

    t1_path, inputs = validate_recording(config, subject, session)
    precomputed_tissue_t1 = None
    p3_override = None
    brain_mask_override = None
    if config.tissue_backend == "ants-atropos":
        layout = BIDSLayout(config.bids_dir)
        raw_t1 = layout.raw_t1(subject, session)
        if raw_t1 is None:
            raise FileNotFoundError(f"Missing raw T1w required for ANTs Atropos segmentation: sub-{subject} ses-{session}")
        brain_mask = layout.brain_mask(subject, session)
        if brain_mask is None and Path(t1_path).resolve() == Path(raw_t1).resolve():
            raise FileNotFoundError(
                f"Missing skull-stripped T1w or brain mask required to seed ANTs Atropos: sub-{subject} ses-{session}"
            )
        atropos_mask = build_brain_csf_seed_mask(config, subject, session, t1_path, raw_t1, brain_mask)
        precomputed_tissue_t1 = segment_t1_atropos(config, subject, session, raw_t1, atropos_mask)
        p3_override = atropos_pve_path(config, subject, session, 3)
    elif config.tissue_backend == "freesurfer":
        layout = BIDSLayout(config.bids_dir)
        raw_t1 = layout.raw_t1(subject, session)
        if raw_t1 is None:
            raise FileNotFoundError(f"Missing raw T1w required for FreeSurfer recon-all: sub-{subject} ses-{session}")
        precomputed_tissue_t1 = segment_t1_freesurfer(config, subject, session, raw_t1, raw_t1_path=raw_t1)
        t1_path = freesurfer_brain_path(config, subject, session)
        brain_mask_override = freesurfer_brain_mask_path(config, subject, session)
        p3_override = freesurfer_pve_path(config, subject, session, 3)

    anat = prepare_anatomical(config, subject, session, t1_path, p3_override=p3_override, brain_mask_override=brain_mask_override)
    mrsi = run_mrsi_workflow(config, subject, session, inputs)
    registration = run_registration_workflow(config, subject, session, mrsi.reference, anat.registration_t1w, anat.registration_mask)
    tissue = run_tissue_workflow(
        config,
        subject,
        session,
        anat.registration_t1w,
        anat.registration_mask,
        mrsi.reference,
        registration.mrsi_to_t1.inverse,
        precomputed_tissue_t1=precomputed_tissue_t1,
    )

    corrected_maps = mrsi.preproc_maps
    tissue_4d = None
    if not config.no_pvc:
        tissue_4d = create_tissue_4d(config, subject, session, tissue.mrsi, mrsi.reference)
        corrected_maps = run_pvc(config, subject, session, mrsi.preproc_maps, tissue_4d, mrsi.brainmask)
        mrsi.corrected_maps = corrected_maps

    transformed = transform_mrsi_maps(
        config,
        subject,
        session,
        corrected_maps,
        registration.mrsi_to_t1.forward,
        registration.t1_to_mni.forward if registration.t1_to_mni else None,
        anat.registration_t1w,
    )

    parcels = run_parcellation_workflow(config, subject, session, mrsi.reference, registration)
    regional = extract_regional_metabolites(
        config,
        subject,
        session,
        corrected_maps,
        parcels,
        mrsi.qcmasks,
        mrsi.snr_map,
        mrsi.linewidth_map,
        mrsi.crlb_maps,
        tissue.mrsi,
    )
    connectivity = run_connectivity_workflow(config, subject, session, regional, parcels)
    outputs = {
        "t1w": anat.t1w,
        "registration_t1w": anat.registration_t1w,
        "mrsi_reference": mrsi.reference,
        "qc_summary": mrsi.qc_summary,
        "tissue_4d": tissue_4d,
        "atlas_mrsi": parcels.atlas_mrsi,
        "regional_table": regional,
        "connectivity": connectivity,
        "transformed_maps": transformed,
    }
    report = run_reports_workflow(config, subject, session, outputs)
    outputs["report"] = report
    outputs["provenance"] = write_provenance(
        config,
        config.derivative_dir / f"sub-{subject}" / (f"ses-{session}" if session else "ses-none") / "mrsiprep_provenance.json",
        {"subject": subject, "session": session, "outputs": outputs},
    )
    return outputs
