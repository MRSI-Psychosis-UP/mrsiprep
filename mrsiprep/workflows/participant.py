"""Participant workflow orchestration."""

from __future__ import annotations

import os
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path

from rich import box
from rich.table import Table

from mrsiprep.io.bids import BIDSLayout, Recording
from mrsiprep.io.derivatives import init_derivative
from mrsiprep.io.loaders import load_mrsi_inputs
from mrsiprep.io.validators import ValidationError, validate_recording
from mrsiprep.mrsi.pvc import create_tissue_4d, run_pvc
from mrsiprep.utils.debug import Debug
from mrsiprep.mrsi.resampling import transform_mrsi_maps
from mrsiprep.parcellation.extraction import extract_regional_metabolites
from mrsiprep.parcellation.metprofiles import export_metprofile_npz
from mrsiprep.parcellation.synthseg import run_synthseg_parcellation
from mrsiprep.reports.parcel_qc import write_parcel_qc
from mrsiprep.reports.qc_combine import combine_qc_reports
from mrsiprep.reports.connectivity_overview import write_connectivity_qc_report
from mrsiprep.reports.mrsi_preproc import write_mrsi_preproc_qc_report
from mrsiprep.reports.parcellation_overview import write_parcellation_qc_report
from mrsiprep.reports.registration_overview import write_registration_overview_report
from mrsiprep.reports.tissue import write_tissue_qc_report
from mrsiprep.tissue.synthseg_fast import (
    extract_t1_synthseg,
    segment_t1_synthseg_fast,
    synthseg_native_labels_path,
    synthseg_fast_brain_path,
    synthseg_fast_brain_mask_path,
    synthseg_fast_csf_probseg_path,
)
from mrsiprep.utils.images import resolve_mni_resolution
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


def _gather_input_availability(config, subject: str, session: str | None) -> dict:
    layout = BIDSLayout(config.bids_dir)
    recording_id = f"sub-{subject}"
    if session:
        recording_id += f"_ses-{session}"

    t1_path = layout.t1(subject, session, config.t1_pattern)
    t1_status = bool(t1_path and t1_path.exists())

    inputs = load_mrsi_inputs(layout, subject, session, config.metabolites)
    total_expected = len(config.metabolites)
    found_count = len(inputs.metabolite_maps)

    brainmask_status = bool(inputs.brainmask)
    if config.tissue_backend == "existing":
        tissue_statuses = [
            bool(layout.cat12_probseg(subject, session, idx))
            for idx in (1, 2, 3)
        ]
        tissue_label = " ".join(
            f"[{'green' if status else 'red'}]p{idx}[/{'green' if status else 'red'}]"
            for idx, status in enumerate(tissue_statuses, 1)
        )
    else:
        tissue_label = "[cyan]AUTO[/cyan]"

    transforms = {}
    for stage in ("mrsi", "anat", "t1-template", "template-mni"):
        stage_paths = layout.transform(subject, session, stage)
        transforms[stage] = bool(stage_paths and all(path.exists() for path in stage_paths))

    return {
        "recording_id": recording_id,
        "subject": subject,
        "session": session,
        "t1": t1_status,
        "mrsi_found": found_count,
        "mrsi_expected": total_expected,
        "brainmask": brainmask_status,
        "tissue": tissue_label,
        "transforms": transforms,
    }


def _render_preflight_table(config, summaries: list[dict], debug: Debug) -> None:
    CHECK_MARK = "[green]✔[/green]"
    CROSS_MARK = "[red]X[/red]"
    PROC_MARK = "[orange3]PROC[/orange3]"
    NA_MARK = "[grey58]N/A[/grey58]"
    transform_columns = [
        ("mrsi", "MRSI→T1"),
        ("anat", "T1→MNI"),
        ("t1-template", "T1→Template"),
        ("template-mni", "Template→MNI"),
    ]

    table = Table(box=box.SIMPLE_HEAVY, show_lines=False, title="Input availability summary")
    table.add_column("Recording", style="cyan", no_wrap=True)
    table.add_column("T1w ref", justify="center")
    table.add_column("MRSI files", justify="center")
    table.add_column("Brainmask", justify="center")
    table.add_column("Tissue files", justify="center")
    for _, label in transform_columns:
        table.add_column(label, justify="center")

    missing_count = 0
    total_missing_files = 0
    missing_recordings = []
    subject_session_counts: dict[str, int] = {}
    for row in summaries:
        subject_session_counts[row["subject"]] = subject_session_counts.get(row["subject"], 0) + (1 if row["session"] else 1)

    for row in summaries:
        mrsi_color = "green" if row["mrsi_found"] == row["mrsi_expected"] else "red"
        mrsi_cell = f"[{mrsi_color}]{row['mrsi_found']}/{row['mrsi_expected']}[/{mrsi_color}]"
        t1_cell = CHECK_MARK if row["t1"] else CROSS_MARK
        brainmask_cell = CHECK_MARK if row["brainmask"] else PROC_MARK

        row_cells = [
            row["recording_id"],
            t1_cell,
            mrsi_cell,
            brainmask_cell,
            row["tissue"],
        ]

        for stage_key, _ in transform_columns:
            if not row["session"] and stage_key in {"t1-template", "template-mni"}:
                row_cells.append(NA_MARK)
            else:
                row_cells.append(CHECK_MARK if row["transforms"].get(stage_key, False) else PROC_MARK)

        table.add_row(*row_cells)

        missing_items = []
        if not row["t1"]:
            missing_items.append("T1w")
        if row["mrsi_found"] != row["mrsi_expected"]:
            missing_items.append(f"MRSI {row['mrsi_found']}/{row['mrsi_expected']}")
        if config.processing_mode == "full" and config.tissue_backend == "existing" and "red" in row["tissue"]:
            missing_items.append("Tissue")

        if missing_items:
            missing_count += 1
            total_missing_files += len(missing_items)
            missing_recordings.append(row["recording_id"])

    debug.separator()
    debug.title("Preflight input availability")
    debug.console.print(table)
    if missing_count:
        debug.error(
            f"Detected {total_missing_files} missing or incomplete file categories across {missing_count}/{len(summaries)} recordings. "
            f"Affected recordings: {', '.join(missing_recordings)}"
        )
    else:
        debug.success("All required inputs are available for the selected recordings.")

    nproc, nthreads, cpu_warning = config.resolve_cpu_budget()
    if cpu_warning:
        debug.always(f"[warning]WARNING:[/warning] {cpu_warning}")
    else:
        debug.always(f"CPU budget: --nproc {nproc} x --nthreads {nthreads} = {nproc * nthreads} threads (of {os.cpu_count()} available).")


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


def _process_one_recording(config, subject: str, session: str | None) -> RecordingStatus:
    debug = Debug(verbose=config.verbose)
    msg = f"sub-{subject}" + (f" ses-{session}" if session else "")
    debug.always(f"[proc]START[/proc] {msg}")
    LOGGER.info("START %s", msg)
    start = time.monotonic()
    try:
        outputs = run_single_recording(config, subject, session)
        elapsed = time.monotonic() - start
        debug.always(f"[success]FINISHED[/success] {msg} in {_format_elapsed(elapsed)}")
        LOGGER.info("FINISHED %s in %s", msg, _format_elapsed(elapsed))
        return RecordingStatus(subject, session, "success", outputs=outputs)
    except Exception as exc:  # batch-safe failure
        elapsed = time.monotonic() - start
        debug.always(f"[failure]FAILED[/failure] {msg} after {_format_elapsed(elapsed)}: {exc}")
        LOGGER.error("FAILED %s after %s: %s", msg, _format_elapsed(elapsed), exc)
        if config.verbose >= 2:
            LOGGER.error(traceback.format_exc())
        return RecordingStatus(subject, session, "failed", error=str(exc))


def _format_elapsed(seconds: float) -> str:
    minutes, secs = divmod(int(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{secs:02d}s"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def run_participant_workflow(config) -> list[RecordingStatus]:
    ensure_work_dirs(config)
    init_derivative(config.derivative_dir)
    debug = Debug(verbose=config.verbose)

    recordings = collect_recordings(config)
    ready: list[Recording] = []
    statuses: list[RecordingStatus] = []
    for recording in recordings:
        subject = normalize_subject(recording.subject)
        session = normalize_session(recording.session)
        try:
            validate_recording(config, subject, session)
            _validate_backend_inputs(config, subject, session)
            ready.append(Recording(subject, session))
        except (ValidationError, FileNotFoundError) as exc:
            msg = f"sub-{subject}" + (f" ses-{session}" if session else "")
            debug.error("SKIP", msg, str(exc))
            statuses.append(RecordingStatus(subject, session, "skipped", error=str(exc)))

    if not ready:
        return statuses

    if config.nproc <= 1:
        for recording in ready:
            statuses.append(_process_one_recording(config, recording.subject, recording.session))
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(max_workers=config.nproc) as executor:
            futures = {
                executor.submit(_process_one_recording, config, recording.subject, recording.session): recording
                for recording in ready
            }
            for future in as_completed(futures):
                statuses.append(future.result())
    return statuses


def validate_participant_inputs(config) -> list[RecordingStatus]:
    debug = Debug(verbose=config.verbose)
    recordings = collect_recordings(config)
    summaries = [
        _gather_input_availability(config, normalize_subject(rec.subject), normalize_session(rec.session))
        for rec in recordings
    ]
    _render_preflight_table(config, summaries, debug)

    statuses: list[RecordingStatus] = []
    for recording in recordings:
        subject = normalize_subject(recording.subject)
        session = normalize_session(recording.session)
        try:
            t1_path, inputs = validate_recording(config, subject, session)
            _validate_backend_inputs(config, subject, session)
            outputs = {
                "t1w": t1_path,
                "metabolites": sorted(inputs.metabolite_maps),
                "snr": inputs.snr_map,
                "linewidth": inputs.linewidth_map,
                "brainmask": inputs.brainmask,
            }
            statuses.append(RecordingStatus(subject, session, "success", outputs=outputs))
            debug.success("VALID", f"sub-{subject}", f"ses-{session}" if session else "")
        except Exception as exc:
            statuses.append(RecordingStatus(subject, session, "failed", error=str(exc)))
            debug.error("INVALID", f"sub-{subject}", f"ses-{session}" if session else "", str(exc))
    return statuses


def _validate_backend_inputs(config, subject: str, session: str | None) -> None:
    layout = BIDSLayout(config.bids_dir)
    raw_t1 = layout.raw_t1(subject, session)
    if config.processing_mode == "light" and raw_t1 is None:
        raise FileNotFoundError(f"Missing raw T1w required for light-mode SynthSeg parcellation: sub-{subject} ses-{session}")
    if config.tissue_backend == "synthseg-fast" and raw_t1 is None:
        raise FileNotFoundError(f"Missing raw T1w required for {config.tissue_backend}: sub-{subject} ses-{session}")
    if config.parcellation_mode == "mni" and config.atlas == "custom":
        if not config.custom_atlas or not config.custom_atlas.exists():
            raise FileNotFoundError("--custom-atlas is required for --parcellation-mode mni --atlas custom")
        if not config.custom_atlas_lut or not config.custom_atlas_lut.exists():
            raise FileNotFoundError("--custom-atlas-lut is required for --parcellation-mode mni --atlas custom")


def run_single_recording(config, subject: str, session: str | None) -> dict:
    subject = normalize_subject(subject)
    session = normalize_session(session)
    LOGGER.info("Processing sub-%s%s", subject, f" ses-{session}" if session else "")
    debug = Debug(verbose=config.verbose)

    t1_path, inputs = validate_recording(config, subject, session)
    layout = BIDSLayout(config.bids_dir)
    raw_t1 = layout.raw_t1(subject, session)
    if raw_t1 is None:
        raise FileNotFoundError(f"Missing raw T1w required for MRSIPrep: sub-{subject} ses-{session}")
    precomputed_tissue_t1 = None
    p3_override = None
    brain_mask_override = None
    with debug.step("Tissue segmentation"):
        if config.processing_mode == "light":
            synthseg_brain, synthseg_mask = extract_t1_synthseg(config, subject, session, raw_t1)
            if config.registration_t1_target == "brain":
                t1_path = synthseg_brain
                brain_mask_override = synthseg_mask
        elif config.processing_mode == "full" and config.tissue_backend == "synthseg-fast":
            if raw_t1 is None:
                raise FileNotFoundError(f"Missing raw T1w required for SynthSeg+FAST segmentation: sub-{subject} ses-{session}")
            precomputed_tissue_t1 = segment_t1_synthseg_fast(config, subject, session, raw_t1)
            t1_path = synthseg_fast_brain_path(config, subject, session)
            brain_mask_override = synthseg_fast_brain_mask_path(config, subject, session)
            p3_override = synthseg_fast_csf_probseg_path(config, subject, session)

    with debug.step("Anatomical preparation"):
        anat = prepare_anatomical(config, subject, session, t1_path, p3_override=p3_override, brain_mask_override=brain_mask_override)

    with debug.step("MRSI preprocessing"):
        mrsi = run_mrsi_workflow(config, subject, session, inputs)
        qc_report_mrsi_preproc = write_mrsi_preproc_qc_report(config, subject, session, mrsi.raw_maps, mrsi.preproc_maps)

    with debug.step("MRSI-T1w-MNI registration"):
        registration = run_registration_workflow(config, subject, session, mrsi.reference, anat.registration_t1w, anat.registration_mask)

    tissue = None
    if config.processing_mode == "full":
        with debug.step("Tissue probability maps in MRSI space"):
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

    dseg_for_qc = synthseg_native_labels_path(config, subject, session) if config.processing_mode == "full" and config.tissue_backend == "synthseg-fast" else None
    qc_report_tissue = write_tissue_qc_report(config, subject, session, raw_t1, dseg_for_qc, tissue.t1 if tissue is not None else None)

    corrected_maps = mrsi.preproc_maps
    tissue_4d = None
    if config.processing_mode == "full" and not config.no_pvc:
        assert tissue is not None
        with debug.step("Partial volume correction"):
            tissue_4d = create_tissue_4d(config, subject, session, tissue.mrsi, mrsi.reference)
            corrected_maps = run_pvc(config, subject, session, mrsi.preproc_maps, tissue_4d, mrsi.brainmask)
            mrsi.corrected_maps = corrected_maps

    with debug.step("Resampling MRSI maps to T1w/MNI space"):
        transformed = transform_mrsi_maps(
            config,
            subject,
            session,
            corrected_maps,
            registration.mrsi_to_t1.forward,
            registration.t1_to_mni.forward if registration.t1_to_mni else None,
            anat.registration_t1w,
            mrsi_reference=mrsi.reference,
            crlb_maps=mrsi.crlb_maps,
            snr_map=mrsi.snr_map,
            linewidth_map=mrsi.linewidth_map,
        )
        mni_resolution = resolve_mni_resolution(config.mni_resolution, anat.registration_t1w, mrsi.reference) if registration.t1_to_mni else None
        qc_report_registration = write_registration_overview_report(
            config,
            subject,
            session,
            raw_t1,
            transformed.get("T1w", {}).get(config.ref_met),
            transformed.get("MNI152NLin2009cAsym", {}).get(config.ref_met),
            mni_resolution=mni_resolution,
            orig_ref_map_path=corrected_maps.get(config.ref_met),
            mrsi_to_t1_transforms=registration.mrsi_to_t1.forward,
        )

    with debug.step("SynthSeg parcellation and QC"):
        preliminary_parcels = run_synthseg_parcellation(
            config,
            subject,
            session,
            raw_t1,
            mrsi.reference,
            registration.mrsi_to_t1.inverse,
        )
        parcel_qc = write_parcel_qc(
            config,
            subject,
            session,
            preliminary_parcels,
            raw_t1,
            mrsi.brainmask,
            registration.mrsi_to_t1.forward,
            mrsi.crlb_maps,
            mrsi.qcmasks,
        )
    parcels = preliminary_parcels
    qc_report_parcellation = None
    if config.processing_mode == "full":
        with debug.step("Parcellation"):
            parcels = run_parcellation_workflow(
                config,
                subject,
                session,
                mrsi.reference,
                registration,
                raw_t1=raw_t1,
                t1_reference=anat.registration_t1w,
            )
            qc_report_parcellation = write_parcellation_qc_report(config, subject, session, raw_t1, parcels.atlas_t1, parcels.labels)

    with debug.step("Regional metabolite extraction"):
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
            tissue.mrsi if tissue is not None else {},
        )

    with debug.step("Connectivity", live=False):
        connectivity = run_connectivity_workflow(
            config,
            subject,
            session,
            regional,
            parcels,
            corrected_maps,
            mrsi.crlb_maps,
            mrsi.brainmask,
            gm_fraction_path=tissue.mrsi.get("GM") if tissue is not None else None,
        )
        qc_report_connectivity = write_connectivity_qc_report(config, subject, session, connectivity.get("matrix_tsv"))

    metprofiles = None
    if config.processing_mode == "full":
        metprofiles = export_metprofile_npz(
            config,
            subject,
            session,
            corrected_maps,
            mrsi.water_map,
            parcels,
            regional,
            anat.registration_mask,
        )
    outputs = {
        "t1w": anat.t1w,
        "registration_t1w": anat.registration_t1w,
        "mrsi_reference": mrsi.reference,
        "qc_summary": mrsi.qc_summary,
        "parcel_qc": parcel_qc,
        "tissue_4d": tissue_4d,
        "atlas_mrsi": parcels.atlas_mrsi,
        "preliminary_atlas_mrsi": preliminary_parcels.atlas_mrsi,
        "regional_table": regional,
        "metprofiles": metprofiles,
        "connectivity": connectivity,
        "transformed_maps": transformed,
        "qc_report_tissue": qc_report_tissue,
        "qc_report_mrsi_preproc": qc_report_mrsi_preproc,
        "qc_report_registration": qc_report_registration,
        "qc_report_parcellation": qc_report_parcellation,
        "qc_report_connectivity": qc_report_connectivity,
    }
    with debug.step("Reports"):
        report = run_reports_workflow(config, subject, session, outputs)
        outputs["report"] = report
        outputs["qc_report_combined"] = combine_qc_reports(config, subject, session)
        for key in ("qc_report_tissue", "qc_report_mrsi_preproc", "qc_report_registration", "qc_report_parcellation", "qc_report_connectivity"):
            outputs.pop(key, None)
    outputs["provenance"] = write_provenance(
        config,
        config.derivative_dir / f"sub-{subject}" / (f"ses-{session}" if session else "ses-none") / "mrsiprep_provenance.json",
        {"subject": subject, "session": session, "outputs": outputs},
    )
    return outputs
