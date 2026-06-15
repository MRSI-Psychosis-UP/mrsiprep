"""CLI parser for MRSIPrep."""

from __future__ import annotations

import argparse
from pathlib import Path

from mrsiprep.config.defaults import METABOLITES_3T, METABOLITES_7T, QUALITY_DEFAULTS
from mrsiprep.config.settings import MRSIPrepConfig


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mrsiprep", description="Preprocess quantified whole-brain MRSI derivatives.")
    parser.add_argument("bids_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("analysis_level", choices=["participant"])
    parser.add_argument("--participant-label", nargs="+", default=[])
    parser.add_argument("--session-label", nargs="+", default=[])
    parser.add_argument("--participants", type=Path, default=None, help="TSV/CSV subject-session list.")
    parser.add_argument("--b0", type=float, default=3.0, choices=[3.0, 7.0])
    parser.add_argument("--metabolites", nargs="+", default=None)
    parser.add_argument("--quality-metrics", nargs="+", default=["snr", "linewidth", "crlb"])
    parser.add_argument("--snr-min", type=float, default=QUALITY_DEFAULTS["snr_min"])
    parser.add_argument("--linewidth-max", type=float, default=QUALITY_DEFAULTS["linewidth_max"])
    parser.add_argument("--crlb-max", type=float, default=QUALITY_DEFAULTS["crlb_max"])

    parser.add_argument("--tissue-backend", choices=["freesurfer", "existing", "ants-atropos", "fast"], default="freesurfer")
    parser.add_argument("--registration-backend", choices=["ants"], default="ants")
    parser.add_argument("--normalization", choices=["simple", "ants-syn", "existing"], default="simple")
    parser.add_argument("--output-spaces", nargs="+", default=["T1w", "MNI152NLin2009cAsym"])
    parser.add_argument("--registration-t1-target", choices=["brain-csf", "brain", "raw"], default="brain-csf")
    parser.add_argument("--csf-pv-threshold", type=float, default=0.95)
    parser.add_argument("--atropos-mask-dilation-mm", type=float, default=4.0)
    parser.add_argument("--ref-met", default="CrPCr")
    parser.add_argument("--t1", dest="t1_pattern", default="desc-brain_T1w")

    parser.add_argument("--parcellation-mode", choices=["chimera", "mni"], default="chimera")
    parser.add_argument("--chimera-scheme", default="LFMIHIFIS")
    parser.add_argument("--chimera-scale", type=_parse_scale, default=3)
    parser.add_argument("--chimera-grow", type=int, default=2)
    parser.add_argument("--atlas", default="schaefer200")
    parser.add_argument("--custom-atlas", type=Path, default=None)
    parser.add_argument("--custom-atlas-lut", type=Path, default=None)
    parser.add_argument("--fs-subjects-dir", type=Path, default=None)
    parser.add_argument("--extraction-mode", choices=["hard", "soft"], default="hard")

    parser.add_argument("--write-connectivity", action="store_true")
    parser.add_argument("--connectivity-method", choices=["pearson", "spearman", "cosine", "euclidean_distance"], default="spearman")
    parser.add_argument("--connectivity-space", choices=["MRSI", "T1w", "MNI"], default="MRSI")
    parser.add_argument("--regional-summary", choices=["mean", "median", "weighted_mean"], default="mean")

    parser.add_argument("--transform", default="mni-origres")
    parser.add_argument("--no-filter", action="store_true")
    parser.add_argument("--spikepc", type=float, default=99.0)
    parser.add_argument("--no-pvc", action="store_true")
    parser.add_argument("--proc-mnilong", action="store_true")
    parser.add_argument("--nthreads", type=int, default=4)
    parser.add_argument("--work-dir", type=Path, default=None)

    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--overwrite-filt", action="store_true")
    parser.add_argument("--overwrite-pve", action="store_true")
    parser.add_argument("--overwrite-t1-reg", action="store_true")
    parser.add_argument("--overwrite-mni-reg", action="store_true")
    parser.add_argument("--overwrite-transform", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> MRSIPrepConfig:
    args = build_parser().parse_args(argv)
    metabolites = args.metabolites
    if metabolites is None:
        metabolites = list(METABOLITES_7T if args.b0 == 7.0 else METABOLITES_3T)
    return MRSIPrepConfig(
        bids_dir=args.bids_dir,
        output_dir=args.output_dir,
        analysis_level=args.analysis_level,
        participant_label=args.participant_label,
        session_label=args.session_label,
        participants_file=args.participants,
        b0=args.b0,
        metabolites=metabolites,
        quality_metrics=args.quality_metrics,
        snr_min=args.snr_min,
        linewidth_max=args.linewidth_max,
        crlb_max=args.crlb_max,
        tissue_backend=args.tissue_backend,
        registration_backend=args.registration_backend,
        normalization=args.normalization,
        output_spaces=args.output_spaces,
        registration_t1_target=args.registration_t1_target,
        csf_pv_threshold=args.csf_pv_threshold,
        atropos_mask_dilation_mm=args.atropos_mask_dilation_mm,
        ref_met=args.ref_met,
        t1_pattern=args.t1_pattern,
        parcellation_mode=args.parcellation_mode,
        chimera_scheme=args.chimera_scheme,
        chimera_scale=args.chimera_scale,
        chimera_grow=args.chimera_grow,
        atlas=args.atlas,
        custom_atlas=args.custom_atlas,
        custom_atlas_lut=args.custom_atlas_lut,
        fs_subjects_dir=args.fs_subjects_dir,
        extraction_mode=args.extraction_mode,
        write_connectivity=args.write_connectivity,
        connectivity_method=args.connectivity_method,
        connectivity_space=args.connectivity_space,
        regional_summary=args.regional_summary,
        transform=args.transform,
        filter_biharmonic=not args.no_filter,
        spike_percentile=args.spikepc,
        no_pvc=args.no_pvc,
        proc_mnilong=args.proc_mnilong,
        nthreads=args.nthreads,
        work_dir=args.work_dir,
        overwrite=args.overwrite,
        overwrite_filt=args.overwrite_filt,
        overwrite_pve=args.overwrite_pve,
        overwrite_t1_reg=args.overwrite_t1_reg,
        overwrite_mni_reg=args.overwrite_mni_reg,
        overwrite_transform=args.overwrite_transform,
        verbose=args.verbose,
    )


def _parse_scale(value) -> int:
    text = str(value)
    if text.startswith("scale"):
        text = text[len("scale") :]
    return int(text)
