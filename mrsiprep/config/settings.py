"""Runtime configuration objects."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

from .defaults import METABOLITES_3T, QUALITY_DEFAULTS


@dataclass
class MRSIPrepConfig:
    bids_dir: Path
    output_dir: Path
    analysis_level: str
    participant_label: list[str] = field(default_factory=list)
    session_label: list[str] = field(default_factory=list)
    participants_file: Path | None = None
    b0: float = 3.0
    metabolites: list[str] = field(default_factory=lambda: list(METABOLITES_3T))
    quality_metrics: list[str] = field(default_factory=lambda: ["snr", "linewidth", "crlb"])
    snr_min: float = QUALITY_DEFAULTS["snr_min"]
    linewidth_max: float = QUALITY_DEFAULTS["linewidth_max"]
    crlb_max: float = QUALITY_DEFAULTS["crlb_max"]
    tissue_backend: str = "freesurfer"
    registration_backend: str = "ants"
    normalization: str = "simple"
    output_spaces: list[str] = field(default_factory=lambda: ["T1w", "MNI152NLin2009cAsym"])
    registration_t1_target: str = "brain-csf"
    csf_pv_threshold: float = 0.95
    atropos_mask_dilation_mm: float = 4.0
    parcellation_mode: str = "chimera"
    chimera_scheme: str = "LFMIHIFIS"
    chimera_scale: int = 3
    chimera_grow: int = 2
    atlas: str = "schaefer200"
    custom_atlas: Path | None = None
    custom_atlas_lut: Path | None = None
    fs_subjects_dir: Path | None = None
    write_connectivity: bool = False
    connectivity_method: str = "spearman"
    connectivity_space: str = "MRSI"
    regional_summary: str = "mean"
    extraction_mode: str = "hard"
    nthreads: int = 4
    ref_met: str = "CrPCr"
    t1_pattern: str = "desc-brain_T1w"
    transform: str = "mni-origres"
    filter_biharmonic: bool = True
    spike_percentile: float = 99.0
    no_pvc: bool = False
    proc_mnilong: bool = False
    overwrite: bool = False
    overwrite_filt: bool = False
    overwrite_pve: bool = False
    overwrite_t1_reg: bool = False
    overwrite_mni_reg: bool = False
    overwrite_transform: bool = False
    overwrite_freesurfer: bool = False
    work_dir: Path | None = None
    verbose: bool = False

    def __post_init__(self) -> None:
        self.bids_dir = Path(self.bids_dir).resolve()
        self.output_dir = Path(self.output_dir).resolve()
        if self.work_dir is None:
            self.work_dir = self.output_dir / "work"
        else:
            self.work_dir = Path(self.work_dir).resolve()
        if self.fs_subjects_dir is not None:
            self.fs_subjects_dir = Path(self.fs_subjects_dir).resolve()

    @property
    def derivative_dir(self) -> Path:
        return self.output_dir / "mrsiprep" if self.output_dir.name == "derivatives" else self.output_dir

    @property
    def source_derivatives_dir(self) -> Path:
        return self.bids_dir / "derivatives"

    @property
    def freesurfer_dir(self) -> Path:
        if self.fs_subjects_dir is not None:
            return self.fs_subjects_dir
        return self.output_dir / "freesurfer"

    def to_dict(self) -> dict:
        out = asdict(self)
        for key, value in list(out.items()):
            if isinstance(value, Path):
                out[key] = str(value)
            elif isinstance(value, list):
                out[key] = [str(item) if isinstance(item, Path) else item for item in value]
        return out
