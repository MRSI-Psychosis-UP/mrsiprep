"""Input validation."""

from __future__ import annotations

from pathlib import Path

from mrsiprep.config.settings import MRSIPrepConfig
from mrsiprep.io.bids import BIDSLayout
from mrsiprep.io.loaders import MRSIInputs, load_mrsi_inputs
from mrsiprep.utils.images import assert_same_grid


class ValidationError(RuntimeError):
    """Raised when a recording cannot be processed."""


def validate_recording(config: MRSIPrepConfig, subject: str, session: str | None) -> tuple[Path, MRSIInputs]:
    layout = BIDSLayout(config.bids_dir)
    t1 = layout.t1(subject, session, config.t1_pattern)
    if not t1 or not t1.exists():
        raise ValidationError(f"Missing T1w reference for sub-{subject} ses-{session}: pattern {config.t1_pattern}")

    inputs = load_mrsi_inputs(layout, subject, session, config.metabolites)
    missing = [met for met in config.metabolites if met not in inputs.metabolite_maps]
    if missing:
        raise ValidationError(f"Missing metabolite maps for sub-{subject} ses-{session}: {', '.join(missing)}")

    if config.registration_t1_target == "brain-csf":
        p3 = layout.cat12_probseg(subject, session, 3)
        if not p3 or not p3.exists():
            raise ValidationError(
                f"Missing CAT12 p3 CSF map required for --registration-t1-target brain-csf: sub-{subject} ses-{session}"
            )

    assert_same_grid(list(inputs.metabolite_maps.values()), "metabolite maps")
    return t1, inputs
