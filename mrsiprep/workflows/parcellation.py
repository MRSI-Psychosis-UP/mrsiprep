"""Parcellation workflow."""

from __future__ import annotations

from mrsiprep.parcellation.chimera_native import run_chimera_parcellation
from mrsiprep.parcellation.mni_atlas import run_mni_parcellation


def run_parcellation_workflow(config, subject, session, mrsi_reference, registration_result):
    if config.parcellation_mode == "chimera":
        return run_chimera_parcellation(config, subject, session, mrsi_reference, registration_result.mrsi_to_t1.inverse)
    if config.parcellation_mode == "mni":
        if registration_result.t1_to_mni is None:
            raise RuntimeError("MNI parcellation requires T1-to-MNI normalization.")
        return run_mni_parcellation(
            config,
            subject,
            session,
            mrsi_reference,
            registration_result.t1_to_mni.inverse,
            registration_result.mrsi_to_t1.inverse,
        )
    raise ValueError(f"Unsupported parcellation mode: {config.parcellation_mode}")
