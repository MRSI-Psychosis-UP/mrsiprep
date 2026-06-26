"""Chimera parcellation wrapper."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path


class ChimeraError(RuntimeError):
    """Raised when Chimera parcellation cannot be created."""


def check_chimera() -> None:
    if not shutil.which("chimera"):
        raise ChimeraError("Chimera command not found on PATH.")


def run_chimera(
    bids_dir: str | Path,
    derivatives_dir: str | Path,
    fs_subjects_dir: str | Path,
    t1_path: str | Path,
    subject: str,
    session: str | None,
    scheme: str,
    scale: int,
    grow: int,
    nthreads: int = 4,
    verbose: bool = False,
) -> Path:
    check_chimera()
    derivatives_dir = Path(derivatives_dir)
    ids_line = f"{Path(t1_path).name}\n"
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
        handle.write(ids_line)
        ids_path = Path(handle.name)
    try:
        cmd = [
            "chimera",
            "-b",
            str(bids_dir),
            "-d",
            str(derivatives_dir),
            "--freesurferdir",
            str(fs_subjects_dir),
            "-p",
            scheme,
            "-g",
            str(grow),
            "-s",
            str(scale),
            "-ids",
            str(ids_path),
            "--nthreads",
            str(nthreads),
        ]
        subprocess.run(cmd, check=True, stdout=None if verbose else subprocess.PIPE, stderr=None if verbose else subprocess.PIPE, text=True)
    finally:
        ids_path.unlink(missing_ok=True)
    pattern = f"sub-{subject}"
    if session:
        pattern += f"_ses-{session}"
    candidates = sorted((derivatives_dir / "chimera-atlases").rglob(f"{pattern}*atlas-chimera{scheme}*scale{scale}grow{grow}mm_dseg.nii*"))
    if not candidates:
        raise ChimeraError("Chimera completed but expected parcellation was not found.")
    return candidates[0]
