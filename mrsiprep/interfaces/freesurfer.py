"""FreeSurfer helpers."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


class FreeSurferError(RuntimeError):
    """Raised for FreeSurfer-related failures."""


def check_license() -> None:
    import os

    license_path = os.environ.get("FS_LICENSE")
    if not license_path or not Path(license_path).exists():
        raise FreeSurferError("FreeSurfer license not found. Set FS_LICENSE to use FreeSurfer/Chimera mode.")


def subject_dir_valid(fs_subjects_dir: str | Path, subject: str) -> bool:
    root = Path(fs_subjects_dir) / subject
    return (root / "mri").exists() and (root / "surf").exists()


def run_recon_all(t1_path: str | Path, fs_subjects_dir: str | Path, subject: str, force: bool = False, nthreads: int = 4) -> Path:
    if not shutil.which("recon-all"):
        raise FreeSurferError("recon-all command not found on PATH.")
    check_license()
    fs_subjects_dir = Path(fs_subjects_dir)
    if subject_dir_valid(fs_subjects_dir, subject) and not force:
        return fs_subjects_dir / subject
    fs_subjects_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["recon-all", "-s", subject, "-i", str(t1_path), "-all", "-sd", str(fs_subjects_dir), "-openmp", str(nthreads)]
    subprocess.run(cmd, check=True)
    return fs_subjects_dir / subject
