"""FreeSurfer helpers."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


class FreeSurferError(RuntimeError):
    """Raised for FreeSurfer-related failures."""


def check_license() -> None:
    import os

    candidates = []
    if os.environ.get("FS_LICENSE"):
        candidates.append(Path(os.environ["FS_LICENSE"]))
    if os.environ.get("FREESURFER_HOME"):
        candidates.append(Path(os.environ["FREESURFER_HOME"]) / "license.txt")
    candidates.append(Path("/opt/freesurfer/license.txt"))
    for license_path in candidates:
        if license_path.exists():
            os.environ["FS_LICENSE"] = str(license_path)
            return
    raise FreeSurferError("FreeSurfer license not found. Set FS_LICENSE to use FreeSurfer/Chimera mode.")


def require_command(command: str) -> str:
    path = shutil.which(command)
    if not path:
        raise FreeSurferError(f"{command} command not found on PATH.")
    return path


def freesurfer_subject_id(t1_path: str | Path) -> str:
    path = Path(t1_path)
    name = path.name
    if name.endswith(".nii.gz"):
        return name[:-7]
    if name.endswith(".mgz"):
        return name[:-4]
    return path.stem


def subject_dir_valid(fs_subjects_dir: str | Path, subject: str) -> bool:
    root = Path(fs_subjects_dir) / subject
    required = [
        root / "mri" / "brain.mgz",
        root / "mri" / "aseg.mgz",
        root / "mri" / "orig.mgz",
        root / "surf" / "lh.white",
        root / "surf" / "rh.white",
        root / "surf" / "lh.pial",
        root / "surf" / "rh.pial",
    ]
    return all(path.exists() for path in required)


def run_recon_all(t1_path: str | Path, fs_subjects_dir: str | Path, subject: str, force: bool = False, nthreads: int = 4, verbose: bool = False) -> Path:
    require_command("recon-all")
    check_license()
    fs_subjects_dir = Path(fs_subjects_dir)
    subject_root = fs_subjects_dir / subject
    if force and subject_root.exists():
        shutil.rmtree(subject_root)
    if subject_dir_valid(fs_subjects_dir, subject) and not force:
        return subject_root
    fs_subjects_dir.mkdir(parents=True, exist_ok=True)
    if subject_root.exists() and not force:
        cmd = ["recon-all", "-s", subject, "-all", "-sd", str(fs_subjects_dir), "-openmp", str(nthreads)]
    else:
        cmd = ["recon-all", "-s", subject, "-i", str(t1_path), "-all", "-sd", str(fs_subjects_dir), "-openmp", str(nthreads)]
    subprocess.run(cmd, check=True, stdout=None if verbose else subprocess.PIPE, stderr=None if verbose else subprocess.PIPE, text=True)
    if not subject_dir_valid(fs_subjects_dir, subject):
        raise FreeSurferError(
            f"recon-all finished but required outputs are missing for {subject}: "
            f"{fs_subjects_dir / subject}"
        )
    return fs_subjects_dir / subject
