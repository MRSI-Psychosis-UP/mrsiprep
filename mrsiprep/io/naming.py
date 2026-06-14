"""BIDS-like naming helpers."""

from __future__ import annotations

from pathlib import Path

from mrsiprep.utils.misc import normalize_session, normalize_subject


def prefix(subject: str, session: str | None = None) -> str:
    parts = [f"sub-{normalize_subject(subject)}"]
    ses = normalize_session(session)
    if ses:
        parts.append(f"ses-{ses}")
    return "_".join(parts)


def subject_session_dir(root: Path, subject: str, session: str | None, suffix: str) -> Path:
    sub_dir = root / f"sub-{normalize_subject(subject)}"
    ses = normalize_session(session)
    return sub_dir / f"ses-{ses}" / suffix if ses else sub_dir / suffix


def anat_derivative(root: Path, subject: str, session: str | None, **entities) -> Path:
    return _derivative(root, subject, session, "anat", "T1w", **entities)


def mrsi_derivative(root: Path, subject: str, session: str | None, **entities) -> Path:
    return _derivative(root, subject, session, "mrsi", "mrsi", **entities)


def parcellation_derivative(root: Path, subject: str, session: str | None, **entities) -> Path:
    return _derivative(root, subject, session, "parcellations", "dseg", **entities)


def connectome_derivative(root: Path, subject: str, session: str | None, suffix: str, **entities) -> Path:
    return _derivative(root, subject, session, "connectomics", suffix, **entities)


def figure_derivative(root: Path, subject: str, session: str | None, extension: str = "svg", **entities) -> Path:
    return _derivative(root, subject, session, "figures", extension, **entities)


def _derivative(root: Path, subject: str, session: str | None, folder: str, suffix: str, **entities) -> Path:
    out_dir = subject_session_dir(root, subject, session, folder)
    suffix = entities.pop("suffix_override", suffix)
    name_parts = [prefix(subject, session)]
    order = [
        "space",
        "res",
        "met",
        "label",
        "atlas",
        "scale",
        "from_",
        "to",
        "mode",
        "desc",
    ]
    for key in order:
        if key not in entities:
            continue
        value = entities[key]
        if value is None:
            continue
        entity = "from" if key == "from_" else key
        name_parts.append(f"{entity}-{value}")
    ext = suffix if suffix.startswith(".") else f"_{suffix}.nii.gz" if suffix in {"T1w", "mrsi", "dseg", "mask", "probseg"} else f".{suffix}"
    if ext.startswith("_"):
        filename = "_".join(name_parts) + ext
    else:
        filename = "_".join(name_parts) + ext
    return out_dir / filename


def transform_path(root: Path, subject: str, session: str | None, from_space: str, to_space: str, desc: str | None = None, ext: str = ".h5") -> Path:
    out_dir = subject_session_dir(root, subject, session, "anat" if from_space == "T1w" else "mrsi")
    parts = [prefix(subject, session), f"from-{from_space}", f"to-{to_space}"]
    if desc:
        parts.append(f"desc-{desc}")
    parts.append("xfm")
    return out_dir / ("_".join(parts) + ext)
