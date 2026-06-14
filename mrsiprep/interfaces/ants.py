"""ANTs/antspyx interface."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib


class ANTsError(RuntimeError):
    """Raised when ANTs cannot complete a requested operation."""


def _import_ants():
    try:
        import ants  # type: ignore
    except Exception as exc:
        raise ANTsError("antspyx is required for this operation. Install mrsiprep[ants] or provide existing transforms.") from exc
    return ants


def _load_ants_image(image):
    ants = _import_ants()
    if isinstance(image, (str, Path)):
        path = Path(image)
        if not path.exists():
            raise ANTsError(f"Image path does not exist: {path}")
        return ants.image_read(str(path))
    if isinstance(image, nib.Nifti1Image):
        with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            nib.save(image, str(tmp_path))
            return ants.image_read(str(tmp_path))
        finally:
            tmp_path.unlink(missing_ok=True)
    if image.__class__.__name__ == "ANTsImage":
        return image
    raise ANTsError("Image must be a path, nibabel image, or ANTsImage.")


def register(fixed, moving, out_prefix: str | Path, transform: str = "sr", fixed_mask=None, moving_mask=None, verbose: bool = False) -> dict:
    ants = _import_ants()
    fixed_img = _load_ants_image(fixed)
    moving_img = _load_ants_image(moving)
    fixed_mask_img = _load_ants_image(fixed_mask) if fixed_mask is not None else None
    moving_mask_img = _load_ants_image(moving_mask) if moving_mask is not None else None
    tx = ants.registration(
        fixed=fixed_img,
        moving=moving_img,
        fixed_mask=fixed_mask_img,
        moving_mask=moving_mask_img,
        type_of_transform=f"antsRegistrationSyN[{transform}]",
        verbose=verbose,
    )
    return save_all_transforms(tx, out_prefix)


def save_all_transforms(ants_tx: dict, out_prefix: str | Path) -> dict[str, list[Path]]:
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {"forward": [], "inverse": []}
    for path in ants_tx.get("fwdtransforms", []):
        path = Path(path)
        if "Warp" in path.name:
            out = out_prefix.with_suffix(".syn.nii.gz")
        elif "Affine" in path.name:
            out = out_prefix.with_suffix(".affine.mat")
        else:
            out = out_prefix.parent / f"{out_prefix.name}.{path.name}"
        shutil.copy2(path, out)
        outputs["forward"].append(out)
    for path in ants_tx.get("invtransforms", []):
        path = Path(path)
        if "InverseWarp" in path.name or "Warp" in path.name:
            out = out_prefix.with_suffix(".syn_inv.nii.gz")
        elif "Affine" in path.name:
            out = out_prefix.with_suffix(".affine_inv.mat")
        else:
            out = out_prefix.parent / f"{out_prefix.name}.{path.name}"
        shutil.copy2(path, out)
        outputs["inverse"].append(out)
    return outputs


def apply_transforms(fixed, moving, transforms: list[str | Path], out_path: str | Path | None = None, interpolation: str = "linear"):
    ants = _import_ants()
    existing = [str(Path(path)) for path in transforms if Path(path).exists()]
    if not existing:
        raise ANTsError(f"No transform files exist in list: {transforms}")
    fixed_img = _load_ants_image(fixed)
    moving_img = _load_ants_image(moving)
    warped = ants.apply_transforms(fixed=fixed_img, moving=moving_img, transformlist=existing, interpolator=interpolation)
    if out_path is None:
        return warped
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ants.image_write(warped, str(out_path))
    return out_path


def require_cli(command: str) -> str:
    path = shutil.which(command)
    if not path:
        raise ANTsError(f"Required ANTs command not found on PATH: {command}")
    return path


def run_cli(cmd: list[str], verbose: bool = False) -> None:
    subprocess.run(cmd, check=True, stdout=None if verbose else subprocess.PIPE, stderr=None if verbose else subprocess.PIPE, text=True)
