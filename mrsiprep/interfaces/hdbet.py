"""HD-BET skull-stripping helper."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib

from mrsiprep.utils.misc import parse_bids_entities


def prepare_hd_bet_env() -> dict[str, str]:
    env = os.environ.copy()
    env["MKL_THREADING_LAYER"] = "GNU"
    return env


def find_t1_files(bids_dir: Path, pattern: str = "*T1w.nii*") -> list[Path]:
    return sorted(path for path in bids_dir.rglob(pattern) if path.is_file() and "derivatives" not in path.parts)


def build_outputs(bids_dir: Path, t1_path: Path) -> tuple[Path, Path]:
    meta = parse_bids_entities(t1_path)
    sub = meta.get("sub") or "unknown"
    ses = meta.get("ses") or ""
    acq = meta.get("acq")
    parts = [f"sub-{sub}"]
    if ses:
        parts.append(f"ses-{ses}")
    if meta.get("run"):
        parts.append(f"run-{meta['run']}")
    if acq:
        parts.append(f"acq-{acq}")
    base = "_".join(parts)
    out_dir = bids_dir / "derivatives" / "skullstrip" / f"sub-{sub}" / (f"ses-{ses}" if ses else "")
    return out_dir / f"{base}_desc-brain_T1w.nii.gz", out_dir / f"{base}_desc-brainmask_T1w.nii.gz"


def run_hd_bet(t1_path: Path, out_brain: Path, hd_bet_bin: str = "hd-bet", device: str = "cuda", disable_tta: bool = False, verbose: bool = False) -> None:
    out_brain.parent.mkdir(parents=True, exist_ok=True)
    cmd = [hd_bet_bin, "-i", str(t1_path), "-o", str(out_brain), "-device", device]
    if disable_tta:
        cmd.append("--disable_tta")
    subprocess.run(
        cmd,
        check=True,
        env=prepare_hd_bet_env(),
        stdout=None if verbose else subprocess.PIPE,
        stderr=None if verbose else subprocess.PIPE,
        text=True,
    )


def create_brainmask(out_brain: Path, out_mask: Path, verbose: bool = False) -> None:
    tmp_mask = out_mask.with_name("tmp_" + out_mask.name)
    try:
        subprocess.run(["fslmaths", str(out_brain), "-bin", str(tmp_mask)], check=True, stdout=None if verbose else subprocess.DEVNULL)
        subprocess.run(["fslmaths", str(tmp_mask), "-fillh", str(out_mask)], check=True, stdout=None if verbose else subprocess.DEVNULL)
    finally:
        if tmp_mask.exists():
            tmp_mask.unlink()


def process_t1_file(
    bids_dir: Path,
    t1_path: Path,
    hd_bet_bin: str = "hd-bet",
    device: str = "cuda",
    disable_tta: bool = False,
    verbose: bool = False,
    overwrite: bool = False,
) -> dict:
    out_brain, out_mask = build_outputs(bids_dir, t1_path)
    if out_brain.exists() and out_mask.exists() and not overwrite:
        return {"brain": out_brain, "mask": out_mask, "status": "skipped"}

    hdbet_input = t1_path
    temp_dir = None
    tmp_input = None
    if t1_path.name.endswith(".nii.gz"):
        temp_dir = tempfile.TemporaryDirectory(prefix="mrsiprep_hdbet_")
        tmp_input = Path(temp_dir.name) / t1_path.name[:-7]
        nib.save(nib.load(str(t1_path)), str(tmp_input))
        hdbet_input = tmp_input
    try:
        run_hd_bet(hdbet_input, out_brain, hd_bet_bin, device, disable_tta, verbose)
    finally:
        if tmp_input and tmp_input.exists():
            tmp_input.unlink()
        if temp_dir:
            temp_dir.cleanup()
    create_brainmask(out_brain, out_mask, verbose)
    return {"brain": out_brain, "mask": out_mask, "status": "processed"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Batch HD-BET over BIDS T1w files.")
    parser.add_argument("bids_dir", type=Path)
    parser.add_argument("--pattern", default="*T1w.nii*")
    parser.add_argument("--hd-bet-bin", default="hd-bet")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--disable-tta", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if not args.bids_dir.is_dir():
        raise SystemExit(f"BIDS dir not found: {args.bids_dir}")
    files = find_t1_files(args.bids_dir, args.pattern)
    if not files:
        raise SystemExit("No T1 files found.")
    for path in files:
        result = process_t1_file(args.bids_dir, path, args.hd_bet_bin, args.device, args.disable_tta, args.verbose, args.overwrite)
        print(f"{result['status']}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
