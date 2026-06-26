"""MRSI resampling helpers."""

from __future__ import annotations

from pathlib import Path

from nilearn import datasets

from mrsiprep.io.naming import mrsi_derivative
from mrsiprep.registration.transforms import apply_image_transform
from mrsiprep.utils.images import resolve_mni_resolution


def transform_mrsi_maps(
    config,
    subject: str,
    session: str | None,
    maps: dict[str, Path],
    mrsi_to_t1: list[Path],
    t1_to_mni: list[Path] | None,
    t1_reference: Path,
    mrsi_reference: Path | None = None,
    crlb_maps: dict[str, Path] | None = None,
    snr_map: Path | None = None,
    linewidth_map: Path | None = None,
) -> dict[str, dict[str, Path]]:
    crlb_maps = crlb_maps or {}
    threads = config.nthreads

    def _spike_mask_path(met: str) -> Path | None:
        path = mrsi_derivative(config.derivative_dir, subject, session, space="MRSI", met=met, desc="spikemask", suffix_override="mask")
        return path if path.exists() else None

    def _quality_items():
        items = [("crlb", met, path) for met, path in crlb_maps.items() if met in maps]
        if snr_map is not None:
            items.append(("snr", None, snr_map))
        if linewidth_map is not None:
            items.append(("fwhm", None, linewidth_map))
        if config.transform_spikemask:
            for met in maps:
                spike_path = _spike_mask_path(met)
                if spike_path is not None:
                    items.append(("spikemask", met, spike_path))
        return items

    def _resample_space(space: str, fixed, transforms: list[Path], res: int | None = None) -> dict[str, Path]:
        space_outputs: dict[str, Path] = {}
        for met, path in maps.items():
            out = mrsi_derivative(config.derivative_dir, subject, session, space=space, res=res, met=met, desc="signal", suffix_override="mrsi")
            space_outputs[met] = (
                apply_image_transform(fixed, path, transforms, out, interpolation="linear", threads=threads)
                if not out.exists() or config.overwrite_transform
                else out
            )
        for desc, met, path in _quality_items():
            interpolation = "genericLabel" if desc == "spikemask" else "linear"
            out = mrsi_derivative(config.derivative_dir, subject, session, space=space, res=res, met=met, desc=desc, suffix_override="mrsi")
            space_outputs[f"{desc}{f'-{met}' if met else ''}"] = (
                apply_image_transform(fixed, path, transforms, out, interpolation=interpolation, threads=threads)
                if not out.exists() or config.overwrite_transform
                else out
            )
        return space_outputs

    outputs: dict[str, dict[str, Path]] = {}
    if "T1w" in config.output_spaces or "t1w" in config.transform:
        outputs["T1w"] = _resample_space("T1w", t1_reference, mrsi_to_t1)
    if ("MNI152NLin2009cAsym" in config.output_spaces or "mni" in config.transform) and t1_to_mni:
        resolution = resolve_mni_resolution(config.mni_resolution, t1_reference, mrsi_reference)
        template = datasets.load_mni152_template(resolution)
        transforms = list(t1_to_mni) + list(mrsi_to_t1)
        outputs["MNI152NLin2009cAsym"] = _resample_space("MNI152NLin2009cAsym", template, transforms, res=resolution)
    return outputs
