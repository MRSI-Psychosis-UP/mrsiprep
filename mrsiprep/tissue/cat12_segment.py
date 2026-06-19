"""CAT12-like tissue segmentation scaffold.

This module is the high-level Python entry point for the staged CAT12 port:

1. build initial native-space GM/WM/CSF maps with Python-native tools;
2. optionally refine those maps with the vendored CAT12 AMAP C core once the
   compiled extension is available.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

from mrsiprep.io.naming import anat_derivative
from mrsiprep.tissue.ants_atropos import segment_t1_atropos
from mrsiprep.tissue.cat12_amap import AmapParameters, run_cat12_amap
from mrsiprep.tissue.cat12_cleanup import (
    Cat12CleanupParameters,
    clean_final_with_atlas,
    clean_gwc,
    correct_outer_csf_with_atlas,
)
from mrsiprep.tissue.tpm_gmm import (
    CAT12_TPM_LABELS,
    Cat12TemplateAssets,
    DEFAULT_ANTS_THREADS,
    TpmGmmParameters,
    initialize_tissues_with_cat12_tpm,
    initial_labels_from_probabilities as _cat12_initial_labels_from_probabilities,
)
from mrsiprep.utils.images import save_nifti
from mrsiprep.utils.misc import normalize_session, normalize_subject


TISSUE_LABELS = ("GM", "WM", "CSF")


@dataclass(frozen=True)
class TissueSegmentationResult:
    """Native-space tissue probability maps."""

    probabilities: dict[str, Path]
    backend: str
    refined_with_amap: bool = False


@dataclass(frozen=True)
class Cat12PreAmapParameters:
    """Controls for CAT's MATLAB-side preparation before ``cat_amap``."""

    enabled: bool = True
    support_threshold: float = 0.05
    csf_floor: bool = True
    csf_floor_value: float = 0.33
    csf_floor_scale: float = 0.8
    add_csf_noise: bool = True
    random_seed: int = 0
    auto_mrf: bool = True
    use_class_probabilities: bool = False
    seed_min_probability: float = 0.0
    csf_seed_min_probability: float = 0.0


@dataclass(frozen=True)
class Cat12PreAmapResult:
    normalized: np.ndarray
    prepared: np.ndarray
    labels: np.ndarray
    brain_mask: np.ndarray
    brain_mask0: np.ndarray
    csf_floor: np.ndarray
    mrf_weight: float
    parameters: AmapParameters


@dataclass
class CAT12LikeTissueSegmenter:
    """Standalone callable segmenter for raw T1w images."""

    output_dir: Path
    overwrite: bool = False
    verbose: bool = False
    initializer: str = "cat12-tpm-gmm"
    refine_with_amap: bool = False
    amap_parameters: AmapParameters = field(default_factory=AmapParameters)
    pre_amap_parameters: Cat12PreAmapParameters = field(default_factory=Cat12PreAmapParameters)
    cleanup_parameters: Cat12CleanupParameters = field(default_factory=Cat12CleanupParameters)
    tpm_gmm_parameters: TpmGmmParameters = field(default_factory=TpmGmmParameters)
    cat12_assets: Cat12TemplateAssets = field(default_factory=Cat12TemplateAssets)

    def segment_t1(
        self,
        t1_path: str | Path,
        subject: str = "standalone",
        session: str | None = None,
        brain_mask: str | Path | None = None,
    ) -> TissueSegmentationResult:
        """Segment a raw T1w image into native-space GM/WM/CSF probability maps."""

        t1_path = Path(t1_path)
        brain_mask_path = Path(brain_mask) if brain_mask is not None else None
        cfg = SimpleNamespace(
            derivative_dir=Path(self.output_dir),
            overwrite=self.overwrite,
            verbose=self.verbose,
        )
        subject = normalize_subject(subject)
        session = normalize_session(session)
        if self.initializer in {"cat12-tpm-gmm", "tpm-gmm"}:
            initial_result = initialize_tissues_with_cat12_tpm(
                t1_path=t1_path,
                output_dir=Path(self.output_dir),
                subject=subject,
                session=session,
                brain_mask=brain_mask_path,
                assets=self.cat12_assets,
                parameters=self.tpm_gmm_parameters,
                overwrite=self.overwrite,
                verbose=self.verbose,
            )
            initial = initial_result.probabilities
            initial_classes = initial_result.all_probabilities
            initial_atlas = initial_result.atlas
            initial_backend = "cat12-tpm-gmm"
        elif self.initializer == "atropos":
            initial = segment_t1_atropos(cfg, subject, session, t1_path, brain_mask_path)
            initial_classes = None
            initial_atlas = None
            initial_backend = "ants-atropos"
        else:
            raise ValueError(f"Unsupported tissue initializer: {self.initializer}")
        if not self.refine_with_amap:
            return TissueSegmentationResult(initial, backend=initial_backend, refined_with_amap=False)

        refined = refine_initial_maps_with_amap(
            t1_path=t1_path,
            initial_probabilities=initial,
            class_probabilities=initial_classes,
            output_dir=Path(self.output_dir),
            subject=subject,
            session=session,
            parameters=self.amap_parameters,
            pre_amap_parameters=self.pre_amap_parameters,
            cleanup_parameters=self.cleanup_parameters,
            atlas_path=initial_atlas,
            output_desc="cat12TpmAmap" if initial_backend == "cat12-tpm-gmm" else "cat12Amap",
            write_intermediates=True,
        )
        return TissueSegmentationResult(refined, backend=f"{initial_backend}+cat12-amap", refined_with_amap=True)


def refine_initial_maps_with_amap(
    t1_path: str | Path,
    initial_probabilities: dict[str, Path],
    output_dir: str | Path,
    subject: str,
    session: str | None,
    class_probabilities: dict[str, Path] | None = None,
    parameters: AmapParameters | None = None,
    pre_amap_parameters: Cat12PreAmapParameters | None = None,
    cleanup_parameters: Cat12CleanupParameters | None = None,
    atlas_path: str | Path | None = None,
    normalization: str = "cat-las",
    output_desc: str | None = None,
    write_intermediates: bool = False,
) -> dict[str, Path]:
    """Refine preliminary maps with CAT12 AMAP and write GM/WM/CSF maps."""

    t1_img = nib.load(str(t1_path))
    t1 = np.nan_to_num(t1_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    probs = _load_probability_stack(initial_probabilities)
    class_probs = _load_probability_stack(class_probabilities, labels=CAT12_TPM_LABELS) if class_probabilities is not None else None
    src = _normalize_t1_for_amap(t1, probs, mode=normalization)
    voxel_size = tuple(float(v) for v in t1_img.header.get_zooms()[:3])
    pre_params = pre_amap_parameters or Cat12PreAmapParameters()
    pre_amap = prepare_cat12_amap_input(
        normalized=src,
        probabilities=probs,
        voxel_size=voxel_size,
        amap_parameters=parameters,
        pre_amap_parameters=pre_params,
        class_probabilities=class_probs,
    )
    amap = run_cat12_amap(pre_amap.prepared, pre_amap.labels, voxel_size, pre_amap.parameters)
    cleanup = cleanup_parameters or Cat12CleanupParameters()
    gwc_probabilities = clean_gwc(amap.probabilities, voxel_size=voxel_size, parameters=cleanup) if cleanup.enabled else amap.probabilities
    atlas_native_path = Path(atlas_path) if atlas_path is not None else None
    atlas_cleanup_applied = False
    outer_csf_applied = False
    if cleanup.enabled and atlas_native_path is not None and atlas_native_path.exists():
        atlas = np.nan_to_num(nib.load(str(atlas_native_path)).get_fdata(dtype=np.float32).squeeze(), copy=False)
        if cleanup.atlas_cleanup:
            final_probabilities = clean_final_with_atlas(
                gwc_probabilities,
                normalized_t1=src * 3.0,
                atlas_labels=atlas,
                voxel_size=voxel_size,
                parameters=cleanup,
            )
            atlas_cleanup_applied = True
            outer_csf_applied = cleanup.outer_csf_correction
        elif cleanup.outer_csf_correction:
            final_probabilities = correct_outer_csf_with_atlas(
                gwc_probabilities,
                normalized_t1=src * 3.0,
                atlas_labels=atlas,
                voxel_size=voxel_size,
                parameters=cleanup,
            )
            outer_csf_applied = True
        else:
            final_probabilities = gwc_probabilities
    else:
        final_probabilities = gwc_probabilities

    if write_intermediates:
        desc_prefix = output_desc or "cat12Amap"
        input_path = anat_derivative(
            Path(output_dir),
            subject,
            session,
            space="T1w",
            desc=f"{desc_prefix}Input{normalization.replace('-', '')}",
            suffix_override="T1w",
        )
        prepared_path = anat_derivative(
            Path(output_dir),
            subject,
            session,
            space="T1w",
            desc=f"{desc_prefix}Ymib",
            suffix_override="T1w",
        )
        seed_path = anat_derivative(
            Path(output_dir),
            subject,
            session,
            space="T1w",
            desc=f"{desc_prefix}Seed",
            suffix_override="dseg",
        )
        save_nifti(src.astype(np.float32), t1_img, input_path, dtype=np.float32)
        save_nifti(pre_amap.prepared.astype(np.float32), t1_img, prepared_path, dtype=np.float32)
        save_nifti(pre_amap.labels.astype(np.uint8), t1_img, seed_path, dtype=np.uint8)
        save_nifti(
            pre_amap.brain_mask.astype(np.uint8),
            t1_img,
            anat_derivative(Path(output_dir), subject, session, space="T1w", desc=f"{desc_prefix}Yb", suffix_override="mask"),
            dtype=np.uint8,
        )
        save_nifti(
            pre_amap.brain_mask0.astype(np.uint8),
            t1_img,
            anat_derivative(Path(output_dir), subject, session, space="T1w", desc=f"{desc_prefix}Yb0", suffix_override="mask"),
            dtype=np.uint8,
        )
        save_nifti(
            pre_amap.csf_floor.astype(np.float32),
            t1_img,
            anat_derivative(Path(output_dir), subject, session, space="T1w", desc=f"{desc_prefix}CSFFloor", suffix_override="T1w"),
            dtype=np.float32,
        )
        for idx, label_name in enumerate(TISSUE_LABELS):
            raw_out = anat_derivative(
                Path(output_dir),
                subject,
                session,
                space="T1w",
                label=label_name,
                desc=f"{desc_prefix}Raw",
                suffix_override="probseg",
            )
            save_nifti(amap.probabilities[..., idx].astype(np.float32), t1_img, raw_out, dtype=np.float32)
            gwc_out = anat_derivative(
                Path(output_dir),
                subject,
                session,
                space="T1w",
                label=label_name,
                desc=f"{desc_prefix}Gwc",
                suffix_override="probseg",
            )
            save_nifti(gwc_probabilities[..., idx].astype(np.float32), t1_img, gwc_out, dtype=np.float32)
        if outer_csf_applied:
            outer_csf_mask = (
                (final_probabilities[..., 2] > gwc_probabilities[..., 2] + 1e-6)
                & (final_probabilities[..., 0] + final_probabilities[..., 1] < gwc_probabilities[..., 0] + gwc_probabilities[..., 1])
            )
            outer_csf_path = anat_derivative(
                Path(output_dir),
                subject,
                session,
                space="T1w",
                desc=f"{desc_prefix}OuterCSF",
                suffix_override="mask",
            )
            save_nifti(outer_csf_mask.astype(np.uint8), t1_img, outer_csf_path, dtype=np.uint8)
        cleanup_path = anat_derivative(Path(output_dir), subject, session, desc=f"{desc_prefix}Cleanup", suffix_override="json")
        cleanup_path.parent.mkdir(parents=True, exist_ok=True)
        cleanup_path.write_text(
            json.dumps(
                {
                    "enabled": cleanup.enabled,
                    "cleanup_strength": cleanup.cleanup_strength,
                    "cleanup_level": cleanup.level(voxel_size),
                    "extra_cleanup": cleanup.extra_cleanup,
                    "atlas_cleanup_enabled": cleanup.atlas_cleanup,
                    "atlas_cleanup_applied": atlas_cleanup_applied,
                    "atlas_native": str(atlas_native_path) if atlas_native_path is not None else None,
                    "outer_csf_applied": outer_csf_applied,
                    "final_cleanup_strength": cleanup.final_cleanup_strength(voxel_size),
                    "final_cleanup_distance": cleanup.final_cleanup_distance(),
                    "outer_csf_correction": cleanup.outer_csf_correction,
                    "outer_csf_distance_mm": cleanup.outer_csf_distance_mm,
                    "outer_csf_intensity": cleanup.outer_csf_intensity,
                    "source": "cat_main_clean_gwc",
                    "atlas_source": "cat_main_cleanup" if atlas_cleanup_applied else None,
                    "pre_amap": {
                        "enabled": pre_params.enabled,
                        "support_threshold": pre_params.support_threshold,
                        "csf_floor": pre_params.csf_floor,
                        "csf_floor_value": pre_params.csf_floor_value,
                        "csf_floor_scale": pre_params.csf_floor_scale,
                        "add_csf_noise": pre_params.add_csf_noise,
                        "auto_mrf": pre_params.auto_mrf,
                        "use_class_probabilities": pre_params.use_class_probabilities,
                        "seed_min_probability": pre_params.seed_min_probability,
                        "csf_seed_min_probability": pre_params.csf_seed_min_probability,
                        "mrf_weight": pre_amap.mrf_weight,
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    outputs = {}
    for idx, label_name in enumerate(TISSUE_LABELS):
        out = anat_derivative(
            Path(output_dir),
            subject,
            session,
            space="T1w",
            label=label_name,
            desc=output_desc,
            suffix_override="probseg",
        )
        outputs[label_name] = save_nifti(final_probabilities[..., idx].astype(np.float32), t1_img, out, dtype=np.float32)
    return outputs


def run_amap_from_cat_reference(
    t1_path: str | Path,
    p1_path: str | Path,
    p2_path: str | Path,
    p3_path: str | Path,
    output_dir: str | Path,
    subject: str,
    session: str | None,
    p0_path: str | Path | None = None,
    parameters: AmapParameters | None = None,
    pre_amap_parameters: Cat12PreAmapParameters | None = None,
    cleanup_parameters: Cat12CleanupParameters | None = None,
    normalization: str = "cat-las",
) -> tuple[dict[str, Path], dict, Path]:
    """Run AMAP using CAT p maps as seed/reference and write validation metrics."""

    t1_path = Path(t1_path)
    output_dir = Path(output_dir)
    t1_img = nib.load(str(t1_path))
    t1 = np.nan_to_num(t1_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    references = {
        "GM": Path(p1_path),
        "WM": Path(p2_path),
        "CSF": Path(p3_path),
    }
    probs = _load_probability_stack(references)
    src = _normalize_t1_for_amap(t1, probs, mode=normalization)
    seed_labels = None
    if p0_path is not None:
        seed_labels = _initial_labels_from_p0(Path(p0_path), t1_img.shape[:3])
    voxel_size = tuple(float(v) for v in t1_img.header.get_zooms()[:3])
    pre_params = pre_amap_parameters or Cat12PreAmapParameters()
    pre_amap = prepare_cat12_amap_input(
        normalized=src,
        probabilities=probs,
        voxel_size=voxel_size,
        amap_parameters=parameters,
        pre_amap_parameters=pre_params,
        seed_labels=seed_labels,
    )
    amap = run_cat12_amap(pre_amap.prepared, pre_amap.labels, voxel_size, pre_amap.parameters)
    cleanup = cleanup_parameters or Cat12CleanupParameters()
    final_probabilities = (
        clean_gwc(amap.probabilities, voxel_size=voxel_size, parameters=cleanup)
        if cleanup.enabled
        else amap.probabilities
    )

    outputs = {}
    input_path = anat_derivative(output_dir, subject, session, space="T1w", desc=f"catAmapInput{normalization.replace('-', '')}", suffix_override="T1w")
    prepared_path = anat_derivative(output_dir, subject, session, space="T1w", desc="catAmapYmib", suffix_override="T1w")
    seed_path = anat_derivative(output_dir, subject, session, space="T1w", desc="catAmapSeed", suffix_override="dseg")
    save_nifti(src.astype(np.float32), t1_img, input_path, dtype=np.float32)
    save_nifti(pre_amap.prepared.astype(np.float32), t1_img, prepared_path, dtype=np.float32)
    save_nifti(pre_amap.labels.astype(np.uint8), t1_img, seed_path, dtype=np.uint8)
    save_nifti(
        pre_amap.brain_mask.astype(np.uint8),
        t1_img,
        anat_derivative(output_dir, subject, session, space="T1w", desc="catAmapYb", suffix_override="mask"),
        dtype=np.uint8,
    )
    save_nifti(
        pre_amap.brain_mask0.astype(np.uint8),
        t1_img,
        anat_derivative(output_dir, subject, session, space="T1w", desc="catAmapYb0", suffix_override="mask"),
        dtype=np.uint8,
    )
    save_nifti(
        pre_amap.csf_floor.astype(np.float32),
        t1_img,
        anat_derivative(output_dir, subject, session, space="T1w", desc="catAmapCSFFloor", suffix_override="T1w"),
        dtype=np.float32,
    )

    metrics = {
        "backend": "cat12-amap",
        "reference": "cat12-pmaps",
        "normalization": normalization,
        "cleanup": {
            "enabled": cleanup.enabled,
            "cleanup_strength": cleanup.cleanup_strength,
            "cleanup_level": cleanup.level(voxel_size),
            "extra_cleanup": cleanup.extra_cleanup,
            "source": "cat_main_clean_gwc",
        },
        "amap_input": str(input_path),
        "amap_prepared": str(prepared_path),
        "seed_labels": str(seed_path),
        "pre_amap": {
            "enabled": pre_params.enabled,
            "support_threshold": pre_params.support_threshold,
            "csf_floor": pre_params.csf_floor,
            "csf_floor_value": pre_params.csf_floor_value,
            "csf_floor_scale": pre_params.csf_floor_scale,
            "add_csf_noise": pre_params.add_csf_noise,
            "auto_mrf": pre_params.auto_mrf,
            "use_class_probabilities": pre_params.use_class_probabilities,
            "seed_min_probability": pre_params.seed_min_probability,
            "csf_seed_min_probability": pre_params.csf_seed_min_probability,
            "mrf_weight": pre_amap.mrf_weight,
        },
        "classes": {},
    }
    for idx, label_name in enumerate(TISSUE_LABELS):
        out = anat_derivative(output_dir, subject, session, space="T1w", label=label_name, desc="catAmapSeeded", suffix_override="probseg")
        raw_out = anat_derivative(output_dir, subject, session, space="T1w", label=label_name, desc="catAmapSeededRaw", suffix_override="probseg")
        save_nifti(amap.probabilities[..., idx].astype(np.float32), t1_img, raw_out, dtype=np.float32)
        outputs[label_name] = save_nifti(final_probabilities[..., idx].astype(np.float32), t1_img, out, dtype=np.float32)
        metrics["classes"][label_name] = _probability_metrics(final_probabilities[..., idx], probs[..., idx], voxel_size)

    metrics["means"] = [float(v) for v in amap.means]
    metrics["stds"] = [float(v) for v in amap.stds]
    metrics_path = anat_derivative(output_dir, subject, session, desc="catAmapValidation", suffix_override="json")
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return outputs, metrics, metrics_path


def prepare_cat12_amap_input(
    normalized: np.ndarray,
    probabilities: np.ndarray,
    voxel_size: tuple[float, float, float],
    amap_parameters: AmapParameters | None = None,
    pre_amap_parameters: Cat12PreAmapParameters | None = None,
    seed_labels: np.ndarray | None = None,
    class_probabilities: np.ndarray | None = None,
) -> Cat12PreAmapResult:
    """Prepare ``Ymib`` and ``Yp0b`` similarly to CAT's ``cat_main_amap``."""

    params = pre_amap_parameters or Cat12PreAmapParameters()
    amap_params = amap_parameters or AmapParameters()
    src = np.nan_to_num(np.asarray(normalized, dtype=np.float64).squeeze(), copy=False)
    probs = _validate_probability_array(probabilities)
    if src.shape != probs.shape[:3]:
        raise ValueError(f"Normalized image shape {src.shape} does not match probability shape {probs.shape[:3]}.")
    class_probs = (
        _validate_class_probability_array(class_probabilities, src.shape)
        if class_probabilities is not None and params.use_class_probabilities
        else None
    )
    labels_seed = _validate_seed_labels(seed_labels, src.shape) if seed_labels is not None else None
    if not params.enabled:
        labels = labels_seed if labels_seed is not None else _initial_labels_for_cat12_amap(probs, params, class_probs)
        prepared = np.clip(src, 0.0, 2.0)
        prepared[labels == 0] = 0.0
        return Cat12PreAmapResult(
            normalized=src.astype(np.float64),
            prepared=np.round(prepared * 1e4) / 1e4,
            labels=labels,
            brain_mask=labels > 0,
            brain_mask0=labels > 0,
            csf_floor=np.zeros(src.shape, dtype=np.float32),
            mrf_weight=float(amap_params.mrf_weight),
            parameters=amap_params,
        )

    if class_probs is not None and labels_seed is None:
        brain_mask0 = _cat12_tissue_support_from_classes(class_probs, params)
    else:
        brain_mask0 = np.sum(probs, axis=-1) > params.support_threshold
    brain_mask = _cat12_brain_mask(brain_mask0)
    labels = labels_seed if labels_seed is not None else _initial_labels_for_cat12_amap(probs, params, class_probs)
    labels = labels.astype(np.uint8, copy=False)
    labels[~brain_mask] = 0

    prepared = np.clip(src, 0.0, 2.0)
    prepared[~brain_mask] = 0.0
    prepared = np.round(prepared * 1e4) / 1e4

    csf_floor = np.zeros(prepared.shape, dtype=np.float32)
    if params.csf_floor:
        csf_floor = _cat12_csf_floor(prepared, labels, brain_mask, voxel_size, params)
        prepared = np.maximum(prepared, csf_floor)
        prepared[~brain_mask] = 0.0
        prepared = np.round(np.clip(prepared, 0.0, 2.0) * 1e4) / 1e4

    mrf_weight = float(amap_params.mrf_weight)
    if params.auto_mrf and mrf_weight == 0.0:
        mrf_weight = _estimate_cat12_auto_mrf(prepared, probs, voxel_size)
    resolved_amap = replace(amap_params, mrf_weight=mrf_weight)
    return Cat12PreAmapResult(
        normalized=src.astype(np.float64),
        prepared=prepared.astype(np.float64),
        labels=labels,
        brain_mask=brain_mask,
        brain_mask0=brain_mask0,
        csf_floor=csf_floor.astype(np.float32),
        mrf_weight=mrf_weight,
        parameters=resolved_amap,
    )


def _load_probability_stack(probabilities: dict[str, Path], labels: tuple[str, ...] = TISSUE_LABELS) -> np.ndarray:
    missing = [label for label in labels if label not in probabilities]
    if missing:
        raise ValueError(f"Missing initial tissue probability maps: {', '.join(missing)}")
    arrays = []
    ref_shape = None
    for label in labels:
        data = np.nan_to_num(nib.load(str(probabilities[label])).get_fdata(dtype=np.float32).squeeze(), copy=False)
        if ref_shape is None:
            ref_shape = data.shape
        elif data.shape != ref_shape:
            raise ValueError("Initial tissue probability maps do not share the same shape.")
        arrays.append(np.clip(data, 0.0, 1.0))
    return np.stack(arrays, axis=-1)


def _validate_probability_array(probabilities: np.ndarray) -> np.ndarray:
    probs = np.asarray(probabilities, dtype=np.float32)
    if probs.ndim != 4 or probs.shape[-1] != 3:
        raise ValueError(f"Expected GM/WM/CSF probability stack, got {probs.shape}")
    return np.clip(np.nan_to_num(probs, copy=False), 0.0, 1.0)


def _validate_class_probability_array(probabilities: np.ndarray | None, expected_shape: tuple[int, int, int]) -> np.ndarray:
    if probabilities is None:
        raise ValueError("class probabilities cannot be None")
    probs = np.asarray(probabilities, dtype=np.float32)
    if probs.ndim != 4 or probs.shape[-1] < len(CAT12_TPM_LABELS):
        raise ValueError(f"Expected six-class CAT12 probability stack, got {probs.shape}")
    if probs.shape[:3] != expected_shape:
        raise ValueError(f"Class probability shape {probs.shape[:3]} does not match expected shape {expected_shape}.")
    return np.clip(np.nan_to_num(probs[..., : len(CAT12_TPM_LABELS)], copy=False), 0.0, 1.0)


def _cat12_reordered_classes(class_probabilities: np.ndarray) -> np.ndarray:
    # CAT uses max(Ycls(:,:,:,[3,1,2,4:Kb2])) before mapping to AMAP labels.
    return np.stack(
        [
            class_probabilities[..., 2],  # CSF
            class_probabilities[..., 0],  # GM
            class_probabilities[..., 1],  # WM
            class_probabilities[..., 3],  # bone -> CSF-like label
            class_probabilities[..., 4],  # soft tissue -> background
            class_probabilities[..., 5],  # background
        ],
        axis=-1,
    )


def _initial_labels_from_cat12_classes(class_probabilities: np.ndarray, params: Cat12PreAmapParameters) -> np.ndarray:
    reordered = _cat12_reordered_classes(class_probabilities)
    max_prob = np.max(reordered, axis=-1)
    max_idx = np.argmax(reordered, axis=-1)
    cat_label_map = np.asarray([1, 2, 3, 1, 0, 0], dtype=np.uint8)
    labels = cat_label_map[max_idx]
    labels[max_prob <= params.support_threshold] = 0
    return labels.astype(np.uint8, copy=False)


def _cat12_tissue_support_from_classes(class_probabilities: np.ndarray, params: Cat12PreAmapParameters) -> np.ndarray:
    return _initial_labels_from_cat12_classes(class_probabilities, params) > 0


def _initial_labels_for_cat12_amap(
    probabilities: np.ndarray,
    params: Cat12PreAmapParameters,
    class_probabilities: np.ndarray | None = None,
) -> np.ndarray:
    if class_probabilities is not None:
        labels = _initial_labels_from_cat12_classes(class_probabilities, params)
    else:
        labels = _initial_labels_from_probabilities(probabilities, min_support=params.support_threshold)
    seed_min = max(0.0, float(params.seed_min_probability))
    csf_seed_min = max(0.0, float(params.csf_seed_min_probability))
    if seed_min == 0.0 and csf_seed_min == 0.0:
        return labels

    probs = _validate_probability_array(probabilities)
    cat_order = np.stack([probs[..., 2], probs[..., 0], probs[..., 1]], axis=-1)
    if seed_min > 0.0:
        labels[np.max(cat_order, axis=-1) < seed_min] = 0
    if csf_seed_min > 0.0:
        labels[(labels == 1) & (probs[..., 2] < csf_seed_min)] = 0
    return labels


def _validate_seed_labels(seed_labels: np.ndarray | None, expected_shape: tuple[int, int, int]) -> np.ndarray:
    if seed_labels is None:
        raise ValueError("seed_labels cannot be None")
    labels = np.asarray(seed_labels, dtype=np.uint8).squeeze()
    if labels.shape != expected_shape:
        raise ValueError(f"Seed label shape {labels.shape} does not match expected shape {expected_shape}.")
    if np.any(labels > 3):
        raise ValueError("Seed labels must use 0=background, 1=CSF, 2=GM, 3=WM.")
    return labels.copy()


def _cat12_brain_mask(support: np.ndarray) -> np.ndarray:
    mask = np.asarray(support, dtype=bool)
    if not np.any(mask):
        return mask
    structure = ndi.generate_binary_structure(3, 2)
    mask = ndi.binary_closing(mask, structure=structure, iterations=1)
    mask = ndi.binary_fill_holes(mask)
    labels, n_labels = ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))
    if n_labels == 0:
        return mask
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return labels == int(np.argmax(counts))


def _cat12_csf_floor(
    prepared: np.ndarray,
    labels: np.ndarray,
    brain_mask: np.ndarray,
    voxel_size: tuple[float, float, float],
    params: Cat12PreAmapParameters,
) -> np.ndarray:
    ycsf = params.csf_floor_value * brain_mask.astype(np.float32)
    ycsf = _smooth_fwhm(ycsf, tuple(0.6 * float(v) for v in voxel_size), voxel_size)
    if params.add_csf_noise:
        rng = np.random.default_rng(params.random_seed)
        noise1, noise2 = _estimate_csf_noise(prepared, labels, voxel_size)
        rand = rng.standard_normal(prepared.shape).astype(np.float32)
        ycsf = ycsf + _smooth_fwhm(rand, (0.5, 0.5, 0.5), voxel_size) * max(0.005, min(0.2, noise1 / 4.0))
        rand = rng.standard_normal(prepared.shape).astype(np.float32)
        ycsf = ycsf + _smooth_fwhm(rand, (1.0, 1.0, 1.0), voxel_size) * max(0.005, min(0.2, noise2))
    support_smooth = _smooth_fwhm((ycsf > 0).astype(np.float32), (2.0, 2.0, 2.0), voxel_size)
    floor = np.maximum(0.0, ycsf * params.csf_floor_scale * support_smooth)
    floor[~brain_mask] = 0.0
    return floor.astype(np.float32)


def _estimate_csf_noise(prepared: np.ndarray, labels: np.ndarray, voxel_size: tuple[float, float, float]) -> tuple[float, float]:
    wm = labels == 3
    if np.count_nonzero(wm) < 100:
        return 0.02, 0.02
    local_mean = _masked_gaussian_mean(prepared.astype(np.float32), wm, (1.0, 1.0, 1.0), voxel_size)
    residual = prepared[wm] - local_mean[wm]
    noise = _robust_std(residual) / max(float(np.mean(voxel_size)), np.finfo(np.float32).eps)
    return float(noise), float(noise)


def _estimate_cat12_auto_mrf(prepared: np.ndarray, probabilities: np.ndarray, voxel_size: tuple[float, float, float]) -> float:
    std_values = []
    for idx in (0, 1):  # CAT estimates from high-confidence GM and WM.
        mask = probabilities[..., idx] > 0.94
        if np.count_nonzero(mask) < 100:
            mask = probabilities[..., idx] > 0.75
        if np.count_nonzero(mask) < 100:
            continue
        local_mean = _masked_gaussian_mean(prepared.astype(np.float32), mask, (2.0, 2.0, 2.0), voxel_size)
        local_var = _masked_gaussian_mean(((prepared - local_mean) ** 2).astype(np.float32), mask, (2.0, 2.0, 2.0), voxel_size)
        vals = np.sqrt(np.maximum(local_var[mask], 0.0))
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size:
            std_values.append(float(np.nanmean(vals)))
    if not std_values:
        return 0.0
    return float(min(0.15, 3.0 * float(np.nanmean(std_values))) * 0.5)


def _masked_gaussian_mean(
    values: np.ndarray,
    mask: np.ndarray,
    fwhm_mm: tuple[float, float, float],
    voxel_size: tuple[float, float, float],
) -> np.ndarray:
    weights = mask.astype(np.float32)
    numerator = _smooth_fwhm(np.nan_to_num(values * weights, copy=False), fwhm_mm, voxel_size)
    denominator = _smooth_fwhm(weights, fwhm_mm, voxel_size)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator, dtype=np.float32), where=denominator > 1e-6)


def _smooth_fwhm(
    values: np.ndarray,
    fwhm_mm: tuple[float, float, float],
    voxel_size: tuple[float, float, float],
) -> np.ndarray:
    fwhm = np.asarray(fwhm_mm, dtype=np.float32)
    voxel = np.asarray(voxel_size, dtype=np.float32)
    sigma = fwhm / np.maximum(voxel, np.finfo(np.float32).eps) / np.sqrt(8.0 * np.log(2.0))
    return ndi.gaussian_filter(values.astype(np.float32, copy=False), sigma=sigma, mode="nearest")


def _robust_std(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=np.float32)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0
    med = float(np.nanmedian(vals))
    mad = float(np.nanmedian(np.abs(vals - med)))
    return 1.4826 * mad


def _initial_labels_from_probabilities(probabilities: np.ndarray, min_support: float = 0.05) -> np.ndarray:
    """Build CAT AMAP labels: 0=background, 1=CSF, 2=GM, 3=WM."""

    return _cat12_initial_labels_from_probabilities(probabilities, min_support=min_support)


def _initial_labels_from_p0(p0_path: Path, expected_shape: tuple[int, int, int]) -> np.ndarray:
    p0 = np.nan_to_num(nib.load(str(p0_path)).get_fdata(dtype=np.float32).squeeze(), copy=False)
    if p0.shape != expected_shape:
        raise ValueError(f"p0 shape {p0.shape} does not match T1 shape {expected_shape}.")
    return np.clip(np.rint(p0), 0, 3).astype(np.uint8)


def _normalize_t1_for_amap(t1: np.ndarray, probabilities: np.ndarray, mode: str = "cat-las") -> np.ndarray:
    """Approximate CAT's intensity-normalized AMAP input.

    ``cat-global`` follows the threshold mapping style in ``cat_main_gintnorm``.
    ``cat-las`` adds a lightweight local approximation of
    ``cat_main_LASsimple``.
    """

    if mode == "anchors":
        return _normalize_t1_with_class_anchors(t1, probabilities)
    peaks = _estimate_tissue_peaks(t1, probabilities)
    if mode == "cat-global":
        return _cat_global_intensity_normalization(t1, peaks)
    if mode == "cat-las":
        return _cat_las_simple_normalization(t1, probabilities, peaks)
    raise ValueError(f"Unsupported AMAP normalization mode: {mode}")


def _normalize_t1_with_class_anchors(t1: np.ndarray, probabilities: np.ndarray) -> np.ndarray:
    """Simple class-median interpolation kept as a baseline."""

    gm = probabilities[..., 0] > 0.5
    wm = probabilities[..., 1] > 0.5
    csf = probabilities[..., 2] > 0.5
    anchors = []
    targets = []
    for mask, target in ((csf, 1.0 / 3.0), (gm, 2.0 / 3.0), (wm, 1.0)):
        values = t1[mask & np.isfinite(t1)]
        if values.size:
            anchors.append(float(np.nanmedian(values)))
            targets.append(target)
    if len(anchors) >= 2:
        order = np.argsort(anchors)
        anchors_arr = np.asarray(anchors, dtype=np.float32)[order]
        targets_arr = np.asarray(targets, dtype=np.float32)[order]
        normalized = np.interp(t1, anchors_arr, targets_arr, left=0.0, right=1.2)
    else:
        finite = t1[np.isfinite(t1)]
        lo, hi = np.percentile(finite, [1, 99]) if finite.size else (0.0, 1.0)
        scale = max(float(hi - lo), np.finfo(np.float32).eps)
        normalized = (t1 - lo) / scale
    return np.clip(np.nan_to_num(normalized, copy=False), 0.0, 2.0).astype(np.float64)


def _estimate_tissue_peaks(t1: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(t1)
    peaks = {}
    for idx, label in enumerate(TISSUE_LABELS):
        prob = probabilities[..., idx]
        mask = finite & (prob > 0.75) & (np.argmax(probabilities, axis=-1) == idx)
        if np.count_nonzero(mask) < 1000:
            mask = finite & (prob > 0.5)
        if np.count_nonzero(mask) < 100:
            mask = finite & (prob > 0.1)
        values = t1[mask]
        if values.size == 0:
            values = t1[finite]
        peaks[label] = float(np.nanmedian(values)) if values.size else 0.0
    return peaks


def _cat_global_intensity_normalization(t1: np.ndarray, peaks: dict[str, float]) -> np.ndarray:
    finite = np.isfinite(t1)
    if not np.any(finite):
        return np.zeros_like(t1, dtype=np.float64)
    bg_min = float(np.nanpercentile(t1[finite], 0.1))
    bg_con = max(bg_min, peaks["CSF"] - float(np.mean(np.diff([peaks["CSF"], peaks["GM"], peaks["WM"]]))))
    bm_th = max(bg_min, min(bg_con, peaks["CSF"] - (peaks["GM"] - peaks["CSF"])))
    bm_csf_th = min(bg_con, 0.5 * (bm_th + peaks["CSF"]))
    wm_plus = peaks["WM"] + 0.5 * (peaks["WM"] - peaks["CSF"])
    hi = max(wm_plus, float(np.nanpercentile(t1[finite], 99.9)))
    thresholds = np.asarray([bg_min, bm_th, bm_csf_th, peaks["CSF"], peaks["GM"], peaks["WM"], wm_plus, hi], dtype=np.float64)
    targets = np.asarray([0.0, 0.02, 0.05, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64) / 3.0
    return _piecewise_map(t1, thresholds, targets)


def _cat_las_simple_normalization(t1: np.ndarray, probabilities: np.ndarray, peaks: dict[str, float], las_strength: float = 0.5) -> np.ndarray:
    from scipy import ndimage as ndi

    finite = np.isfinite(t1)
    if not np.any(finite):
        return np.zeros_like(t1, dtype=np.float64)
    sigma = 4.0 * max(0.0, 1.0 - las_strength)
    local = {}
    for idx, label in enumerate(TISSUE_LABELS):
        prob = np.clip(probabilities[..., idx].astype(np.float32), 0.0, 1.0)
        weights = np.clip((prob - 0.35) / 0.65, 0.0, 1.0) ** 2
        numerator = ndi.gaussian_filter(np.nan_to_num(t1 * weights, copy=False), sigma=sigma, mode="nearest")
        denominator = ndi.gaussian_filter(weights, sigma=sigma, mode="nearest")
        lab = np.divide(numerator, denominator, out=np.full_like(t1, peaks[label], dtype=np.float64), where=denominator > 1e-4)
        confident = probabilities[..., idx] > 0.5
        if np.any(confident):
            med = float(np.nanmedian(lab[confident]))
            if np.isfinite(med) and abs(med) > np.finfo(np.float32).eps:
                lab = lab / med * peaks[label]
        local[label] = lab

    bg = float(np.nanpercentile(t1[finite], 0.1))
    gm = local["GM"]
    wm = np.maximum(local["WM"], gm + np.finfo(np.float32).eps)
    csf = np.minimum(local["CSF"], gm - np.finfo(np.float32).eps)
    out = np.zeros_like(t1, dtype=np.float64)

    high = t1 >= wm
    out[high] = 3.0 + (t1[high] - wm[high]) / np.maximum(np.finfo(np.float32).eps, wm[high] - csf[high])
    mid = (t1 >= gm) & (t1 < wm)
    out[mid] = 2.0 + (t1[mid] - gm[mid]) / np.maximum(np.finfo(np.float32).eps, wm[mid] - gm[mid])
    low = (t1 >= csf) & (t1 < gm)
    out[low] = 1.0 + (t1[low] - csf[low]) / np.maximum(np.finfo(np.float32).eps, gm[low] - csf[low])
    background = t1 < csf
    out[background] = (t1[background] - bg) / np.maximum(np.finfo(np.float32).eps, csf[background] - bg)
    return np.clip(np.nan_to_num(out / 3.0, copy=False), 0.0, 2.0).astype(np.float64)


def _piecewise_map(values: np.ndarray, thresholds: np.ndarray, targets: np.ndarray) -> np.ndarray:
    order = np.argsort(thresholds)
    x = thresholds[order]
    y = targets[order]
    keep = np.r_[True, np.diff(x) > np.finfo(np.float32).eps]
    x = x[keep]
    y = y[keep]
    if x.size < 2:
        return np.zeros_like(values, dtype=np.float64)
    out = np.interp(values, x, y)
    above = values > x[-1]
    if np.any(above):
        slope = (y[-1] - y[-2]) / max(np.finfo(np.float32).eps, x[-1] - x[-2])
        out[above] = y[-1] + (values[above] - x[-1]) * slope
    below = values < x[0]
    if np.any(below):
        slope = (y[1] - y[0]) / max(np.finfo(np.float32).eps, x[1] - x[0])
        out[below] = y[0] + (values[below] - x[0]) * slope
    return np.clip(np.nan_to_num(out, copy=False), 0.0, 2.0).astype(np.float64)


def _probability_metrics(candidate: np.ndarray, reference: np.ndarray, voxel_size: tuple[float, float, float]) -> dict[str, float]:
    cand = np.asarray(candidate, dtype=np.float32)
    ref = np.asarray(reference, dtype=np.float32)
    diff = cand - ref
    valid = np.isfinite(cand) & np.isfinite(ref)
    if not np.any(valid):
        return {"mae": float("nan"), "rmse": float("nan"), "corr": float("nan")}
    cand_v = cand[valid]
    ref_v = ref[valid]
    voxel_ml = float(np.prod(voxel_size) / 1000.0)
    if np.std(cand_v) > 0 and np.std(ref_v) > 0:
        corr = float(np.corrcoef(cand_v, ref_v)[0, 1])
    else:
        corr = float("nan")
    return {
        "mae": float(np.mean(np.abs(diff[valid]))),
        "rmse": float(np.sqrt(np.mean(diff[valid] ** 2))),
        "corr": corr,
        "candidate_volume_ml": float(np.sum(cand_v) * voxel_ml),
        "reference_volume_ml": float(np.sum(ref_v) * voxel_ml),
        "volume_delta_ml": float((np.sum(cand_v) - np.sum(ref_v)) * voxel_ml),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mrsiprep-cat12-tissue", description="Standalone CAT12-like T1 tissue segmentation.")
    parser.add_argument("t1", type=Path, help="Raw T1w NIfTI image.")
    parser.add_argument("output_dir", type=Path, help="Output derivative directory.")
    parser.add_argument("--subject", default="standalone")
    parser.add_argument("--session", default=None)
    parser.add_argument("--brain-mask", type=Path, default=None)
    parser.add_argument("--initializer", choices=["cat12-tpm-gmm", "atropos"], default="cat12-tpm-gmm")
    parser.add_argument("--tpm", type=Path, default=None, help="CAT12/SPM-like TPM NIfTI. Defaults to the vendored CAT12 TPM.")
    parser.add_argument("--template-t1", type=Path, default=None, help="Template T1 used for TPM-to-native registration.")
    parser.add_argument("--template-mask", type=Path, default=None, help="Template mask used for TPM-to-native registration.")
    parser.add_argument("--cat-atlas", type=Path, default=None, help="CAT major-region atlas. Defaults to the vendored CAT cat.nii.gz.")
    parser.add_argument(
        "--template-transform",
        default=TpmGmmParameters().registration_transform,
        help=f"ANTs shorthand passed to antsRegistrationSyN. Default: {TpmGmmParameters().registration_transform}.",
    )
    parser.add_argument(
        "--ants-threads",
        type=int,
        default=DEFAULT_ANTS_THREADS,
        help=f"Threads used by ANTs registration/N4/transform calls. Default: {DEFAULT_ANTS_THREADS}.",
    )
    parser.add_argument("--no-n4", action="store_true", help="Skip N4 bias correction before TPM/GMM fitting.")
    parser.add_argument("--use-amap", action="store_true", help="Refine initial maps with the compiled CAT12 AMAP extension.")
    parser.add_argument(
        "--amap-mrf-weight",
        type=float,
        default=AmapParameters().mrf_weight,
        help="Explicit CAT AMAP MRF weight. Keep 0 to allow auto-MRF unless --no-amap-auto-mrf is set.",
    )
    parser.add_argument("--no-pre-amap", action="store_true", help="Skip CAT-style Ymib/Yp0b preparation before AMAP.")
    parser.add_argument(
        "--pre-amap-support-threshold",
        type=float,
        default=Cat12PreAmapParameters().support_threshold,
        help=f"TPM/posterior support threshold used for pre-AMAP masks. Default: {Cat12PreAmapParameters().support_threshold}.",
    )
    parser.add_argument("--no-csf-floor", action="store_true", help="Skip CAT-style low-intensity CSF floor before AMAP.")
    parser.add_argument(
        "--csf-floor-value",
        type=float,
        default=Cat12PreAmapParameters().csf_floor_value,
        help=f"CAT-normalized CSF floor value before scaling. Default: {Cat12PreAmapParameters().csf_floor_value}.",
    )
    parser.add_argument(
        "--csf-floor-scale",
        type=float,
        default=Cat12PreAmapParameters().csf_floor_scale,
        help=f"Scale applied to the CSF floor before AMAP. Default: {Cat12PreAmapParameters().csf_floor_scale}.",
    )
    parser.add_argument("--no-csf-noise", action="store_true", help="Skip CAT-style smoothed CSF-floor noise.")
    parser.add_argument("--no-amap-auto-mrf", action="store_true", help="Do not estimate CAT-style auto-MRF before AMAP.")
    parser.add_argument(
        "--use-six-class-amap-seed",
        action="store_true",
        help="Generate AMAP seed labels from six CAT TPM/GMM classes when available. Off by default; useful for diagnostics.",
    )
    parser.add_argument(
        "--amap-seed-min-probability",
        type=float,
        default=Cat12PreAmapParameters().seed_min_probability,
        help="Minimum winning GM/WM/CSF probability required to keep an auto-generated AMAP seed label. Default: 0.",
    )
    parser.add_argument(
        "--amap-csf-seed-min-probability",
        type=float,
        default=Cat12PreAmapParameters().csf_seed_min_probability,
        help="Minimum CSF probability required to keep an auto-generated CSF AMAP seed label. Default: 0.",
    )
    parser.add_argument("--no-cleanup", action="store_true", help="Skip CAT12-style post-AMAP clean_gwc cleanup.")
    parser.add_argument(
        "--use-atlas-cleanup",
        action="store_true",
        help="Enable the full atlas-guided CAT12 final cleanup approximation after clean_gwc.",
    )
    parser.add_argument("--no-atlas-cleanup", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--no-outer-csf-correction",
        action="store_true",
        help="Skip the outer-boundary WM-to-CSF correction inside atlas cleanup.",
    )
    parser.add_argument(
        "--cleanup-strength",
        type=float,
        default=Cat12CleanupParameters().cleanup_strength,
        help=f"CAT12 cleanupstr value used after AMAP. Default: {Cat12CleanupParameters().cleanup_strength}.",
    )
    parser.add_argument(
        "--cleanup-extra",
        type=int,
        choices=[0, 1, 2, 3],
        default=Cat12CleanupParameters().extra_cleanup,
        help="Enable CAT12 clean_gwc extra cleanup branch 0/1/2/3. Default: 0.",
    )
    parser.add_argument(
        "--outer-csf-distance-mm",
        type=float,
        default=Cat12CleanupParameters().outer_csf_distance_mm,
        help=f"Distance from current tissue support used for outer CSF correction. Default: {Cat12CleanupParameters().outer_csf_distance_mm}.",
    )
    parser.add_argument(
        "--outer-csf-intensity",
        type=float,
        default=Cat12CleanupParameters().outer_csf_intensity,
        help=f"Maximum CAT-normalized intensity eligible for outer CSF correction. Default: {Cat12CleanupParameters().outer_csf_intensity}.",
    )
    parser.add_argument("--seed-p0", type=Path, default=None, help="CAT p0 label map used as AMAP seed labels.")
    parser.add_argument("--seed-p1", type=Path, default=None, help="CAT p1/GM probability map used as seed/reference.")
    parser.add_argument("--seed-p2", type=Path, default=None, help="CAT p2/WM probability map used as seed/reference.")
    parser.add_argument("--seed-p3", type=Path, default=None, help="CAT p3/CSF probability map used as seed/reference.")
    parser.add_argument("--amap-normalization", choices=["anchors", "cat-global", "cat-las"], default="cat-las")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    amap_parameters = AmapParameters(mrf_weight=args.amap_mrf_weight)
    pre_amap_parameters = Cat12PreAmapParameters(
        enabled=not args.no_pre_amap,
        support_threshold=args.pre_amap_support_threshold,
        csf_floor=not args.no_csf_floor,
        csf_floor_value=args.csf_floor_value,
        csf_floor_scale=args.csf_floor_scale,
        add_csf_noise=not args.no_csf_noise,
        auto_mrf=not args.no_amap_auto_mrf,
        use_class_probabilities=args.use_six_class_amap_seed,
        seed_min_probability=args.amap_seed_min_probability,
        csf_seed_min_probability=args.amap_csf_seed_min_probability,
    )
    cat_seed_paths = [args.seed_p1, args.seed_p2, args.seed_p3]
    if any(cat_seed_paths):
        if not all(cat_seed_paths):
            raise SystemExit("--seed-p1, --seed-p2, and --seed-p3 must be provided together.")
        if not args.use_amap:
            raise SystemExit("CAT seed maps are only used with --use-amap.")
        outputs, _, metrics_path = run_amap_from_cat_reference(
            t1_path=args.t1,
            p1_path=args.seed_p1,
            p2_path=args.seed_p2,
            p3_path=args.seed_p3,
            p0_path=args.seed_p0,
            output_dir=args.output_dir,
            subject=args.subject,
            session=args.session,
            parameters=amap_parameters,
            pre_amap_parameters=pre_amap_parameters,
            cleanup_parameters=Cat12CleanupParameters(
                enabled=not args.no_cleanup,
                cleanup_strength=args.cleanup_strength,
                extra_cleanup=args.cleanup_extra,
                atlas_cleanup=args.use_atlas_cleanup and not args.no_atlas_cleanup,
                outer_csf_correction=not args.no_outer_csf_correction,
                outer_csf_distance_mm=args.outer_csf_distance_mm,
                outer_csf_intensity=args.outer_csf_intensity,
            ),
            normalization=args.amap_normalization,
        )
        for label in TISSUE_LABELS:
            print(f"{label}: {outputs[label]}")
        print(f"metrics: {metrics_path}")
        return 0

    segmenter = CAT12LikeTissueSegmenter(
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        verbose=args.verbose,
        initializer=args.initializer,
        refine_with_amap=args.use_amap,
        amap_parameters=amap_parameters,
        pre_amap_parameters=pre_amap_parameters,
        cleanup_parameters=Cat12CleanupParameters(
            enabled=not args.no_cleanup,
            cleanup_strength=args.cleanup_strength,
            extra_cleanup=args.cleanup_extra,
            atlas_cleanup=args.use_atlas_cleanup and not args.no_atlas_cleanup,
            outer_csf_correction=not args.no_outer_csf_correction,
            outer_csf_distance_mm=args.outer_csf_distance_mm,
            outer_csf_intensity=args.outer_csf_intensity,
        ),
        cat12_assets=Cat12TemplateAssets(
            tpm=args.tpm or Cat12TemplateAssets().tpm,
            template=args.template_t1 or Cat12TemplateAssets().template,
            template_mask=args.template_mask if args.template_mask is not None else Cat12TemplateAssets().template_mask,
            atlas=args.cat_atlas if args.cat_atlas is not None else Cat12TemplateAssets().atlas,
        ),
        tpm_gmm_parameters=TpmGmmParameters(
            registration_transform=args.template_transform,
            ants_threads=args.ants_threads,
            n4_bias_correct=not args.no_n4,
        ),
    )
    result = segmenter.segment_t1(args.t1, subject=args.subject, session=args.session, brain_mask=args.brain_mask)
    for label in TISSUE_LABELS:
        print(f"{label}: {result.probabilities[label]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
