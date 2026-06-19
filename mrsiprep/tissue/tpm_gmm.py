"""CAT12 TPM-seeded Gaussian tissue initialization."""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.interfaces.ants import ANTsError, Registration, require_cli, run_cli
from mrsiprep.io.naming import anat_derivative
from mrsiprep.utils.images import save_nifti
from mrsiprep.utils.misc import normalize_session, normalize_subject


TISSUE_LABELS = ("GM", "WM", "CSF")
CAT12_TPM_LABELS = ("GM", "WM", "CSF", "BONE", "SOFT", "BG")
CAT12_GAUSSIANS_PER_CLASS = (1, 1, 2, 3, 4, 2)
CAT12_TEMPLATE_DIR = Path(__file__).resolve().parents[1] / "data" / "templates" / "cat12"
DEFAULT_ANTS_THREADS = max(16, os.cpu_count() or 16)


@dataclass(frozen=True)
class Cat12TemplateAssets:
    """Local CAT12 template/TPM assets used by the Python initializer."""

    tpm: Path = CAT12_TEMPLATE_DIR / "TPM_Age11.5.nii.gz"
    template: Path = CAT12_TEMPLATE_DIR / "T1.nii.gz"
    template_mask: Path | None = CAT12_TEMPLATE_DIR / "brainmask.nii.gz"
    atlas: Path | None = CAT12_TEMPLATE_DIR / "cat.nii.gz"

    def validated(self) -> "Cat12TemplateAssets":
        missing = [
            path
            for path in (self.tpm, self.template, self.template_mask, self.atlas)
            if path is not None and not Path(path).exists()
        ]
        if missing:
            raise FileNotFoundError("Missing CAT12 template assets: " + ", ".join(str(path) for path in missing))
        return Cat12TemplateAssets(
            tpm=Path(self.tpm),
            template=Path(self.template),
            template_mask=Path(self.template_mask) if self.template_mask is not None else None,
            atlas=Path(self.atlas) if self.atlas is not None else None,
        )


@dataclass(frozen=True)
class TpmGmmParameters:
    """Controls for CAT12 TPM-seeded intensity posterior estimation."""

    registration_transform: str = "a"
    ants_threads: int = DEFAULT_ANTS_THREADS
    n4_bias_correct: bool = True
    em_iterations: int = 8
    gaussians_per_class: tuple[int, ...] = CAT12_GAUSSIANS_PER_CLASS
    prior_weight: float = 1.0
    prior_floor: float = 1e-5
    support_threshold: float = 0.03
    robust_percentiles: tuple[float, float] = (0.5, 99.5)
    variance_floor_fraction: float = 1e-4


@dataclass(frozen=True)
class TpmGmmResult:
    """Native-space TPM priors and GMM posteriors."""

    probabilities: dict[str, Path]
    all_probabilities: dict[str, Path]
    priors: dict[str, Path]
    support_mask: Path
    corrected_t1: Path
    seed_labels: Path
    metrics: Path
    transforms: dict[str, list[str]]
    atlas: Path | None = None


def initialize_tissues_with_cat12_tpm(
    t1_path: str | Path,
    output_dir: str | Path,
    subject: str,
    session: str | None,
    brain_mask: str | Path | None = None,
    assets: Cat12TemplateAssets | None = None,
    parameters: TpmGmmParameters | None = None,
    overwrite: bool = False,
    verbose: bool = False,
) -> TpmGmmResult:
    """Create GM/WM/CSF probability maps from warped CAT12 TPM priors.

    ``brain_mask`` is used only as an optional fixed-image registration mask.
    The posterior support remains driven by CAT12 intracranial TPM priors so
    CSF is not removed by a skull-stripping mask.
    """

    t1_path = Path(t1_path)
    output_dir = Path(output_dir)
    subject = normalize_subject(subject)
    session = normalize_session(session)
    params = parameters or TpmGmmParameters()
    asset_paths = (assets or Cat12TemplateAssets()).validated()

    priors, transforms = warp_cat12_tpm_priors_to_t1(
        t1_path=t1_path,
        output_dir=output_dir,
        subject=subject,
        session=session,
        brain_mask=brain_mask,
        assets=asset_paths,
        transform=params.registration_transform,
        ants_threads=params.ants_threads,
        overwrite=overwrite,
        verbose=verbose,
    )
    atlas_path = warp_cat12_atlas_to_t1(
        t1_path=t1_path,
        output_dir=output_dir,
        subject=subject,
        session=session,
        assets=asset_paths,
        transforms=transforms,
        ants_threads=params.ants_threads,
        overwrite=overwrite,
        verbose=verbose,
    )

    t1_img = nib.load(str(t1_path))
    t1 = np.nan_to_num(t1_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    prior_stack = _load_probability_stack(priors, labels=CAT12_TPM_LABELS)
    brain_support = np.sum(prior_stack[..., :3], axis=-1) > params.support_threshold

    support_path = anat_derivative(
        output_dir,
        subject,
        session,
        space="T1w",
        desc="cat12TpmSupport",
        suffix_override="mask",
    )
    save_nifti(brain_support.astype(np.uint8), t1_img, support_path, dtype=np.uint8)

    corrected_path = anat_derivative(output_dir, subject, session, space="T1w", desc="cat12N4", suffix_override="T1w")
    if params.n4_bias_correct:
        _n4_bias_correct(t1_path, support_path, corrected_path, overwrite=overwrite, verbose=verbose, threads=params.ants_threads)
        fit_img = nib.load(str(corrected_path))
        fit_data = np.nan_to_num(fit_img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    else:
        save_nifti(t1.astype(np.float32), t1_img, corrected_path, dtype=np.float32)
        fit_data = t1

    posterior, stats = fit_tpm_weighted_gmm(fit_data, prior_stack, params)
    voxel_ml = float(np.prod(t1_img.header.get_zooms()[:3]) / 1000.0)
    probabilities = {}
    all_probabilities = {}
    for idx, label in enumerate(CAT12_TPM_LABELS):
        out = anat_derivative(
            output_dir,
            subject,
            session,
            space="T1w",
            label=label,
            desc="cat12TpmGmm",
            suffix_override="probseg",
        )
        all_probabilities[label] = save_nifti(posterior[..., idx].astype(np.float32), t1_img, out, dtype=np.float32)
        stats[label]["volume_ml"] = float(np.sum(posterior[..., idx]) * voxel_ml)
    probabilities = {label: all_probabilities[label] for label in TISSUE_LABELS}

    seed = initial_labels_from_probabilities(posterior[..., :3])
    seed_path = anat_derivative(output_dir, subject, session, space="T1w", desc="cat12TpmGmmSeed", suffix_override="dseg")
    save_nifti(seed.astype(np.uint8), t1_img, seed_path, dtype=np.uint8)

    metrics = {
        "backend": "cat12-tpm-gmm",
        "tpm": str(asset_paths.tpm),
        "template": str(asset_paths.template),
        "template_mask": str(asset_paths.template_mask) if asset_paths.template_mask is not None else None,
        "atlas_template": str(asset_paths.atlas) if asset_paths.atlas is not None else None,
        "atlas_native": str(atlas_path) if atlas_path is not None else None,
        "registration_transform": params.registration_transform,
        "registration_direction": "T1w_to_cat12Template",
        "ants_threads": params.ants_threads,
        "ants_reproducible": True,
        "n4_bias_correct": params.n4_bias_correct,
        "gaussians_per_class": list(params.gaussians_per_class),
        "tpm_classes": list(CAT12_TPM_LABELS),
        "transforms": transforms,
        "priors": {label: str(path) for label, path in priors.items()},
        "probabilities": {label: str(path) for label, path in all_probabilities.items()},
        "tissue_probabilities": {label: str(path) for label, path in probabilities.items()},
        "support_mask": str(support_path),
        "corrected_t1": str(corrected_path),
        "seed_labels": str(seed_path),
        "classes": stats,
    }
    metrics_path = anat_derivative(output_dir, subject, session, desc="cat12TpmGmm", suffix_override="json")
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return TpmGmmResult(
        probabilities=probabilities,
        all_probabilities=all_probabilities,
        priors=priors,
        support_mask=support_path,
        corrected_t1=corrected_path,
        seed_labels=seed_path,
        metrics=metrics_path,
        transforms=transforms,
        atlas=atlas_path,
    )


def warp_cat12_tpm_priors_to_t1(
    t1_path: str | Path,
    output_dir: str | Path,
    subject: str,
    session: str | None,
    brain_mask: str | Path | None = None,
    assets: Cat12TemplateAssets | None = None,
    transform: str = "a",
    ants_threads: int = DEFAULT_ANTS_THREADS,
    overwrite: bool = False,
    verbose: bool = False,
) -> tuple[dict[str, Path], dict[str, list[str]]]:
    """Warp CAT12 GM/WM/CSF TPM volumes into native T1 space.

    The registration is estimated as native T1w -> CAT12 template. Template
    priors are then pulled back to native space with the inverse transform set.
    """

    t1_path = Path(t1_path)
    output_dir = Path(output_dir)
    subject = normalize_subject(subject)
    session = normalize_session(session)
    asset_paths = (assets or Cat12TemplateAssets()).validated()
    prefix = _subject_to_template_transform_prefix(output_dir, subject, session)

    existing_transforms = _subject_to_template_transform_set(prefix)
    if existing_transforms["forward"] and not overwrite:
        transforms = existing_transforms
    else:
        transforms = _register_subject_to_template(
            fixed=asset_paths.template,
            moving=t1_path,
            out_prefix=prefix,
            transform=transform,
            threads=ants_threads,
            fixed_mask=asset_paths.template_mask,
            moving_mask=brain_mask,
            verbose=verbose,
        )

    outputs = {
        label: anat_derivative(
            output_dir,
            subject,
            session,
            space="T1w",
            label=label,
            desc="cat12TpmPrior",
            suffix_override="probseg",
        )
        for label in CAT12_TPM_LABELS
    }
    if all(path.exists() for path in outputs.values()) and not overwrite:
        return outputs, transforms

    tpm_img = nib.load(str(asset_paths.tpm))
    tpm = tpm_img.get_fdata(dtype=np.float32)
    if tpm.ndim != 4 or tpm.shape[-1] < len(CAT12_TPM_LABELS):
        raise ValueError(
            f"CAT12 TPM must contain six GM/WM/CSF/bone/soft/background volumes, got shape {tpm.shape}: {asset_paths.tpm}"
        )
    scale = 255.0 if np.nanmax(tpm[..., : len(CAT12_TPM_LABELS)]) > 2.0 else 1.0

    for idx, label in enumerate(CAT12_TPM_LABELS):
        data = np.clip(np.asarray(tpm[..., idx], dtype=np.float32) / scale, 0.0, 1.0)
        moving = _template_prior_path(output_dir, subject, session, label)
        save_nifti(data, _nifti_like(tpm_img, data), moving, dtype=np.float32)
        _warp_template_image_to_t1(
            t1_path,
            moving,
            transforms["inverse"],
            outputs[label],
            interpolation="linear",
            verbose=verbose,
            threads=ants_threads,
        )
        _clip_probability_image(outputs[label])
    return outputs, transforms


def warp_cat12_atlas_to_t1(
    t1_path: str | Path,
    output_dir: str | Path,
    subject: str,
    session: str | None,
    assets: Cat12TemplateAssets | None,
    transforms: dict[str, list[str]],
    ants_threads: int = DEFAULT_ANTS_THREADS,
    overwrite: bool = False,
    verbose: bool = False,
) -> Path | None:
    """Warp CAT's major-region atlas into native T1w space."""

    asset_paths = (assets or Cat12TemplateAssets()).validated()
    if asset_paths.atlas is None:
        return None
    output = anat_derivative(
        Path(output_dir),
        subject,
        session,
        space="T1w",
        desc="cat12Atlas",
        suffix_override="dseg",
    )
    if output.exists() and not overwrite:
        return output
    _warp_template_image_to_t1(
        Path(t1_path),
        asset_paths.atlas,
        transforms["inverse"],
        output,
        interpolation="nearest",
        verbose=verbose,
        threads=ants_threads,
    )
    _round_label_image(output)
    return output


def fit_tpm_weighted_gmm(
    t1: np.ndarray,
    priors: np.ndarray,
    parameters: TpmGmmParameters | None = None,
) -> tuple[np.ndarray, dict[str, dict[str, float]]]:
    """Fit class-conditional Gaussian mixtures with voxelwise TPM priors."""

    params = parameters or TpmGmmParameters()
    image = np.asarray(t1, dtype=np.float32)
    prior_stack = _normalize_prior_stack(priors, params)
    n_classes = prior_stack.shape[-1]
    class_labels = _class_labels(n_classes)
    gaussian_counts = _gaussian_counts(params, n_classes)
    component_classes = np.repeat(np.arange(n_classes, dtype=np.int16), gaussian_counts)
    support = np.sum(np.asarray(priors, dtype=np.float32), axis=-1) > params.support_threshold
    finite = np.isfinite(image) & support
    if np.count_nonzero(finite) < 10:
        raise ValueError("Insufficient finite voxels inside CAT12 TPM support for GMM initialization.")

    lo, hi = np.percentile(image[finite], params.robust_percentiles)
    clipped = np.clip(image, float(lo), float(hi)).astype(np.float32)
    intensity_range = max(float(hi - lo), np.finfo(np.float32).eps)
    var_floor = (params.variance_floor_fraction * intensity_range) ** 2

    values = clipped[finite].astype(np.float32)
    prior_values = prior_stack[finite, :].astype(np.float32)
    means, variances, mixture_weights = _initialize_gmm_components(
        values,
        prior_values,
        gaussian_counts,
        component_classes,
        var_floor,
    )
    for _ in range(max(1, params.em_iterations)):
        means, variances, mixture_weights = _gmm_maximization_step(
            values,
            prior_values,
            means,
            variances,
            mixture_weights,
            component_classes,
            gaussian_counts,
            params,
            var_floor,
        )

    posterior = np.zeros((*image.shape, n_classes), dtype=np.float32)
    posterior[finite, :] = _gmm_class_posteriors(
        values,
        prior_values,
        means,
        variances,
        mixture_weights,
        component_classes,
        n_classes,
        params,
    )
    stats = {}
    for idx, label in enumerate(class_labels):
        class_component_indices = np.flatnonzero(component_classes == idx)
        weights = mixture_weights[class_component_indices].astype(np.float64)
        if weights.size and np.sum(weights) > 0:
            weights = weights / np.sum(weights)
            class_means = means[class_component_indices]
            class_vars = variances[class_component_indices]
            class_mean = float(np.sum(weights * class_means))
            second_moment = float(np.sum(weights * (class_vars + class_means**2)))
            class_std = float(np.sqrt(max(second_moment - class_mean**2, var_floor)))
        else:
            class_mean = float("nan")
            class_std = float("nan")
        stats[label] = {
            "mean": class_mean,
            "std": class_std,
            "gaussians": int(gaussian_counts[idx]),
            "volume_ml": float("nan"),
        }
    return posterior, stats


def _class_labels(n_classes: int) -> tuple[str, ...]:
    if n_classes <= len(CAT12_TPM_LABELS):
        return CAT12_TPM_LABELS[:n_classes]
    return tuple(f"CLASS{idx + 1}" for idx in range(n_classes))


def _gaussian_counts(parameters: TpmGmmParameters, n_classes: int) -> np.ndarray:
    counts = np.asarray(parameters.gaussians_per_class, dtype=np.int16)
    if counts.size < n_classes:
        counts = np.pad(counts, (0, n_classes - counts.size), constant_values=1)
    counts = counts[:n_classes]
    if np.any(counts < 1):
        raise ValueError(f"Gaussian counts must be positive, got {tuple(int(v) for v in counts)}")
    return counts


def _initialize_gmm_components(
    values: np.ndarray,
    priors: np.ndarray,
    gaussian_counts: np.ndarray,
    component_classes: np.ndarray,
    variance_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_components = int(np.sum(gaussian_counts))
    means = np.zeros(n_components, dtype=np.float64)
    variances = np.full(n_components, variance_floor, dtype=np.float64)
    mixture_weights = np.zeros(n_components, dtype=np.float64)
    global_var = max(float(np.var(values.astype(np.float64))), variance_floor)
    order = np.argsort(values)
    sorted_values = values[order].astype(np.float64)

    for class_idx, n_gaussians in enumerate(gaussian_counts):
        component_indices = np.flatnonzero(component_classes == class_idx)
        weights = priors[:, class_idx].astype(np.float64)
        total = float(np.sum(weights))
        if total <= np.finfo(np.float32).eps:
            qs = np.linspace(0.2, 0.8, int(n_gaussians))
            class_means = np.quantile(values, qs)
            class_var = global_var
        else:
            sorted_weights = weights[order]
            qs = (np.arange(int(n_gaussians), dtype=np.float64) + 0.5) / float(n_gaussians)
            class_means = _weighted_quantiles_sorted(sorted_values, sorted_weights, qs)
            class_mean = float(np.sum(weights * values) / total)
            class_var = max(float(np.sum(weights * (values - class_mean) ** 2) / total), variance_floor)
        for idx, component_idx in enumerate(component_indices):
            means[component_idx] = float(class_means[idx])
            variances[component_idx] = class_var
            mixture_weights[component_idx] = 1.0 / float(n_gaussians)
    return means, variances, mixture_weights


def _weighted_quantiles_sorted(sorted_values: np.ndarray, sorted_weights: np.ndarray, quantiles: np.ndarray) -> np.ndarray:
    total = float(np.sum(sorted_weights))
    if total <= np.finfo(np.float32).eps:
        return np.quantile(sorted_values, quantiles)
    cdf = np.cumsum(sorted_weights)
    return np.interp(quantiles * total, cdf, sorted_values, left=sorted_values[0], right=sorted_values[-1])


def _gmm_maximization_step(
    values: np.ndarray,
    priors: np.ndarray,
    means: np.ndarray,
    variances: np.ndarray,
    mixture_weights: np.ndarray,
    component_classes: np.ndarray,
    gaussian_counts: np.ndarray,
    parameters: TpmGmmParameters,
    variance_floor: float,
    chunk_size: int = 500_000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_components = means.size
    weight_sum = np.zeros(n_components, dtype=np.float64)
    weighted_sum = np.zeros(n_components, dtype=np.float64)
    weighted_sumsq = np.zeros(n_components, dtype=np.float64)
    for start in range(0, values.size, chunk_size):
        stop = min(values.size, start + chunk_size)
        chunk_values = values[start:stop]
        resp = _component_responsibilities(
            chunk_values,
            priors[start:stop, :],
            means,
            variances,
            mixture_weights,
            component_classes,
            parameters,
        ).astype(np.float64, copy=False)
        chunk_values64 = chunk_values.astype(np.float64, copy=False)
        weight_sum += np.sum(resp, axis=0)
        weighted_sum += resp.T @ chunk_values64
        weighted_sumsq += resp.T @ (chunk_values64**2)

    new_means = means.copy()
    new_variances = variances.copy()
    for component_idx in range(n_components):
        if weight_sum[component_idx] <= np.finfo(np.float32).eps:
            continue
        mean = weighted_sum[component_idx] / weight_sum[component_idx]
        var = weighted_sumsq[component_idx] / weight_sum[component_idx] - mean**2
        new_means[component_idx] = float(mean)
        new_variances[component_idx] = max(float(var), variance_floor)

    new_mixture_weights = mixture_weights.copy()
    for class_idx, _ in enumerate(gaussian_counts):
        component_indices = np.flatnonzero(component_classes == class_idx)
        total = float(np.sum(weight_sum[component_indices]))
        if total <= np.finfo(np.float32).eps:
            new_mixture_weights[component_indices] = 1.0 / float(component_indices.size)
        else:
            new_mixture_weights[component_indices] = weight_sum[component_indices] / total
    return new_means, new_variances, new_mixture_weights


def _gmm_class_posteriors(
    values: np.ndarray,
    priors: np.ndarray,
    means: np.ndarray,
    variances: np.ndarray,
    mixture_weights: np.ndarray,
    component_classes: np.ndarray,
    n_classes: int,
    parameters: TpmGmmParameters,
    chunk_size: int = 500_000,
) -> np.ndarray:
    class_posteriors = np.zeros((values.size, n_classes), dtype=np.float32)
    for start in range(0, values.size, chunk_size):
        stop = min(values.size, start + chunk_size)
        resp = _component_responsibilities(
            values[start:stop],
            priors[start:stop, :],
            means,
            variances,
            mixture_weights,
            component_classes,
            parameters,
        )
        for class_idx in range(n_classes):
            class_posteriors[start:stop, class_idx] = np.sum(resp[:, component_classes == class_idx], axis=1)
    return class_posteriors


def _component_responsibilities(
    values: np.ndarray,
    priors: np.ndarray,
    means: np.ndarray,
    variances: np.ndarray,
    mixture_weights: np.ndarray,
    component_classes: np.ndarray,
    parameters: TpmGmmParameters,
) -> np.ndarray:
    log_resp = np.empty((values.size, means.size), dtype=np.float32)
    for component_idx, class_idx in enumerate(component_classes):
        var = max(float(variances[component_idx]), np.finfo(np.float32).eps)
        log_likelihood = -0.5 * ((values - means[component_idx]) ** 2 / var + np.log(2.0 * np.pi * var))
        log_prior = parameters.prior_weight * np.log(np.clip(priors[:, class_idx], parameters.prior_floor, 1.0))
        log_mix = np.log(max(float(mixture_weights[component_idx]), parameters.prior_floor))
        log_resp[:, component_idx] = (log_likelihood + log_prior + log_mix).astype(np.float32)
    return _softmax_last_axis(log_resp)


def initial_labels_from_probabilities(probabilities: np.ndarray, min_support: float = 0.05) -> np.ndarray:
    """Build CAT AMAP labels: 0=background, 1=CSF, 2=GM, 3=WM."""

    if probabilities.ndim != 4 or probabilities.shape[-1] != 3:
        raise ValueError(f"Expected GM/WM/CSF probability stack, got {probabilities.shape}")
    gm = probabilities[..., 0]
    wm = probabilities[..., 1]
    csf = probabilities[..., 2]
    stack = np.stack([csf, gm, wm], axis=-1)
    labels = np.argmax(stack, axis=-1).astype(np.uint8) + 1
    labels[np.sum(stack, axis=-1) <= min_support] = 0
    return labels


def _load_probability_stack(probabilities: dict[str, Path], labels: tuple[str, ...] = TISSUE_LABELS) -> np.ndarray:
    arrays = []
    ref_shape = None
    for label in labels:
        if label not in probabilities:
            raise ValueError(f"Missing {label} probability map.")
        data = np.nan_to_num(nib.load(str(probabilities[label])).get_fdata(dtype=np.float32).squeeze(), copy=False)
        if ref_shape is None:
            ref_shape = data.shape
        elif data.shape != ref_shape:
            raise ValueError("Probability maps do not share the same shape.")
        arrays.append(np.clip(data, 0.0, 1.0))
    return np.stack(arrays, axis=-1).astype(np.float32)


def _normalize_prior_stack(priors: np.ndarray, parameters: TpmGmmParameters) -> np.ndarray:
    prior_stack = np.clip(np.asarray(priors, dtype=np.float32), 0.0, 1.0)
    if prior_stack.ndim != 4 or prior_stack.shape[-1] < 1:
        raise ValueError(f"Expected 4D prior stack, got {prior_stack.shape}")
    support = np.sum(prior_stack, axis=-1) > parameters.support_threshold
    prior_stack = np.where(support[..., None], np.maximum(prior_stack, parameters.prior_floor), 0.0)
    denom = np.sum(prior_stack, axis=-1, keepdims=True)
    return np.divide(prior_stack, denom, out=np.zeros_like(prior_stack), where=denom > 0)


def _weighted_gaussian_stats(
    image: np.ndarray,
    responsibilities: np.ndarray,
    finite: np.ndarray,
    variance_floor: float,
) -> tuple[np.ndarray, np.ndarray]:
    means = np.zeros(3, dtype=np.float64)
    variances = np.zeros(3, dtype=np.float64)
    values = image[finite].astype(np.float64)
    for idx in range(3):
        weights = responsibilities[..., idx][finite].astype(np.float64)
        weight_sum = float(np.sum(weights))
        if weight_sum <= np.finfo(np.float32).eps:
            means[idx] = float(np.nanmedian(values))
            variances[idx] = variance_floor
            continue
        mean = float(np.sum(weights * values) / weight_sum)
        var = float(np.sum(weights * (values - mean) ** 2) / weight_sum)
        means[idx] = mean
        variances[idx] = max(var, variance_floor)
    return means, variances


def _softmax_last_axis(log_values: np.ndarray) -> np.ndarray:
    shifted = log_values - np.max(log_values, axis=-1, keepdims=True)
    exp_values = np.exp(shifted, dtype=np.float32)
    denom = np.sum(exp_values, axis=-1, keepdims=True)
    return np.divide(exp_values, denom, out=np.zeros_like(exp_values), where=denom > 0)


def _n4_bias_correct(t1_path: Path, mask_path: Path, output_path: Path, overwrite: bool, verbose: bool, threads: int) -> Path:
    if output_path.exists() and not overwrite:
        return output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    run_cli(
        [
            require_cli("N4BiasFieldCorrection"),
            "-d",
            "3",
            "-i",
            str(t1_path),
            "-x",
            str(mask_path),
            "-o",
            str(output_path),
        ],
        verbose=verbose,
        threads=threads,
    )
    return output_path


def _clip_probability_image(path: Path) -> None:
    img = nib.load(str(path))
    data = np.clip(np.nan_to_num(img.get_fdata(dtype=np.float32).squeeze(), copy=False), 0.0, 1.0)
    save_nifti(data.astype(np.float32), img, path, dtype=np.float32)


def _round_label_image(path: Path) -> None:
    img = nib.load(str(path))
    data = np.nan_to_num(img.get_fdata(dtype=np.float32).squeeze(), copy=False)
    save_nifti(np.clip(np.rint(data), 0, 255).astype(np.uint8), img, path, dtype=np.uint8)


def _nifti_like(reference: nib.Nifti1Image, data: np.ndarray) -> nib.Nifti1Image:
    header = reference.header.copy()
    header.set_data_dtype(np.float32)
    header.set_data_shape(data.shape)
    return nib.Nifti1Image(data.astype(np.float32), reference.affine, header)


def _subject_to_template_transform_prefix(output_dir: Path, subject: str, session: str | None) -> Path:
    sub = f"sub-{normalize_subject(subject)}"
    ses = f"ses-{normalize_session(session)}" if session else "ses-none"
    return output_dir / "transforms" / "cat12_tpm" / sub / ses / "anat" / f"{sub}_{ses}_desc-T1w_to_cat12Template"


def _template_prior_path(output_dir: Path, subject: str, session: str | None, label: str) -> Path:
    return anat_derivative(
        output_dir / "work" / "cat12_tpm",
        subject,
        session,
        space="cat12MNI",
        label=label,
        desc="templatePrior",
        suffix_override="probseg",
    )


def _subject_to_template_transform_set(prefix: Path) -> dict[str, list[str]]:
    forward = []
    inverse = []
    warp = prefix.with_suffix(".syn.nii.gz")
    inverse_warp = prefix.with_suffix(".syn_inv.nii.gz")
    affine = prefix.with_suffix(".affine.mat")
    if not affine.exists():
        return {"forward": forward, "inverse": inverse}
    if warp.exists() or inverse_warp.exists():
        if not (warp.exists() and inverse_warp.exists()):
            return {"forward": [], "inverse": []}
        forward.append(str(warp))
        inverse.append(f"[{affine},1]")
        inverse.append(str(inverse_warp))
    else:
        inverse.append(f"[{affine},1]")
    forward.append(str(affine))
    return {"forward": forward, "inverse": inverse}


def _run_ants_registration_cli(
    fixed: Path,
    moving: Path,
    out_prefix: Path,
    transform: str,
    threads: int,
    fixed_mask: str | Path | None,
    moving_mask: str | Path | None,
    verbose: bool,
) -> dict[str, list[str]]:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    ants_prefix = out_prefix.parent / f"{out_prefix.name}_"
    cmd = [
        require_cli("antsRegistrationSyN.sh"),
        "-d",
        "3",
        "-f",
        str(fixed),
        "-m",
        str(moving),
        "-o",
        str(ants_prefix),
        "-t",
        transform,
        "-n",
        str(max(1, int(threads))),
        "-y",
        "1",
        "-e",
        "1",
    ]
    if fixed_mask is not None or moving_mask is not None:
        fixed_mask_arg = str(fixed_mask) if fixed_mask is not None else "NULL"
        moving_mask_arg = str(moving_mask) if moving_mask is not None else "NULL"
        cmd.extend(["-x", f"{fixed_mask_arg},{moving_mask_arg}"])
    run_cli(cmd, verbose=verbose, threads=threads)

    produced_affine = Path(f"{ants_prefix}0GenericAffine.mat")
    produced_warp = Path(f"{ants_prefix}1Warp.nii.gz")
    produced_inverse_warp = Path(f"{ants_prefix}1InverseWarp.nii.gz")
    affine = out_prefix.with_suffix(".affine.mat")
    warp = out_prefix.with_suffix(".syn.nii.gz")
    inverse_warp = out_prefix.with_suffix(".syn_inv.nii.gz")
    if produced_warp.exists():
        shutil.copy2(produced_warp, warp)
    if not produced_affine.exists():
        raise RuntimeError(f"ANTs registration did not produce an affine transform: {produced_affine}")
    shutil.copy2(produced_affine, affine)
    if produced_inverse_warp.exists():
        shutil.copy2(produced_inverse_warp, inverse_warp)
    return _subject_to_template_transform_set(out_prefix)


def _register_subject_to_template(
    fixed: Path,
    moving: Path,
    out_prefix: Path,
    transform: str,
    threads: int,
    fixed_mask: str | Path | None,
    moving_mask: str | Path | None,
    verbose: bool,
) -> dict[str, list[str]]:
    reg = Registration()
    try:
        tx, _ = reg.register(
            fixed_input=fixed,
            moving_input=moving,
            fixed_mask=fixed_mask,
            moving_mask=moving_mask,
            transform=transform,
            verbose=verbose,
            threads=threads,
        )
        reg.save_all_transforms(tx, out_prefix)
        transforms = _subject_to_template_transform_set(out_prefix)
        if transforms["forward"]:
            return transforms
    except ANTsError:
        pass
    return _run_ants_registration_cli(
        fixed=fixed,
        moving=moving,
        out_prefix=out_prefix,
        transform=transform,
        threads=threads,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        verbose=verbose,
    )


def _warp_template_image_to_t1(
    fixed: Path,
    moving: Path,
    transforms: list[str],
    output: Path,
    interpolation: str,
    verbose: bool,
    threads: int,
) -> Path:
    return _apply_transforms_cli(
        fixed,
        moving,
        transforms,
        output,
        interpolation=interpolation,
        verbose=verbose,
        threads=threads,
    )


def _apply_transforms_cli(
    fixed: Path,
    moving: Path,
    transforms: list[str],
    output: Path,
    interpolation: str,
    verbose: bool,
    threads: int,
) -> Path:
    if not transforms:
        raise RuntimeError("No ANTs transforms were provided for TPM warping.")
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        require_cli("antsApplyTransforms"),
        "-d",
        "3",
        "-i",
        str(moving),
        "-r",
        str(fixed),
        "-o",
        str(output),
        "-n",
        _ants_cli_interpolation(interpolation),
    ]
    for transform in transforms:
        cmd.extend(["-t", transform])
    run_cli(cmd, verbose=verbose, threads=threads)
    return output


def _ants_cli_interpolation(interpolation: str) -> str:
    mapping = {
        "linear": "Linear",
        "nearest": "NearestNeighbor",
        "genericLabel": "GenericLabel",
        "bSpline": "BSpline",
    }
    return mapping.get(interpolation, interpolation)
