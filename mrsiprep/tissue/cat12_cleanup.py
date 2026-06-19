"""CAT12-style post-AMAP cleanup helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi


TISSUE_LABELS = ("GM", "WM", "CSF")


class CATLabels:
    """Major-region IDs from CAT12 ``cat_defaults.m`` for ``cat.nii``."""

    CT = 1
    CB = 3
    BG = 5
    BV = 7
    TH = 9
    ON = 11
    MB = 13
    BS = 13
    VT = 15
    NV = 17
    HC = 19
    HD = 21
    HI = 23
    PH = 25
    LE = 27


@dataclass(frozen=True)
class Cat12CleanupParameters:
    """Controls mirroring CAT12's post-AMAP ``cleanupstr`` handling."""

    enabled: bool = True
    cleanup_strength: float = 0.5
    extra_cleanup: int = 0
    atlas_cleanup: bool = False
    outer_csf_correction: bool = True
    outer_csf_distance_mm: float = 3.0
    outer_csf_intensity: float = 2.55

    def level(self, voxel_size: tuple[float, float, float]) -> float:
        """Return CAT's ``cat_main_clean_gwc`` level for this voxel size."""

        mean_vx = max(float(np.mean(voxel_size)), np.finfo(np.float32).eps)
        return float(min(1.0, max(0.0, self.cleanup_strength * 2.0 / mean_vx)))

    def final_cleanup_strength(self, voxel_size: tuple[float, float, float]) -> float:
        """Return CAT's internal final-cleanup strength."""

        mean_vx = max(1.0, float(np.mean(voxel_size)))
        return float(min(1.0, max(0.0, self.cleanup_strength * 0.5 / mean_vx)))

    def final_cleanup_distance(self) -> float:
        return float(min(2.0, max(0.0, 1.0 + self.cleanup_strength)))


def clean_gwc(
    probabilities: np.ndarray,
    voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0),
    parameters: Cat12CleanupParameters | None = None,
) -> np.ndarray:
    """Approximate CAT12 ``cat_main_clean_gwc`` for GM/WM/CSF probabilities.

    CAT passes AMAP probabilities through ``cat_main_clean_gwc`` before
    writing p1/p2/p3 maps.  The routine uses WM-connected conditional
    dilations to limit GM/WM and then grows a second mask for CSF.  This
    Python port keeps the same class order and threshold logic, while using
    SciPy morphology/convolution in place of SPM/CAT MATLAB primitives.
    """

    params = parameters or Cat12CleanupParameters()
    probs = _validate_probability_stack(probabilities)
    if not params.enabled:
        return probs.copy()

    level = params.level(voxel_size)
    p_uint8 = _to_uint8_probabilities(probs)
    if params.extra_cleanup in (1, 3):
        p_uint8 = _clean_close_to_skull(p_uint8, level, voxel_size)
    if params.extra_cleanup in (2, 3):
        p_uint8 = _clean_unconnected_wm_gm(p_uint8, level, voxel_size)
    p_uint8 = _classic_clean_gwc(p_uint8, level)
    return np.clip(p_uint8.astype(np.float32) / 255.0, 0.0, 1.0)


def clean_final_with_atlas(
    probabilities: np.ndarray,
    normalized_t1: np.ndarray,
    atlas_labels: np.ndarray,
    voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0),
    parameters: Cat12CleanupParameters | None = None,
) -> np.ndarray:
    """Approximate CAT12 ``cat_main_cleanup`` using the native CAT atlas.

    This is the cleanup stage CAT applies after ``cat_main_clean_gwc``.  It
    uses the CAT major-region atlas to protect ventricles/cerebellum/deep
    structures and to convert likely meningeal GM/WM to CSF or background.
    The full MATLAB routine depends on CAT's deformation field and several
    atlas-refinement products; this port implements the three direct
    probability-map corrections in native space.
    """

    params = parameters or Cat12CleanupParameters()
    probs = _validate_probability_stack(probabilities)
    if not (params.enabled and params.atlas_cleanup):
        return probs.copy()
    intensity = np.nan_to_num(np.asarray(normalized_t1, dtype=np.float32).squeeze(), copy=False)
    atlas = np.nan_to_num(np.asarray(atlas_labels, dtype=np.float32).squeeze(), copy=False)
    if intensity.shape != probs.shape[:3]:
        raise ValueError(f"Normalized T1 shape {intensity.shape} does not match probabilities {probs.shape[:3]}.")
    if atlas.shape != probs.shape[:3]:
        raise ValueError(f"Atlas shape {atlas.shape} does not match probabilities {probs.shape[:3]}.")

    p_uint8 = _to_uint8_probabilities(probs)
    support = np.sum(p_uint8, axis=-1) > 0
    p_crop, bbox = _crop_to_support(p_uint8, support)
    if p_crop.size == 0 or bbox is None:
        return probs.copy()
    intensity_crop = intensity[bbox].astype(np.float32, copy=False)
    atlas_crop = np.clip(np.rint(atlas[bbox]), 0, 255).astype(np.uint8)
    atlas_crop = _fill_atlas_holes(atlas_crop, np.sum(p_crop, axis=-1) > 0)

    cleaned = _cat_main_cleanup_uint8(p_crop, intensity_crop, atlas_crop, voxel_size, params)
    return np.clip(_uncrop(cleaned, p_uint8.shape, bbox).astype(np.float32) / 255.0, 0.0, 1.0)


def correct_outer_csf_with_atlas(
    probabilities: np.ndarray,
    normalized_t1: np.ndarray,
    atlas_labels: np.ndarray,
    voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0),
    parameters: Cat12CleanupParameters | None = None,
) -> np.ndarray:
    """Correct WM-dominant outer-boundary voxels that are CSF-like."""

    params = parameters or Cat12CleanupParameters()
    probs = _validate_probability_stack(probabilities)
    if not (params.enabled and params.outer_csf_correction):
        return probs.copy()
    intensity = np.nan_to_num(np.asarray(normalized_t1, dtype=np.float32).squeeze(), copy=False)
    atlas = np.nan_to_num(np.asarray(atlas_labels, dtype=np.float32).squeeze(), copy=False)
    if intensity.shape != probs.shape[:3]:
        raise ValueError(f"Normalized T1 shape {intensity.shape} does not match probabilities {probs.shape[:3]}.")
    if atlas.shape != probs.shape[:3]:
        raise ValueError(f"Atlas shape {atlas.shape} does not match probabilities {probs.shape[:3]}.")

    p_uint8 = _to_uint8_probabilities(probs)
    support = np.sum(p_uint8, axis=-1) > 0
    p_crop, bbox = _crop_to_support(p_uint8, support)
    if p_crop.size == 0 or bbox is None:
        return probs.copy()
    intensity_crop = intensity[bbox].astype(np.float32, copy=False)
    atlas_crop = np.clip(np.rint(atlas[bbox]), 0, 255).astype(np.uint8)
    atlas_crop = _fill_atlas_holes(atlas_crop, np.sum(p_crop, axis=-1) > 0)
    corrected = _correct_outer_csf_boundary_uint8(p_crop, intensity_crop, atlas_crop, voxel_size, params)
    return np.clip(_uncrop(corrected, p_uint8.shape, bbox).astype(np.float32) / 255.0, 0.0, 1.0)


def _validate_probability_stack(probabilities: np.ndarray) -> np.ndarray:
    probs = np.asarray(probabilities, dtype=np.float32)
    if probs.ndim != 4 or probs.shape[-1] != 3:
        raise ValueError(f"Expected a GM/WM/CSF probability stack, got {probs.shape}")
    return np.clip(np.nan_to_num(probs, copy=False), 0.0, 1.0)


def _to_uint8_probabilities(probabilities: np.ndarray) -> np.ndarray:
    return np.clip(np.rint(probabilities * 255.0), 0, 255).astype(np.uint8)


def _classic_clean_gwc(p_uint8: np.ndarray, level: float) -> np.ndarray:
    p = p_uint8.astype(np.float32, copy=True)
    p, bbox = _crop_to_support(p, np.sum(p, axis=-1) > 0)
    if p.size == 0:
        return p_uint8.copy()

    gm = p[..., 0]
    wm = p[..., 1]
    csf = p[..., 2]

    b = wm.copy()
    b *= _largest_components(b > 128, max_components=10, relative_threshold=0.1).astype(np.float32)

    kernel = _cat_cleanup_kernel()
    th1 = 0.2 if level > 1.0 else 0.15
    vxa = float(np.mean(p.shape[:3]) / 256.0)
    niter = max(1, int(np.floor(32.0 * vxa)))
    niter2 = niter

    for j in range(1, niter + 1):
        th = th1 if j > (2.0 * vxa) else 0.6
        b = ((b / 255.0) > th).astype(np.float32) * (wm + gm)
        b = ndi.convolve(np.rint(b), kernel, mode="nearest")

    c = b.copy()
    for _ in range(niter2):
        c = ((c / 255.0) > th).astype(np.float32) * (wm + gm + csf)
        c = ndi.convolve(np.rint(c), kernel, mode="nearest")

    th = 0.05
    gm_f = gm / 255.0
    wm_f = wm / 255.0
    csf_f = csf / 255.0
    brain = (((b / 255.0) > th) * (gm_f + wm_f)) > th
    gm_f *= brain
    wm_f *= brain
    csf_support = (((c / 255.0) > th) * (gm_f + wm_f + csf_f)) > th
    csf_f *= csf_support

    cleaned = _renormalize_probabilities(np.stack([gm_f, wm_f, csf_f], axis=-1))
    return _uncrop(_to_uint8_probabilities(cleaned), p_uint8.shape, bbox)


def _clean_close_to_skull(p_uint8: np.ndarray, level: float, voxel_size: tuple[float, float, float]) -> np.ndarray:
    p = p_uint8.astype(np.float32, copy=True)
    y_p0 = p[..., 2] / 255.0 + p[..., 0] / 255.0 * 2.0 + p[..., 1] / 255.0 * 3.0
    y_be = _morph(_morph(_morph(y_p0 > 0.5, "ldc", 1.0, voxel_size), "de", 5.0 * level, voxel_size), "lc", 1.0, voxel_size)

    y_msk = (
        y_be
        | _morph(_morph(y_p0 > 2.5, "do", 0.5 + level / 2.0, voxel_size), "l", (10, 0.1), voxel_size)
        | _morph(_morph(y_p0 > 2.8, "do", level, voxel_size), "l", (10, 0.1), voxel_size)
    )
    p[..., 1] *= ndi.gaussian_filter(_morph(y_msk, "dd", 2.5, voxel_size).astype(np.float32), sigma=0.5)

    y_p0 = p[..., 2] / 255.0 + p[..., 0] / 255.0 * 2.0 + p[..., 1] / 255.0 * 3.0
    y_wd = _morph(y_p0 > 2.5, "dd", 5.0 - level, voxel_size)
    for _ in range(2):
        y_msk = _morph(y_be | (y_wd & _morph(y_p0 > 1.5, "do", level / 2.0 + 1.5, voxel_size)), "ldc", 1.0, voxel_size)
        p[..., 0] *= ndi.gaussian_filter(_morph(y_msk, "dd", 1.5, voxel_size).astype(np.float32), sigma=0.5)
    y_p0 = p[..., 2] / 255.0 + p[..., 0] / 255.0 * 2.0 + p[..., 1] / 255.0 * 3.0
    for _ in range(2):
        y_msk = _morph(y_be | _morph(y_p0 > 0.95, "do", level + 2.5, voxel_size), "ldc", 1.0, voxel_size)
        p[..., 2] *= ndi.gaussian_filter(_morph(y_msk, "dd", 1.5, voxel_size).astype(np.float32), sigma=0.5)
    return np.clip(np.rint(p), 0, 255).astype(np.uint8)


def _clean_unconnected_wm_gm(p_uint8: np.ndarray, level: float, voxel_size: tuple[float, float, float]) -> np.ndarray:
    p = p_uint8.astype(np.float32, copy=True)
    for _ in range(2):
        wm = p[..., 1] / 255.0
        wm_unconnected = (wm > 0.5) & ~_largest_components(wm > 0.5, max_components=10, relative_threshold=0.1)
        near_csf = _morph(p[..., 2] > 0, "dd", round(level), voxel_size)
        is_gm = ndi.uniform_filter(p[..., 0], size=3, mode="nearest") > ndi.uniform_filter(p[..., 2], size=3, mode="nearest")
        no_wm = (ndi.uniform_filter(wm, size=3, mode="nearest") < 0.5) & near_csf | wm_unconnected
        move = p[..., 1] * (no_wm & is_gm)
        p[..., 0] += move
        p[..., 1] -= move
        move = p[..., 1] * (no_wm & ~is_gm)
        p[..., 2] += move
        p[..., 1] -= move

    for _ in range(2):
        gwm = (p[..., 0] + p[..., 1]) / 255.0
        gwm_unconnected = (gwm > 0.5) & ~_largest_components(gwm > 0.5, max_components=10, relative_threshold=0.1)
        near_csf = _morph(p[..., 2] > 0, "dd", round(2.0 * level), voxel_size)
        is_brain_mask = ndi.uniform_filter(p[..., 2], size=3, mode="nearest") > 127.5
        no_wm = (ndi.uniform_filter(gwm, size=3, mode="nearest") < 0.5) & near_csf | gwm_unconnected
        move = p[..., 0] * (no_wm & is_brain_mask)
        p[..., 2] += move
        p[..., 0] -= move
    return np.clip(np.rint(p), 0, 255).astype(np.uint8)


def _cat_p0_from_uint8(probabilities: np.ndarray) -> np.ndarray:
    p = probabilities.astype(np.float32, copy=False)
    return p[..., 2] / 255.0 + p[..., 0] / 255.0 * 2.0 + p[..., 1] / 255.0 * 3.0


def _ns(labels: np.ndarray, label: int) -> np.ndarray:
    return (labels == label) | (labels == label + 1)


def _distance_to_mask(mask: np.ndarray, voxel_size: tuple[float, float, float]) -> np.ndarray:
    source = np.asarray(mask) > 0
    if not np.any(source):
        return np.full(source.shape, np.inf, dtype=np.float32)
    return ndi.distance_transform_edt(~source, sampling=voxel_size).astype(np.float32)


def _smooth3(values: np.ndarray) -> np.ndarray:
    return ndi.uniform_filter(values.astype(np.float32, copy=False), size=3, mode="nearest")


def _gaussian_smooth(values: np.ndarray, fwhm_mm: float, voxel_size: tuple[float, float, float]) -> np.ndarray:
    voxel = np.asarray(voxel_size, dtype=np.float32)
    sigma = float(fwhm_mm) / np.maximum(voxel, np.finfo(np.float32).eps) / np.sqrt(8.0 * np.log(2.0))
    return ndi.gaussian_filter(values.astype(np.float32, copy=False), sigma=sigma, mode="nearest")


def _fill_atlas_holes(labels: np.ndarray, support: np.ndarray) -> np.ndarray:
    known = labels > 0
    if not np.any(known):
        return labels
    _, indices = ndi.distance_transform_edt(~known, return_indices=True)
    filled = labels[tuple(indices)]
    return np.where((labels > 0) | ~support, labels, filled).astype(np.uint8)


def _cat_main_cleanup_uint8(
    p_uint8: np.ndarray,
    normalized_t1: np.ndarray,
    atlas_labels: np.ndarray,
    voxel_size: tuple[float, float, float],
    parameters: Cat12CleanupParameters,
) -> np.ndarray:
    p = p_uint8.astype(np.float32, copy=True)
    ymb = np.asarray(normalized_t1, dtype=np.float32)
    labels = np.asarray(atlas_labels, dtype=np.uint8)
    cleanupstr = parameters.final_cleanup_strength(voxel_size)
    cleanupdist = parameters.final_cleanup_distance()
    max_vx = max(float(np.max(voxel_size)), np.finfo(np.float32).eps)
    vxv = 1.0 / max_vx

    outer_csf_seed = (
        _outer_csf_boundary_mask(p, ymb, labels, voxel_size, parameters)
        if parameters.outer_csf_correction
        else np.zeros(p.shape[:3], dtype=bool)
    )
    if np.any(outer_csf_seed):
        p = _apply_outer_csf_mask(p, outer_csf_seed).astype(np.float32)
    y_p0 = _cat_p0_from_uint8(p)
    ybd = _distance_to_mask(~_morph(y_p0 > 0, "ldc", vxv, voxel_size), voxel_size)
    ybd = _distance_to_mask(~_morph((y_p0 > 1.5) | (ybd > 8.0), "ldc", vxv, voxel_size), voxel_size)

    yvt = _morph(_ns(labels, CATLabels.VT) | _ns(labels, CATLabels.BG), "dd", vxv * 3.0, voxel_size)
    ycbp = _morph(_ns(labels, CATLabels.CB), "dd", cleanupdist * vxv, voxel_size)
    ycbn = _ns(labels, CATLabels.CB)
    ylhp = _morph((labels == CATLabels.CT) & (y_p0 < 2.1), "dd", cleanupdist * vxv * 2.0, voxel_size)
    yrhp = _morph((labels == CATLabels.CT + 1) & (y_p0 < 2.1), "dd", cleanupdist * vxv * 2.0, voxel_size)
    yroi = (
        (ybd < cleanupdist * 2.0)
        | (~ycbn & ycbp & (ylhp | yrhp))
        | (ylhp & yrhp)
        | _ns(labels, CATLabels.VT)
        | (_ns(labels, CATLabels.BS) & ycbp)
    )
    yrbv = (y_p0 > 0) & (ybd < 6.0) & _morph((ylhp & yrhp) | (~ycbn & ycbp & (ylhp | yrhp)), "dd", 4.0, voxel_size)
    yroi = (yroi | yrbv) & ~_ns(labels, CATLabels.BS) & ~ycbn

    yrw = (y_p0 > 0) & yroi & (ymb > (1.1 + ybd / 20.0)) & ~_ns(labels, CATLabels.CB)
    yrw = yrw | (_smooth3(yrw.astype(np.float32)) > (0.4 - 0.3 * cleanupstr))
    ygw = _morph((y_p0 >= 1.9) & ~yrw, "ldo", 0.0, voxel_size)
    yrw = yrw | ((y_p0 > 1.0) & yroi & ~ygw)
    yrw = yrw & ~yvt & ~_morph(ygw, "dd", 1.5, voxel_size)
    yrw[_smooth3(yrw.astype(np.float32)) < (0.5 + 0.2 * cleanupstr)] = False
    yrw[_smooth3(yrw.astype(np.float32)) < (0.5 - 0.2 * cleanupstr)] = False

    ybb = _morph(((y_p0 > 0) & ~yrw) | (ybd > 2.0), "ldo", 2.0 / vxv, voxel_size)
    ybb[(_gaussian_smooth(ybb.astype(np.float32), 2.0, voxel_size) > 0.4) & ~yrw] = True
    ybb = _morph(ybb | (ybd > 3.0), "ldc", 1.0 / vxv, voxel_size)
    ybb = _gaussian_smooth(ybb.astype(np.float32), 0.6, voxel_size) > (1.0 / 3.0)

    p[..., 0] = np.minimum(p[..., 0], ybb.astype(np.float32) * 255.0)
    p[..., 1] = np.minimum(p[..., 1], ybb.astype(np.float32) * 255.0)
    p[..., 2] = np.minimum(p[..., 2], ybb.astype(np.float32) * 255.0)
    p[..., 0] = np.minimum(p[..., 0], (~(ybb & yrw)).astype(np.float32) * 255.0)
    p[..., 1] = np.minimum(p[..., 1], (~(ybb & yrw)).astype(np.float32) * 255.0)
    p[..., 2] = np.maximum(p[..., 2], (ybb & yrw).astype(np.float32) * 255.0)

    y_p0 = _cat_p0_from_uint8(p)
    ym = _largest_components(
        ((p[..., 0] + p[..., 1]) > (160.0 + 32.0 * cleanupstr))
        & ~_morph((y_p0 > 1.0) & (y_p0 < 1.5 + cleanupstr / 2.0), "do", vxv, voxel_size)
    ).astype(np.float32)
    ym2 = _morph(ym, "do", min(1.0, 0.7 / max_vx), voxel_size)
    ym[_ns(labels, CATLabels.CT) & ~ym2] = 0.0
    ym = _gaussian_smooth(ym, 0.6, voxel_size)
    ym = (ym < (0.1 * cleanupstr)) & ybb & ~yvt & (ymb > 0.25)
    p[..., 0] = np.minimum(p[..., 0], (~ym).astype(np.float32) * 255.0)
    p[..., 1] = np.minimum(p[..., 1], (~ym).astype(np.float32) * 255.0)
    p[..., 2] = np.maximum(p[..., 2], (ym | (ybb & (y_p0 == 0))).astype(np.float32) * 255.0)

    y_p0 = _cat_p0_from_uint8(p)
    ybs = _ns(labels, CATLabels.BS) & (ymb > (2.0 / 3.0))
    ypve_vb = _morph(_ns(labels, CATLabels.VT) | ybs, "dd", 2.0, voxel_size)
    ypve_cc = (
        _morph(labels == CATLabels.CT, "dd", 3.0 * vxv, voxel_size)
        & _morph(labels == CATLabels.CT + 1, "dd", 3.0 * vxv, voxel_size)
        & _morph(_ns(labels, CATLabels.VT), "dd", 2.0, voxel_size)
    )
    ynpve = _smooth3((_ns(labels, CATLabels.BG) | _ns(labels, CATLabels.TH)).astype(np.float32)) > 0.3
    mid_pve = (y_p0 < 3.0) & (y_p0 > 1.0)
    yroi = (
        (ypve_vb | ypve_cc)
        & ~ynpve
        & _morph(y_p0 >= 2.95, "dd", 2.0 * vxv, voxel_size)
        & _morph(y_p0 <= 1.05, "dd", 2.0 * vxv, voxel_size)
        & mid_pve
        & (_smooth3((mid_pve & ~_morph(mid_pve, "do", 1.5 * vxv, voxel_size)).astype(np.float32)) > 0.1)
    )
    yncm = np.clip((3.0 - y_p0) / 2.0, 0.0, 1.0) * yroi.astype(np.float32)
    p[..., 0] = np.minimum(p[..., 0], (~yroi).astype(np.float32) * 255.0)
    p[..., 1] = p[..., 1] * (~yroi).astype(np.float32) + (yroi.astype(np.float32) - yncm) * 255.0
    p[..., 2] = p[..., 2] * (~yroi).astype(np.float32) + yncm * 255.0
    if np.any(outer_csf_seed):
        p = _apply_outer_csf_mask(p, outer_csf_seed)
    elif parameters.outer_csf_correction:
        p = _correct_outer_csf_boundary_uint8(p, ymb, labels, voxel_size, parameters)
    return np.clip(np.rint(p), 0, 255).astype(np.uint8)


def _correct_outer_csf_boundary_uint8(
    probabilities: np.ndarray,
    normalized_t1: np.ndarray,
    atlas_labels: np.ndarray,
    voxel_size: tuple[float, float, float],
    parameters: Cat12CleanupParameters,
) -> np.ndarray:
    p = probabilities.astype(np.float32, copy=True)
    candidate = _outer_csf_boundary_mask(p, normalized_t1, atlas_labels, voxel_size, parameters)
    if not np.any(candidate):
        return np.clip(np.rint(p), 0, 255).astype(np.uint8)
    return _apply_outer_csf_mask(p, candidate)


def _outer_csf_boundary_mask(
    probabilities: np.ndarray,
    normalized_t1: np.ndarray,
    atlas_labels: np.ndarray,
    voxel_size: tuple[float, float, float],
    parameters: Cat12CleanupParameters,
) -> np.ndarray:
    p = probabilities.astype(np.float32, copy=False)
    ymb = np.asarray(normalized_t1, dtype=np.float32)
    labels = np.asarray(atlas_labels, dtype=np.uint8)
    y_p0 = _cat_p0_from_uint8(p)
    support = y_p0 > 0
    if not np.any(support):
        return np.zeros(p.shape[:3], dtype=bool)

    distance_inside = ndi.distance_transform_edt(support, sampling=voxel_size).astype(np.float32)
    outer_edge = support & (distance_inside <= parameters.outer_csf_distance_mm)
    protected = (
        _ns(labels, CATLabels.VT)
        | _ns(labels, CATLabels.BG)
        | _ns(labels, CATLabels.TH)
        | _ns(labels, CATLabels.HC)
        | _ns(labels, CATLabels.BS)
        | _ns(labels, CATLabels.CB)
    )
    wm_dominant = (p[..., 1] > 96.0) & (p[..., 1] >= p[..., 0]) & (p[..., 1] > p[..., 2])
    low_for_wm = ymb < parameters.outer_csf_intensity
    local_csf_or_background = _morph((p[..., 2] > 32.0) | ~support, "dd", 1.5, voxel_size)
    candidate = outer_edge & wm_dominant & low_for_wm & local_csf_or_background & ~protected
    candidate &= _smooth3(candidate.astype(np.float32)) > 0.05
    return candidate


def _apply_outer_csf_mask(probabilities: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    p = probabilities.astype(np.float32, copy=True)
    p[..., 0] = np.where(candidate, 0.0, p[..., 0])
    p[..., 1] = np.where(candidate, 0.0, p[..., 1])
    p[..., 2] = np.where(candidate, 255.0, p[..., 2])
    return np.clip(np.rint(p), 0, 255).astype(np.uint8)


def _renormalize_probabilities(probabilities: np.ndarray) -> np.ndarray:
    denom = np.sum(probabilities, axis=-1, keepdims=True) + np.finfo(np.float32).eps
    return np.divide(probabilities, denom, out=np.zeros_like(probabilities, dtype=np.float32), where=denom > 0)


def _cat_cleanup_kernel() -> np.ndarray:
    k = np.asarray([0.75, 1.0, 0.75], dtype=np.float32)
    sm = float(np.sum(np.kron(np.kron(k, k), k)) ** (1.0 / 3.0))
    k = k / sm
    return np.einsum("i,j,k->ijk", k, k, k).astype(np.float32)


def _crop_to_support(data: np.ndarray, support: np.ndarray) -> tuple[np.ndarray, tuple[slice, slice, slice] | None]:
    coords = np.argwhere(support)
    if coords.size == 0:
        return data[:0, :0, :0].copy(), None
    lo = np.maximum(coords.min(axis=0) - 2, 0)
    hi = np.minimum(coords.max(axis=0) + 3, support.shape)
    bbox = tuple(slice(int(l), int(h)) for l, h in zip(lo, hi))
    return data[bbox].copy(), bbox


def _uncrop(data: np.ndarray, shape: tuple[int, ...], bbox: tuple[slice, slice, slice] | None) -> np.ndarray:
    out = np.zeros(shape, dtype=data.dtype)
    if bbox is not None:
        out[bbox] = data
    return out


def _largest_components(mask: np.ndarray, max_components: int = 1, relative_threshold: float = 0.0) -> np.ndarray:
    labels, n_labels = ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))
    if n_labels == 0:
        return np.zeros(mask.shape, dtype=bool)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    order = np.argsort(counts)[::-1]
    total = float(np.sum(counts))
    keep = np.zeros(n_labels + 1, dtype=bool)
    kept = 0
    for label in order:
        if label == 0 or counts[label] == 0 or kept >= max_components:
            break
        if relative_threshold > 0.0 and (counts[label] / max(total, 1.0)) <= relative_threshold:
            continue
        keep[label] = True
        kept += 1
    return keep[labels]


def _morph(mask: np.ndarray, action: str, radius: float | int | tuple[int, float] = 1.0, voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> np.ndarray:
    binary = np.asarray(mask) > 0.5
    if action in {"l", "lab"}:
        if isinstance(radius, tuple):
            return _largest_components(binary, max_components=int(radius[0]), relative_threshold=float(radius[1]))
        return _largest_components(binary, max_components=int(radius))
    if action in {"dd", "distdilate"}:
        return ndi.binary_dilation(binary, structure=_ball_structure(float(radius), voxel_size))
    if action in {"de", "disterode"}:
        return ndi.binary_erosion(binary, structure=_ball_structure(float(radius), voxel_size), border_value=1)
    if action in {"dc", "distclose"}:
        return ndi.binary_closing(binary, structure=_ball_structure(float(radius), voxel_size))
    if action in {"do", "distopen"}:
        return ndi.binary_opening(binary, structure=_ball_structure(float(radius), voxel_size))
    if action in {"lc", "labclose", "ldc", "labdistclose"}:
        return _largest_components(_morph(binary, "dc", float(radius), voxel_size))
    if action in {"lo", "labopen", "ldo", "labdistopen"}:
        return _largest_components(_morph(binary, "do", float(radius), voxel_size))
    raise ValueError(f"Unsupported CAT morphology action: {action}")


def _ball_structure(radius_mm: float, voxel_size: tuple[float, float, float]) -> np.ndarray:
    radius_mm = max(float(radius_mm), 0.0)
    if radius_mm <= 0.0:
        return np.ones((1, 1, 1), dtype=bool)
    voxel = np.asarray(voxel_size, dtype=np.float32)
    extents = np.maximum(1, np.ceil(radius_mm / np.maximum(voxel, np.finfo(np.float32).eps)).astype(int))
    grids = np.ogrid[
        -extents[0] : extents[0] + 1,
        -extents[1] : extents[1] + 1,
        -extents[2] : extents[2] + 1,
    ]
    dist = (grids[0] * voxel[0]) ** 2 + (grids[1] * voxel[1]) ** 2 + (grids[2] * voxel[2]) ** 2
    return dist <= (radius_mm + 0.5 * float(np.min(voxel))) ** 2
