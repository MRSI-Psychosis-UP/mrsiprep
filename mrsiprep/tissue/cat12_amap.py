"""Python-facing CAT12 AMAP helpers.

The CAT12 C sources are vendored under ``mrsiprep/tissue/src/cat12_amap``.
This module defines the public Python contract for the future compiled
extension while keeping array ordering and class-order conventions explicit.
"""

from __future__ import annotations

import ctypes
import os
import platform
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CAT12_AMAP_SOURCE_DIR = Path(__file__).resolve().parent / "src" / "cat12_amap"


class CatAmapOptions(ctypes.Structure):
    """ctypes mirror of ``CatAmapOptions`` in ``cat_amap_core.h``."""

    _fields_ = [
        ("n_classes", ctypes.c_int),
        ("n_iters", ctypes.c_int),
        ("sub", ctypes.c_int),
        ("pve", ctypes.c_int),
        ("init_kmeans", ctypes.c_int),
        ("mrf_weight", ctypes.c_double),
        ("iters_icm", ctypes.c_int),
        ("bias_fwhm", ctypes.c_double),
        ("verbose", ctypes.c_int),
    ]


@dataclass(frozen=True)
class AmapParameters:
    """CAT12 AMAP/PVE parameters used by ``cat_main_amap.m``."""

    n_classes: int = 3
    n_iters: int = 10
    sub: int | None = None
    pve: int = 5
    init_kmeans: int = 0
    mrf_weight: float = 0.0
    iters_icm: int | None = None
    bias_fwhm: float = 0.0
    verbose: int = 0

    def resolved(self, voxel_size: tuple[float, float, float]) -> "AmapParameters":
        mean_vx = float(np.mean(voxel_size))
        sub = self.sub if self.sub is not None else max(1, int(round(64.0 / mean_vx)))
        iters_icm = self.iters_icm if self.iters_icm is not None else (50 if self.mrf_weight != 0 else 0)
        return AmapParameters(
            n_classes=self.n_classes,
            n_iters=self.n_iters,
            sub=sub,
            pve=self.pve,
            init_kmeans=self.init_kmeans,
            mrf_weight=self.mrf_weight,
            iters_icm=iters_icm,
            bias_fwhm=self.bias_fwhm,
            verbose=self.verbose,
        )


@dataclass(frozen=True)
class AmapResult:
    """Result returned by CAT12 AMAP in CAT/SPM tissue order."""

    probabilities: np.ndarray
    means: np.ndarray
    stds: np.ndarray
    bias_corrected: np.ndarray | None = None


def cat12_output_classes(n_classes: int, pve: int) -> int:
    if pve == 6:
        return n_classes + 3
    if pve == 5:
        return n_classes + 2
    return n_classes


def to_cat_tissue_order(probabilities: np.ndarray) -> np.ndarray:
    """Convert AMAP internal order ``CSF, GM, WM`` to CAT output ``GM, WM, CSF``."""

    if probabilities.ndim != 4 or probabilities.shape[-1] < 3:
        raise ValueError(f"Expected 4D probability array with at least 3 classes, got {probabilities.shape}")
    return probabilities[..., [1, 2, 0]]


def prepare_amap_arrays(src: np.ndarray, label: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate and convert AMAP inputs to CAT-compatible memory layout."""

    src_arr = np.asarray(src, dtype=np.float64)
    label_arr = np.asarray(label, dtype=np.uint8)
    if src_arr.ndim != 3:
        raise ValueError(f"AMAP source must be 3D, got shape {src_arr.shape}")
    if label_arr.shape != src_arr.shape:
        raise ValueError(f"AMAP label shape {label_arr.shape} does not match source shape {src_arr.shape}")
    if np.any(label_arr > 3):
        raise ValueError("AMAP labels must use 0=background, 1=CSF, 2=GM, 3=WM")
    return np.asfortranarray(src_arr), np.asfortranarray(label_arr)


def run_cat12_amap(
    src: np.ndarray,
    label: np.ndarray,
    voxel_size: tuple[float, float, float],
    parameters: AmapParameters | None = None,
) -> AmapResult:
    """Run CAT12 AMAP through a compiled extension or ctypes shared library."""

    src_arr, label_arr = prepare_amap_arrays(src, label)
    params = (parameters or AmapParameters()).resolved(voxel_size)
    try:
        from mrsiprep.tissue import _cat12_amap  # type: ignore
    except Exception as exc:
        result = _run_cat12_amap_ctypes(src_arr, label_arr, voxel_size, params)
    else:
        result = _cat12_amap.run(
            src_arr,
            label_arr,
            tuple(float(v) for v in voxel_size),
            {
                "n_classes": params.n_classes,
                "n_iters": params.n_iters,
                "sub": params.sub,
                "pve": params.pve,
                "init_kmeans": params.init_kmeans,
                "mrf_weight": params.mrf_weight,
                "iters_icm": params.iters_icm,
                "bias_fwhm": params.bias_fwhm,
                "verbose": params.verbose,
            },
        )

    probabilities = to_cat_tissue_order(np.asarray(result["probabilities"], dtype=np.float32) / 255.0)
    return AmapResult(
        probabilities=probabilities,
        means=np.asarray(result["means"], dtype=np.float64),
        stds=np.asarray(result["stds"], dtype=np.float64),
        bias_corrected=np.asarray(result["bias_corrected"], dtype=np.float64)
        if result.get("bias_corrected") is not None
        else None,
    )


def build_cat12_amap_library(build_dir: str | Path | None = None, force: bool = False) -> Path:
    """Build and return a shared library for the non-MEX CAT12 AMAP wrapper."""

    env_lib = os.environ.get("MRSIPREP_CAT12_AMAP_LIB")
    if env_lib:
        lib_path = Path(env_lib)
        if not lib_path.exists():
            raise RuntimeError(f"MRSIPREP_CAT12_AMAP_LIB does not exist: {lib_path}")
        return lib_path

    build_root = Path(
        build_dir
        or os.environ.get("MRSIPREP_CAT12_AMAP_BUILD_DIR")
        or Path(tempfile.gettempdir()) / "mrsiprep_cat12_amap"
    )
    build_root.mkdir(parents=True, exist_ok=True)
    lib_path = build_root / _shared_library_name()
    sources = _cat12_amap_sources()

    if not force and lib_path.exists():
        lib_mtime = lib_path.stat().st_mtime
        if all(path.stat().st_mtime <= lib_mtime for path in sources):
            return lib_path

    cmd = [
        "gcc",
        "-O3",
        "-std=c99",
        "-fPIC",
        "-shared",
        "-I",
        str(CAT12_AMAP_SOURCE_DIR),
        *[str(path) for path in sources],
        "-lm",
        "-o",
        str(lib_path),
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("gcc is required to build the CAT12 AMAP shared library.") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Failed to build CAT12 AMAP shared library:\n{exc.stderr}") from exc
    return lib_path


def _run_cat12_amap_ctypes(
    src: np.ndarray,
    label: np.ndarray,
    voxel_size: tuple[float, float, float],
    parameters: AmapParameters,
) -> dict[str, np.ndarray]:
    lib = _load_cat12_amap_library()
    out_classes = cat12_output_classes(parameters.n_classes, parameters.pve)
    prob = np.zeros((*src.shape, out_classes), dtype=np.uint8, order="F")
    means = np.zeros(parameters.n_classes, dtype=np.float64)
    stds = np.zeros(parameters.n_classes, dtype=np.float64)
    bias_corrected = np.zeros(src.shape, dtype=np.float64, order="F")

    dims = (ctypes.c_int * 3)(*src.shape)
    voxelsize = (ctypes.c_double * 3)(*tuple(float(v) for v in voxel_size))
    options = CatAmapOptions(
        parameters.n_classes,
        parameters.n_iters,
        int(parameters.sub or 1),
        parameters.pve,
        parameters.init_kmeans,
        parameters.mrf_weight,
        int(parameters.iters_icm or 0),
        parameters.bias_fwhm,
        parameters.verbose,
    )
    code = lib.cat_amap_run(
        src.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        label.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        dims,
        voxelsize,
        ctypes.byref(options),
        prob.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        means.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        stds.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        bias_corrected.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    if code != 0:
        raise RuntimeError(f"CAT12 AMAP failed: {_cat12_error_string(lib, code)}")
    return {
        "probabilities": prob,
        "means": means,
        "stds": stds,
        "bias_corrected": bias_corrected,
    }


def _load_cat12_amap_library():
    lib_path = build_cat12_amap_library()
    lib = ctypes.CDLL(str(lib_path))
    lib.cat_amap_run.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_ubyte),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(CatAmapOptions),
        ctypes.POINTER(ctypes.c_ubyte),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.cat_amap_run.restype = ctypes.c_int
    lib.cat_amap_error_string.argtypes = [ctypes.c_int]
    lib.cat_amap_error_string.restype = ctypes.c_char_p
    return lib


def _cat12_error_string(lib, code: int) -> str:
    raw = lib.cat_amap_error_string(code)
    return raw.decode("utf-8", errors="replace") if raw else f"unknown error {code}"


def _cat12_amap_sources() -> list[Path]:
    return [
        CAT12_AMAP_SOURCE_DIR / "cat_amap_core.c",
        CAT12_AMAP_SOURCE_DIR / "Amap.c",
        CAT12_AMAP_SOURCE_DIR / "Kmeans.c",
        CAT12_AMAP_SOURCE_DIR / "MrfPrior.c",
        CAT12_AMAP_SOURCE_DIR / "Pve.c",
        CAT12_AMAP_SOURCE_DIR / "vollib.c",
    ]


def _shared_library_name() -> str:
    system = platform.system().lower()
    if system == "darwin":
        return "libcat_amap_core.dylib"
    if system == "windows":
        return "cat_amap_core.dll"
    return "libcat_amap_core.so"
