"""Registration QC placeholder figures."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from mrsiprep.io.naming import figure_derivative


def write_registration_qc(config, subject: str, session: str | None, fixed: Path, moving: Path) -> Path:
    import matplotlib.pyplot as plt

    out = figure_derivative(config.derivative_dir, subject, session, extension="svg", desc="registrationqc")
    if out.exists() and not config.overwrite:
        return out
    fixed_data = nib.load(str(fixed)).get_fdata()
    moving_data = nib.load(str(moving)).get_fdata()
    zf = fixed_data.shape[2] // 2
    zm = moving_data.shape[2] // 2
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    axes[0].imshow(np.rot90(fixed_data[:, :, zf]), cmap="gray")
    axes[0].set_title("T1 target")
    axes[1].imshow(np.rot90(moving_data[:, :, zm]), cmap="magma")
    axes[1].set_title("MRSI reference")
    for ax in axes:
        ax.axis("off")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out
