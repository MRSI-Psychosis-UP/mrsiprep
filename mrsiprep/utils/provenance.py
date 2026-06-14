"""Provenance helpers."""

from __future__ import annotations

import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from mrsiprep import __version__


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):  # noqa: D401
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, Path):
            return str(obj)
        return super().default(obj)


def software_versions() -> dict[str, str | None]:
    tools = ["antsRegistrationSyN.sh", "Atropos", "N4BiasFieldCorrection", "hd-bet", "fslmaths", "fast", "petpvc", "chimera", "recon-all"]
    return {tool: shutil.which(tool) for tool in tools}


def write_provenance(config, out_path: str | Path, extra: dict | None = None) -> Path:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mrsiprep_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "config": config.to_dict() if hasattr(config, "to_dict") else {},
        "software": software_versions(),
    }
    if extra:
        payload.update(extra)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, cls=NumpyEncoder), encoding="utf-8")
    return out_path
