"""Provenance helpers."""

from __future__ import annotations

import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
from rich import box
from rich.table import Table

from mrsiprep.utils.debug import Debug

try:
    from mrsiprep import __version__
except ImportError:
    try:
        __version__ = version("mrsiprep")
    except PackageNotFoundError:
        __version__ = "unknown"


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


def check_external_software(debug: Debug) -> bool:
    statuses = software_versions()
    table = Table(box=box.SIMPLE_HEAVY, show_lines=False, title="External software availability")
    table.add_column("Tool", style="cyan", no_wrap=True)
    table.add_column("Path", style="white")
    table.add_column("Status", justify="center")

    all_available = True
    for tool, path in statuses.items():
        if path:
            table.add_row(tool, path, "[green]INSTALLED[/green]")
        else:
            table.add_row(tool, "[red]missing[/red]", "[bold red]MISSING[/bold red]")
            all_available = False

    debug.separator()
    debug.title("External software availability")
    debug.console.print(table)
    if all_available:
        debug.success("All required external binaries are available.")
    else:
        debug.failure("One or more required external binaries are missing.")
    return all_available


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
