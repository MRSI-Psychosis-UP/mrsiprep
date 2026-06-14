"""Connectivity workflow."""

from __future__ import annotations

from mrsiprep.connectivity.export import export_connectivity


def run_connectivity_workflow(config, subject, session, regional_table, parcels):
    if not config.write_connectivity:
        return {}
    return export_connectivity(config, subject, session, regional_table, parcels.atlas_name, parcels.scale)
