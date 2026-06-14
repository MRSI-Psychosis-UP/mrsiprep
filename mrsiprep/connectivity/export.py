"""Connectivity export helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mrsiprep.connectivity.edges import build_edges
from mrsiprep.connectivity.matrix import build_regional_matrix
from mrsiprep.connectivity.nodes import build_nodes
from mrsiprep.connectivity.similarity import compute_similarity
from mrsiprep.io.naming import connectome_derivative


def export_connectivity(config, subject: str, session: str | None, regional_table: Path, atlas_name: str, scale: str | None = None) -> dict[str, Path]:
    matrix = build_regional_matrix(regional_table, value_col=config.regional_summary)
    sim = compute_similarity(matrix, method=config.connectivity_method)
    matrix_tsv = connectome_derivative(config.derivative_dir, subject, session, "tsv", atlas=atlas_name, scale=scale, desc="metsim_matrix")
    matrix_npy = connectome_derivative(config.derivative_dir, subject, session, "npy", atlas=atlas_name, scale=scale, desc="metsim_matrix")
    nodes_tsv = connectome_derivative(config.derivative_dir, subject, session, "tsv", atlas=atlas_name, scale=scale, desc="nodes")
    edges_tsv = connectome_derivative(config.derivative_dir, subject, session, "tsv", atlas=atlas_name, scale=scale, desc="edges")
    matrix_tsv.parent.mkdir(parents=True, exist_ok=True)
    sim.to_csv(matrix_tsv, sep="\t")
    np.save(matrix_npy, sim.to_numpy())
    build_nodes(regional_table).to_csv(nodes_tsv, sep="\t", index=False)
    build_edges(sim, config.connectivity_method).to_csv(edges_tsv, sep="\t", index=False)
    return {"matrix_tsv": matrix_tsv, "matrix_npy": matrix_npy, "nodes": nodes_tsv, "edges": edges_tsv}
