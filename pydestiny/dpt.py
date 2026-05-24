"""DPT — Diffusion Pseudotime.

Haghverdi et al. *Nature Methods* 2016. Given a fitted DiffusionMap and one
or more root cells, compute pseudotime as the cumulative-transition-distance
from the root.

Math:
  DPT(i, j) = sqrt( Σ_k (eigenvector_k[i] - eigenvector_k[j])² / (1 - λ_k)² )

Returns a scalar per cell when given one root.
"""

from __future__ import annotations

import numpy as np

from .diffmap import DiffusionMap


def DPT(
    dm: DiffusionMap,
    *,
    root: int = 0,
    n_eigs_dpt: int | None = None,
) -> np.ndarray:
    """Diffusion pseudotime from a root cell.

    Args:
        dm: fitted DiffusionMap
        root: index of the root cell
        n_eigs_dpt: number of eigenvectors to use in the DPT sum (default: all
                    non-trivial eigenvectors of dm)

    Returns:
        (n_cells,) per-cell pseudotime, normalised to [0, 1].
    """
    eigvals = dm.eigenvalues
    eigvecs = dm.eigenvectors
    if n_eigs_dpt is None:
        n_eigs_dpt = eigvecs.shape[1]
    else:
        n_eigs_dpt = min(n_eigs_dpt, eigvecs.shape[1])

    # Use eigenvalues' 1/(1-λ) weighting (the geometric-series sum of transition powers)
    w = 1.0 / np.maximum(1.0 - eigvals[:n_eigs_dpt], 1e-9)
    # Distance from root in weighted-eigenvector space
    diff = eigvecs[:, :n_eigs_dpt] - eigvecs[root, :n_eigs_dpt]
    dist = np.sqrt(np.sum((diff ** 2) * w ** 2, axis=1))
    # Normalise to [0, 1]
    if dist.max() > dist.min():
        return (dist - dist.min()) / (dist.max() - dist.min())
    return np.zeros_like(dist)
