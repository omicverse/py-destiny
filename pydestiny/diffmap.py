"""DiffusionMap — destiny's core algorithm in pure NumPy/SciPy.

Algorithm (Coifman & Lafon 2006; Haghverdi et al. 2016 single-cell adaptation):

1. Compute pairwise distances D
2. Local sigma per cell σ_i = mean distance to k_local nearest neighbours
3. Build symmetric kernel: K_ij = exp(-D_ij² / (σ_i · σ_j))
4. Density-normalise: K̃_ij = K_ij / (deg_i^α · deg_j^α), α=1 (Markov chain on
   sample geometry)
5. Row-normalise: P_ij = K̃_ij / Σ_j K̃_ij  → diffusion transition matrix
6. Compute symmetric P_sym = D_d^{-½} P D_d^{½} (D_d = row-sums diag)
7. Eigendecomp top n_eigs+1 eigenvectors (the constant is excluded)

R destiny implements this with knn-graph sparse matrices + ARPACK. We use
scipy.sparse + eigsh for the same effect.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.sparse import csr_matrix, diags, eye, lil_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class Sigmas:
    """Per-cell local kernel bandwidths."""
    sigmas: np.ndarray
    n_local: int

    @classmethod
    def from_distance(cls, D: np.ndarray, n_local: int = 5) -> "Sigmas":
        """σ_i = mean distance to n_local nearest neighbours."""
        n = D.shape[0]
        sorted_d = np.sort(D, axis=1)
        # skip the diagonal zero
        sigma = sorted_d[:, 1:n_local + 1].mean(axis=1)
        sigma = np.maximum(sigma, 1e-10)
        return cls(sigmas=sigma, n_local=n_local)


def find_dm_k(n: int) -> int:
    """Default k for diffusion map: ceil(min(n - 1, 7 * (log10(n) + 1)))."""
    return int(np.ceil(min(n - 1, 7 * (np.log10(n) + 1))))


@dataclass
class DiffusionMap:
    """Diffusion map result.

    Attributes:
        eigenvalues: top eigenvalues (n_eigs,)
        eigenvectors: top eigenvectors (n_cells × n_eigs)
        transitions: row-normalised transition matrix (sparse)
        sigmas: per-cell local bandwidth used
        d_norms: per-cell density factor d_i = Σ_j K_ij
    """
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    transitions: csr_matrix
    sigmas: Sigmas
    d_norms: np.ndarray
    distance_method: str

    @classmethod
    def fit(
        cls,
        data: np.ndarray,
        *,
        sigma: str | float = "local",
        k: int | None = None,
        n_eigs: int = 20,
        n_local: int = 5,
        density_norm: bool = True,
        distance: str = "euclidean",
    ) -> "DiffusionMap":
        """Fit a diffusion map to `data` (cells × features).

        Args:
            data: (n_cells × n_features) numeric matrix
            sigma: "local" (per-cell adaptive) or "global" (single value) or float
            k: number of nearest neighbours (default: `find_dm_k(n)`)
            n_eigs: number of eigenvectors to return
            n_local: for local sigma, number of nearest neighbours to average
            density_norm: apply Markov-chain density normalisation
            distance: "euclidean" (default), "cosine"
        """
        X = np.asarray(data, dtype=np.float64)
        n = X.shape[0]
        if k is None:
            k = find_dm_k(n)
        # Full pairwise distance matrix
        from scipy.spatial.distance import cdist
        if distance == "euclidean":
            D = cdist(X, X, "euclidean")
        elif distance == "cosine":
            D = cdist(X, X, "cosine")
        elif distance == "l2":
            D = cdist(X, X, "euclidean")
        else:
            raise ValueError(f"unknown distance: {distance}")

        # Sigma
        if sigma == "local":
            sigmas = Sigmas.from_distance(D, n_local=n_local)
        else:
            sig_val = np.median(np.sort(D, axis=1)[:, 1:n_local + 1].mean(axis=1)) \
                      if sigma == "global" else float(sigma)
            sigmas = Sigmas(sigmas=np.full(n, sig_val), n_local=n_local)

        # Build kernel with adaptive bandwidth
        sigma_outer = np.outer(sigmas.sigmas, sigmas.sigmas)
        D2 = D ** 2
        with np.errstate(divide="ignore", over="ignore"):
            K = np.exp(-D2 / sigma_outer)
        np.fill_diagonal(K, 0.0)

        # Truncate to top-k kNN per row for sparsity
        kth_largest = np.partition(K, n - k - 1, axis=1)[:, n - k - 1]
        mask = K >= kth_largest[:, None]
        K = K * mask
        # Symmetrise
        K = np.maximum(K, K.T)
        np.fill_diagonal(K, 0.0)

        # Density normalisation
        d = K.sum(axis=1)
        d = np.maximum(d, 1e-12)
        if density_norm:
            alpha = 1.0
            K = K / np.outer(d ** alpha, d ** alpha)
            d = K.sum(axis=1)
            d = np.maximum(d, 1e-12)

        # Transition matrix P = D^{-1} K (row-stochastic)
        D_inv = 1.0 / d
        P = K * D_inv[:, None]   # row-normalise
        P_sym = (np.sqrt(d)[:, None] * P) / np.sqrt(d)[None, :]
        # eigsh wants symmetric → symmetrise tiny rounding
        P_sym = (P_sym + P_sym.T) / 2

        n_eigs = min(n_eigs, n - 2)
        # Largest eigenvalues
        eigvals, eigvecs = eigsh(P_sym, k=n_eigs + 1, which="LA")
        # Sort descending; the largest is the trivial constant (eigenvalue 1)
        order = np.argsort(eigvals)[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]
        # Discard the trivial first eigenvector
        eigvals = eigvals[1:]
        eigvecs = eigvecs[:, 1:]
        # Transform back to non-symmetric basis: φ_i = ψ_i / sqrt(d_i)
        eigvecs = eigvecs / np.sqrt(d)[:, None]

        return cls(
            eigenvalues=eigvals,
            eigenvectors=eigvecs,
            transitions=csr_matrix(P),
            sigmas=sigmas,
            d_norms=d,
            distance_method=distance,
        )
