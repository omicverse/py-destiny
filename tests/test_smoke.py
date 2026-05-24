import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pydestiny


def test_import():
    assert pydestiny.__version__.startswith("0.1")
    for fn in ("DiffusionMap", "DPT", "Sigmas", "find_dm_k"):
        assert hasattr(pydestiny, fn)


def test_find_dm_k():
    assert pydestiny.find_dm_k(100) > 0
    assert pydestiny.find_dm_k(10000) > pydestiny.find_dm_k(100)


def test_diffusion_map_runs():
    rng = np.random.default_rng(42)
    # Two clusters in 5D, 80 cells
    X = np.vstack([rng.normal(0, 1, (40, 5)), rng.normal(5, 1, (40, 5))])
    dm = pydestiny.DiffusionMap.fit(X, sigma="local", n_eigs=3, n_local=5)
    assert dm.eigenvalues.shape == (3,)
    assert dm.eigenvectors.shape == (80, 3)
    assert dm.eigenvalues[0] <= 1.0 and dm.eigenvalues[0] > 0


def test_dpt_runs():
    rng = np.random.default_rng(42)
    # A 1D linear trajectory in 5D
    X = np.zeros((50, 5))
    X[:, 0] = np.linspace(0, 10, 50)
    X += rng.normal(0, 0.5, X.shape)
    dm = pydestiny.DiffusionMap.fit(X, sigma="local", n_eigs=3, n_local=5)
    dpt = pydestiny.DPT(dm, root=0)
    assert dpt.shape == (50,)
    assert 0 <= dpt.min() and dpt.max() <= 1
    # DPT should monotonically correlate with the underlying x-coordinate
    from scipy.stats import spearmanr
    rho = abs(spearmanr(dpt, X[:, 0])[0])
    assert rho > 0.5, f"DPT-vs-trajectory rho {rho}"
