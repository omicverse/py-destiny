"""pydestiny — pure-Python port of destiny (Angerer et al. Bioinformatics 2016).

Implements:
- DiffusionMap: distance → adaptive-kernel → density normalisation → eigendecomp
- DPT: diffusion pseudotime via accumulated transition probability
- Sigmas: local/global kernel bandwidth selection

v0.1 covers the core algorithmic API. Plotting + RNA-velocity-aware features
deferred to v0.2.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .diffmap import DiffusionMap, Sigmas, find_dm_k
from .dpt import DPT

__all__ = ["DiffusionMap", "DPT", "Sigmas", "find_dm_k", "__version__"]
