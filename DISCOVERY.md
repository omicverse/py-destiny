# Discovery — py-destiny (scaffold; v0.0 not started)

## 1. Is this package already ported?

`python -m engine.discover_omicverse_deps --check destiny` → **No existing omicverse port found.**

## 2. Dependency audit + scope assessment

destiny (Angerer et al. *Bioinformatics* 2016, 341 citations) provides:
- Diffusion maps (DM)
- Diffusion pseudotime (DPT)
- Eigendecomposition with single-cell noise model (dropout-aware)
- Multi-scale eigenvectors (Coifman-style)

**Size**: 40 exports, ~3000 R LOC + significant Rcpp C++ extensions.

**Why this port is hard**:
- destiny's distance / kernel computation is in **C++ via Rcpp** (`RcppExports.R` is auto-generated; algorithmic code is in `src/*.cpp`). Porting needs to either:
  - (a) Re-implement the C++ kernel in Python/Cython/Numba
  - (b) Use `scanpy.tl.dpt` (already a Python implementation of DPT) and accept that it diverges from destiny in the noise-model details
- destiny's `Sigma` (local kernel bandwidth) selection is a custom algorithm not in scanpy.

**Reusable from omicverse**:
- `omicverse.single._diffusionmap` already provides a homebrew diffusion map. NOT exact destiny parity but functionally similar.
- scanpy.tl.dpt provides Diffusion Pseudotime — close to destiny's DPT.

## 3. Decision

**Skip a from-scratch port in this session**. Instead recommend:
- Make `omicverse.single._diffusionmap` more closely match destiny's API (cosmetic shim).
- Document `scanpy.tl.dpt` as the recommended DPT replacement.
- Defer C++ kernel port to a future major undertaking (~4 weeks).

## 4. v0.1 roadmap (for a future session)

1. Port `Sigma` (local bandwidth selection) from C++ to Python.
2. Port `DiffusionMap` constructor without dropout-noise model first; then add dropout layer.
3. Port `DPT` matching scanpy semantics.
4. Three notebooks comparing to scanpy.
