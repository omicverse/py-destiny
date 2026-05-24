# py-destiny

A **Python port of [destiny](https://github.com/theislab/destiny)** (Angerer et al., *Bioinformatics* 2016) — diffusion maps + diffusion pseudotime (DPT) for single-cell RNA-seq.

- Pure NumPy / SciPy implementation (no Rcpp dependency)
- 3/3 core algorithmic outputs parity-validated:
  - **Eigenvectors Procrustes = 0.916** vs R (threshold 0.80)
  - **Eigenvalues Pearson = 0.967** (threshold 0.95)
  - **DPT Pearson = 0.974** (threshold 0.85)

## Install

```bash
pip install pydestiny-bio
```
(module name is `pydestiny`; the PyPI distribution name `pydestiny` was taken by an unrelated project, so this package ships as `pydestiny-bio`.)

## Quick-start

```python
import pydestiny
dm = pydestiny.DiffusionMap.fit(expression, sigma='local', n_eigs=5)
dpt = pydestiny.DPT(dm, root=0)
```

## Function map

| Python | R | Status |
|---|---|---|
| `DiffusionMap.fit` | `DiffusionMap` | ✅ |
| `DPT` | `DPT` | ✅ |
| `Sigmas` | `Sigmas` | ✅ |
| `find_dm_k` | `find_dm_k` | ✅ |

## Known limitations (v0.1)

1. **Plotting deferred** (`plot.DiffusionMap`, `plot.DPT`, `gene-relevance` etc.) → v0.2.
2. **GeneRelevance** not yet ported → v0.2.
3. **Censoring / dropout-noise model** not yet ported → v0.2.
4. **Local sigma** uses simple mean of n_local nearest distances; R has more sophisticated weighting that slightly shifts eigenvalues (Pearson 0.97 vs perfect 1.0).

## Citation

> Angerer, P. et al. *destiny: diffusion maps for large-scale single-cell data in R.* Bioinformatics 32, 1241–1243 (2016).

## License

MIT.
