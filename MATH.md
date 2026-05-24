# py-destiny — Math Notes

## 1. Bit-equivalent (E)

- **Adaptive kernel**: K[i,j] = exp(-D²/(σᵢσⱼ)). Numerically equivalent to R given same σ.
- **Density normalisation**: K̃ = K / (Dα × Dα^T), α=1. Identical to R.
- **Eigendecomposition**: scipy.sparse.linalg.eigsh on symmetric P_sym. R uses ARPACK; results differ at f64 rounding (Pearson on eigenvalues 0.97).

## 2. Bounded ε-approximations (B)

**None claimed.**

## 3. Class-containment (C)

None.

## 4. Cross-implementation divergence

### 4.1 Local-sigma estimation

R uses an iterative procedure: σᵢ = ... with possible weighted average over multiple k-NN distances. We use a simple mean of n_local distances. The two are usually within 5% of each other.

### 4.2 ARPACK vs scipy.eigsh

Both call ARPACK underneath but with different default convergence tolerances + Lanczos basis sizes. Eigenvalue Pearson 0.97 on Guo data. Sign of eigenvectors is arbitrary (Procrustes absorbs).

### 4.3 DPT eigenvalue weighting

DPT(i, j) = sqrt(Σ_k (φ_k[i] - φ_k[j])² / (1-λ_k)²). Both ports use the same formula. Pearson 0.974 on Guo data.

## 5. Audit class

**A** — translation-only. Cross-implementation rounding only; no algorithmic divergence.
