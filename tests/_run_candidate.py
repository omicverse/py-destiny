"""Candidate runner — runs pydestiny on the same Guo data."""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))
import pydestiny


def main():
    fixture_path, output_path = sys.argv[1], sys.argv[2]
    pp = Path(fixture_path)
    expr = pd.read_csv(pp.with_name(pp.stem + "_expression.csv"), index_col=0).to_numpy(dtype=np.float64)
    print(f"[cand] expression shape: {expr.shape}")

    np.random.seed(42)
    dm = pydestiny.DiffusionMap.fit(expr, sigma="local", n_eigs=5, k=None, n_local=5)
    print(f"[cand] DM: {dm.eigenvalues.shape}, eigvecs {dm.eigenvectors.shape}")

    dpt_vec = pydestiny.DPT(dm, root=0)
    print(f"[cand] DPT length: {len(dpt_vec)}")

    out = {
        "eigenvalues": dm.eigenvalues.tolist(),
        "eigenvectors": dm.eigenvectors.tolist(),
        "dpt": dpt_vec.tolist(),
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote {output_path}")


if __name__ == "__main__":
    main()
