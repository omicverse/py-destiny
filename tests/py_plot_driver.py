"""Render py-destiny plots on Guo fixture."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))

import pydestiny
from pydestiny.plotting import plot_diffusion_map, plot_dpt
from ggplot2_py import ggsave


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(_PORT / "data/fixture_guo_expression.csv", index_col=0).to_numpy(dtype=np.float64)
    np.random.seed(42)
    dm = pydestiny.DiffusionMap.fit(expr, sigma="local", n_eigs=5, k=20)
    dpt = pydestiny.DPT(dm, root=0)

    p1 = plot_diffusion_map(dm, dims=(1, 2))
    ggsave(str(out_dir / "Py_diffmap.png"), plot=p1, width=6, height=4, dpi=100)

    p2 = plot_dpt(dm, dpt, dims=(1, 2))
    ggsave(str(out_dir / "Py_dpt.png"), plot=p2, width=6, height=4, dpi=100)
    print("done")


if __name__ == "__main__":
    main()
