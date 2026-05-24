"""Visualisation — ggplot2-python port of destiny's plot.DiffusionMap + plot.DPT.

R's ``plot.DiffusionMap`` is a 2D/3D scatter in eigenvector space, optionally
coloured by a vector. ``plot.DPT`` overlays DPT pseudotime on the same
embedding. We expose both as simple functions returning a ggplot2-python plot.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    geom_point,
    ggplot,
    labs,
    scale_color_gradientn,
    scale_color_manual,
    theme_classic,
)

from .diffmap import DiffusionMap


def _viridis_palette() -> list[str]:
    return [
        "#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
        "#1F9E89", "#35B779", "#6CCE59", "#B4DE2C", "#FDE725",
    ]


def plot_diffusion_map(
    dm: DiffusionMap,
    *,
    dims: tuple[int, int] = (1, 2),
    col=None,
    pal: list[str] | None = None,
    point_size: float = 1.5,
    point_alpha: float = 1.0,
    title: str | None = None,
):
    """1:1 port of ``destiny::plot.DiffusionMap`` (2D variant).

    Args:
        dm: fitted ``DiffusionMap``.
        dims: 1-indexed DC components to plot (default DC1 vs DC2).
        col: optional per-cell vector — numeric → gradient, categorical → discrete.
        pal: palette (list of colour strings).
        point_size, point_alpha: aesthetics.
    """
    i, j = int(dims[0]) - 1, int(dims[1]) - 1
    eigvec = dm.eigenvectors
    df = pd.DataFrame(
        {
            f"DC{i + 1}": eigvec[:, i],
            f"DC{j + 1}": eigvec[:, j],
        }
    )
    x_col, y_col = df.columns

    if col is None:
        p = (
            ggplot(df, aes(x=x_col, y=y_col))
            + geom_point(size=point_size, alpha=point_alpha)
        )
    else:
        col_arr = pd.Series(col)
        df["__col__"] = col_arr.values
        is_numeric = pd.api.types.is_numeric_dtype(col_arr)
        p = (
            ggplot(df, aes(x=x_col, y=y_col, colour="__col__"))
            + geom_point(size=point_size, alpha=point_alpha)
        )
        if is_numeric:
            p = p + scale_color_gradientn(colours=pal or _viridis_palette())
        else:
            unique = sorted(df["__col__"].unique(), key=str)
            colours = pal if pal else _viridis_palette()
            palette_dict = dict(zip(unique, (colours * (len(unique) // len(colours) + 1))[: len(unique)]))
            p = p + scale_color_manual(values=palette_dict)
    p = p + theme_classic() + labs(x=x_col, y=y_col, colour="" if col is None else "value")
    if title is not None:
        p = p + labs(title=title)
    return p


def plot_dpt(
    dm: DiffusionMap,
    dpt_vec: np.ndarray,
    *,
    dims: tuple[int, int] = (1, 2),
    point_size: float = 1.5,
    title: str | None = "DPT pseudotime",
):
    """1:1 port of ``destiny::plot.DPT``."""
    return plot_diffusion_map(
        dm,
        dims=dims,
        col=np.asarray(dpt_vec, dtype=np.float64),
        pal=[
            "#0000FF", "#3366FF", "#66CCFF", "#FFFFFF",
            "#FFCC66", "#FF6633", "#FF0000",
        ],
        point_size=point_size,
        title=title,
    )
