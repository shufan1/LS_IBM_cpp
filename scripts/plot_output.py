#!/usr/bin/env python3
"""Visualize LS_IBM_cpp output.

The C++ solver writes:
  - coordinates.json : {imax, jmax, lx, ly, xu, yu, xv, yv, xp, yp}
  - dataRDE<N>dt.json: {time, U, V, P, phi: {phi_0, phi_1, ...}, psi}

Every field is a FLAT, row-major array (index = i*ny + j), and each lives on
its own staggered grid:
  U        -> (imax,   jmax+1) on (xu, yu)
  V        -> (imax+1, jmax  ) on (xv, yv)
  P/phi/psi-> (imax+1, jmax+1) on (xp, yp)

Usage:
  python plot_output.py <coordinates.json> <dataRDE_N_dt.json> [out.png]

On Sherlock (not the login node -- use sh_dev or a job):
  ml python py-numpy py-matplotlib
  python plot_output.py coordinates.json Output/dataRDE1dt.json
"""
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless: write a PNG, no display needed on the cluster
import matplotlib.pyplot as plt


def load(path):
    with open(path) as f:
        return json.load(f)


def grid2d(flat, nx, ny):
    """Flat row-major (index=i*ny+j) -> 2D array shaped (nx, ny), i.e. [x, y]."""
    a = np.asarray(flat, dtype=float)
    assert a.size == nx * ny, f"expected {nx*ny} values, got {a.size}"
    return a.reshape(nx, ny)


def field_panel(ax, x, y, field_xy, title, cmap="viridis"):
    # field_xy is [x, y]; pcolormesh wants C shaped (len(y), len(x)) -> transpose.
    pc = ax.pcolormesh(x, y, field_xy.T, shading="auto", cmap=cmap)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.colorbar(pc, ax=ax, shrink=0.8)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    coords = load(sys.argv[1])
    data = load(sys.argv[2])
    out_png = sys.argv[3] if len(sys.argv) > 3 else "output.png"

    imax, jmax = coords["imax"], coords["jmax"]
    xu, yu = coords["xu"], coords["yu"]      # len imax,   jmax+1
    xv, yv = coords["xv"], coords["yv"]      # len imax+1, jmax
    xp, yp = coords["xp"], coords["yp"]      # len imax+1, jmax+1

    U = grid2d(data["U"], imax,     jmax + 1)
    V = grid2d(data["V"], imax + 1, jmax)
    P = grid2d(data["P"], imax + 1, jmax + 1)
    psi = grid2d(data["psi"], imax + 1, jmax + 1)
    phi = {k: grid2d(v, imax + 1, jmax + 1) for k, v in data["phi"].items()}

    nphi = len(phi)
    panels = 4 + nphi  # U, V, P, psi, + one per species
    ncol = 3
    nrow = (panels + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 3.2 * nrow))
    axes = np.atleast_1d(axes).ravel()

    field_panel(axes[0], xu, yu, U, "U (x-velocity)", cmap="RdBu_r")
    field_panel(axes[1], xv, yv, V, "V (y-velocity)", cmap="RdBu_r")
    field_panel(axes[2], xp, yp, P, "P (pressure)", cmap="viridis")
    field_panel(axes[3], xp, yp, psi, "psi (level set)", cmap="coolwarm")
    # grain surface = psi==0 contour, overlaid on the level-set panel
    axes[3].contour(xp, yp, psi.T, levels=[0.0], colors="k", linewidths=1.2)

    for k, (name, ph) in enumerate(sorted(phi.items())):
        ax = axes[4 + k]
        field_panel(ax, xp, yp, ph, name, cmap="magma")
        ax.contour(xp, yp, psi.T, levels=[0.0], colors="w", linewidths=0.8)

    for ax in axes[panels:]:
        ax.axis("off")

    fig.suptitle(f"t = {data.get('time', 0.0):.6g}   (grid {imax}x{jmax})")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
