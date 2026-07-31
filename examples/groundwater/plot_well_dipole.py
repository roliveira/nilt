#!/usr/bin/env python3
"""Plot well dipole results (time series and 2D spatial)."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DIR = Path(".")
FMT = "png"
DPI = 600
MARKEVERY = 5

METHODS = [
    ("stehfest", "s", "C0"),
    ("talbot",   "^", "C1"),
    ("dehoog",   "o", "C2"),
]

# --- Time series ---
for csv in sorted(DIR.glob("*_groundwater_well_dipole.csv")):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    t   = data["t"]
    ana = data["analytical"]
    methods = [(c, m, clr) for c, m, clr in METHODS if c in data.dtype.names]

    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    ax.plot(t, ana, "k-", label="analytical", zorder=10)
    for col, marker, color in methods:
        ax.plot(t, data[col], marker, color=color, markevery=MARKEVERY, label=col.title())

    ax.set_xscale("log")
    ax.set_ylabel("Net drawdown [m]")
    ax.set_title(r"Well dipole — observation at $(5, 0)$ m")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in methods:
        ax_err.semilogy(t, np.abs(data[col] - ana) / ana_abs, marker, color=color, markevery=MARKEVERY)

    ax_err.set_xlabel("Time [s]")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)

# --- 2D spatial ---
for csv in sorted(DIR.glob("*_groundwater_well_dipole_spatial.csv")):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    xs = np.unique(data["x"])
    ys = np.unique(data["y"])
    nx, ny = len(xs), len(ys)
    Z = data["talbot"].reshape(ny, nx)
    X, Y = np.meshgrid(xs, ys)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    cf = ax.contourf(X, Y, Z, levels=30)
    fig.colorbar(cf, ax=ax, label="Net drawdown [m]")
    ax.plot(-15, -5, "kv", ms=8, mew=1.5, label="Pumping (0.010 m³/s)")
    ax.plot( 20,  8, "k^", ms=8, mew=1.5, label="Injection (0.007 m³/s)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(r"Well dipole at $t = 2$ h")
    ax.set_aspect("equal")
    ax.legend(loc="lower left")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)
