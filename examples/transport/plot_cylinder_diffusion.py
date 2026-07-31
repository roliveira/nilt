#!/usr/bin/env python3
"""Plot cylinder diffusion results (time series and radial cross-section)."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DIR = Path(__file__).parent / "build"
FMT = "png"
DPI = 600
MARKEVERY = 5

METHODS = [
    ("stehfest", "s", "C0"),
    ("talbot",   "^", "C1"),
    ("dehoog",   "o", "C2"),
]

# --- Time series ---
for csv in sorted(DIR.glob("*_transport_cylinder_diffusion.csv")):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    t   = data["t"]
    ana = data["analytical"]

    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    ax.plot(t, ana, "k-", label="analytical", zorder=10)
    for col, marker, color in METHODS:
        ax.plot(t, data[col], marker, color=color, markevery=MARKEVERY, label=col.title())

    ax.set_ylabel("Average concentration [mol/m³]")
    ax.set_title(r"$\bar{C}_{\mathrm{avg}} = \frac{C_0}{s}\,\frac{2}{q}\,\frac{I_1(q)}{I_0(q)}$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in METHODS:
        ax_err.semilogy(t, np.abs(data[col] - ana) / ana_abs, marker, color=color, markevery=MARKEVERY)

    ax_err.set_xlabel("Time [s]")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)

# --- Radial cross-section ---
for csv in sorted(DIR.glob("*_transport_cylinder_diffusion_spatial.csv")):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    r_data = data["r"]
    C_data = data["talbot"]

    a = r_data[-1]
    n = 301
    x = np.linspace(-a * 1.5, a * 1.5, n)
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X**2 + Y**2)
    Z = np.where(R <= a, np.interp(R.ravel(), r_data, C_data).reshape(n, n), 0.0)

    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    cf = ax.contourf(X * 1e3, Y * 1e3, Z, levels=30)
    fig.colorbar(cf, ax=ax, label="Concentration [mol/m³]")
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(a * 1e3 * np.cos(theta), a * 1e3 * np.sin(theta), "w--")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(r"$C(r, t=500\,\mathrm{s})$")
    ax.set_aspect("equal")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)
