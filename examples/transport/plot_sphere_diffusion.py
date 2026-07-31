#!/usr/bin/env python3
"""Plot sphere diffusion results."""

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

for csv in sorted(DIR.glob("*_transport_sphere_diffusion.csv")):
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
    ax.set_title(r"$\bar{C}_{\mathrm{avg}} = \frac{C_0}{s}\,\frac{3}{q}\left[\coth q - \frac{1}{q}\right]$")
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
