#!/usr/bin/env python3
"""Plot Theis well drawdown results (time series and distance profile)."""

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
for csv in sorted(DIR.glob("*_groundwater_theis_well_time.csv")):
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

    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Drawdown [ft]")
    ax.set_title(r"Theis well — $r = 100$ ft")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in METHODS:
        ax_err.semilogy(t, np.abs(data[col] - ana) / ana_abs, marker, color=color, markevery=MARKEVERY)

    ax_err.set_xlabel("Time [days]")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)

# --- Distance profile ---
for csv in sorted(DIR.glob("*_groundwater_theis_well_distance.csv")):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    r   = data["r"]
    ana = data["analytical"]

    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    ax.plot(r, ana, "k-", label="analytical", zorder=10)
    for col, marker, color in METHODS:
        ax.plot(r, data[col], marker, color=color, markevery=MARKEVERY, label=col.title())

    ax.invert_yaxis()
    ax.set_ylabel("Drawdown [ft]")
    ax.set_title(r"Theis well — $t = 2$ days")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in METHODS:
        ax_err.semilogy(r, np.abs(data[col] - ana) / ana_abs, marker, color=color, markevery=MARKEVERY)

    ax_err.set_xlabel("Distance from well [ft]")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv.with_suffix(f".{FMT}")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)
