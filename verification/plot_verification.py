#!/usr/bin/env python3
"""Plot verification results for all test functions."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DIR = Path(".")
FMT = "png"
DPI = 600

METHODS = [
    ("Stehfest", "s", "C0"),
    ("Talbot",   "^", "C1"),
    ("DeHoog",   "o", "C2"),
]

FUNC_LABELS = {
    1:  r"$f(t) = 1/\sqrt{\pi t}$",
    2:  r"$f(t) = -\gamma - \ln t$",
    3:  r"$f(t) = t^3/6$",
    4:  r"$f(t) = e^{-t}$",
    5:  r"$f(t) = \sin\sqrt{2t}$",
    6:  r"$f(t) = t$",
    7:  r"$f(t) = t\,e^{-t}$",
    8:  r"$f(t) = \sin t$",
    9:  r"$f(t) = \cos t$",
    10: r"$f(t) = e^{-t}\sin t$",
}

for func in range(1, 11):
    fig, (ax_val, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    plotted_analytical = False
    for method, marker, color in METHODS:
        csvs = sorted(DIR.glob(f"*_func{func}_{method}.csv"))
        for csv in csvs:
            data = np.genfromtxt(csv, delimiter=",", names=True)

            if not plotted_analytical:
                ax_val.plot(data["t"], data["fta"], "k-", label="analytical", zorder=10)
                plotted_analytical = True

            prefix = csv.stem.split("_")[0]
            ax_val.plot(data["t"], data["ftn"], marker, color=color, label=f"{method} ({prefix})")
            ax_err.semilogy(data["t"], data["err"], marker, color=color)

    ax_val.set_title(FUNC_LABELS[func])
    ax_val.set_ylabel(r"$f(t)$")
    ax_val.legend(loc="best", fontsize=8)
    ax_val.grid(True, alpha=0.3)

    ax_err.set_xlabel(r"$t$")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = f"verification_func{func}.{FMT}"
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(out)
