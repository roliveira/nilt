#!/usr/bin/env python3
"""Plot benchmark timing results (C++ and Python side by side).

Reads:
  {cpp,py}_benchmark_timing.csv   (algo-parameter sweep, scalar t)
  {cpp,py}_benchmark_array.csv    (array-size sweep, default params)

Emits:
  benchmark_timing.png
  benchmark_array.png             (only if the array CSVs exist)
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DIR = Path(".")
FMT = "png"
DPI = 600

STYLES = {
    ("Stehfest", "cpp"): {"color": "C0", "marker": "s", "ls": "-"},
    ("Stehfest", "py"):  {"color": "C0", "marker": "s", "ls": "--"},
    ("Talbot",   "cpp"): {"color": "C1", "marker": "^", "ls": "-"},
    ("Talbot",   "py"):  {"color": "C1", "marker": "^", "ls": "--"},
    ("DeHoog",   "cpp"): {"color": "C2", "marker": "o", "ls": "-"},
    ("DeHoog",   "py"):  {"color": "C2", "marker": "o", "ls": "--"},
}


def plot_sweep(csv_glob, xcol, xlabel, title, outfile):
    csvs = sorted(DIR.glob(csv_glob))
    if not csvs:
        return
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for csv in csvs:
        prefix = csv.stem.split("_")[0]  # "cpp" or "py"
        data = np.genfromtxt(csv, delimiter=",", dtype=None, names=True, encoding=None)
        for method in ["Stehfest", "Talbot", "DeHoog"]:
            mask = data["method"] == method
            if not mask.any():
                continue
            style = STYLES[(method, prefix)]
            ax.loglog(
                data[xcol][mask], data["time_us"][mask],
                marker=style["marker"], color=style["color"],
                ls=style["ls"], label=f"{method} ({prefix})",
            )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Time per inversion [$\mu$s]")
    ax.set_title(title)
    ax.legend(ncol=2)
    ax.grid(True, alpha=0.3, which="both")
    fig.savefig(outfile, dpi=DPI)
    plt.close(fig)
    print(outfile)


plot_sweep(
    "*_benchmark_timing.csv",
    xcol="param",
    xlabel="Algorithm parameter (N or M)",
    title=r"Scalar timing: $F(s) = 1/(s+1)$ at $t = 1$",
    outfile=f"benchmark_timing.{FMT}",
)

plot_sweep(
    "*_benchmark_array.csv",
    xcol="nt",
    xlabel="Array length (number of t-values)",
    title=r"Array timing at default params: $F(s) = 1/(s+1)$",
    outfile=f"benchmark_array.{FMT}",
)
