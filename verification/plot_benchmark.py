#!/usr/bin/env python3
"""Plot benchmark timing results (C++ and Python side by side)."""

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

fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

for csv in sorted(DIR.glob("*_benchmark_timing.csv")):
    prefix = csv.stem.split("_")[0]  # "cpp" or "py"
    data = np.genfromtxt(csv, delimiter=",", dtype=None, names=True, encoding=None)

    for method in ["Stehfest", "Talbot", "DeHoog"]:
        mask = data["method"] == method
        if not mask.any():
            continue
        style = STYLES[(method, prefix)]
        ax.loglog(
            data["param"][mask], data["time_us"][mask],
            marker=style["marker"], color=style["color"],
            ls=style["ls"], label=f"{method} ({prefix})",
        )

ax.set_xlabel("Algorithm parameter (N or M)")
ax.set_ylabel(r"Time per inversion [$\mu$s]")
ax.set_title(r"Timing: $F(s) = 1/(s+1)$ at $t = 1$")
ax.legend(ncol=2)
ax.grid(True, alpha=0.3, which="both")

out = f"benchmark_timing.{FMT}"
fig.savefig(out, dpi=DPI)
plt.close(fig)
print(out)
