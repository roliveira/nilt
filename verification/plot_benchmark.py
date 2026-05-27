#!/usr/bin/env python3
"""Plot benchmark timing results."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser

TIMING_STYLES = {
    ("Stehfest", "C++"):    {"color": "C0", "marker": "s", "ls": "-"},
    ("Stehfest", "Python"): {"color": "C0", "marker": "s", "ls": "--"},
    ("Talbot",   "C++"):    {"color": "C1", "marker": "^", "ls": "-"},
    ("Talbot",   "Python"): {"color": "C1", "marker": "^", "ls": "--"},
    ("DeHoog",   "C++"):    {"color": "C2", "marker": "o", "ls": "-"},
    ("DeHoog",   "Python"): {"color": "C2", "marker": "o", "ls": "--"},
}


def main() -> None:
    args = create_argument_parser(__doc__, __file__).parse_args()

    cpp_csv = args.dir / "benchmark_timing.csv"
    py_csv  = args.dir / "benchmark_timing_python.csv"

    if not cpp_csv.exists() and not py_csv.exists():
        print("  No benchmark timing CSVs found, skipping")
        return

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    for csv_path, lang in [(cpp_csv, "C++"), (py_csv, "Python")]:
        if not csv_path.exists():
            continue

        df = pl.read_csv(csv_path)

        for method in ["Stehfest", "Talbot", "DeHoog"]:
            sub = df.filter(pl.col("method") == method)
            if sub.is_empty():
                continue

            style = TIMING_STYLES[(method, lang)]
            ax.loglog(
                sub["param"], sub["time_us"],
                marker=style["marker"], color=style["color"],
                ls=style["ls"], label=f"{method} ({lang})",
            )

    ax.set_xlabel("Algorithm parameter (N / n / M)")
    ax.set_ylabel(r"Time per inversion [$\mu$s]")
    ax.set_title(r"Timing: $F(s) = 1/(s+1)$ at $t = 1$")
    ax.legend(ncol=2)
    ax.grid(True, alpha=0.3, which="both")

    out = args.dir / f"benchmark_timing.{args.format}"
    fig.savefig(out, dpi=600)
    plt.close(fig)
    print(f"  {out}")


if __name__ == "__main__":
    main()
