#!/usr/bin/env python3
"""Plot Theis well drawdown results (time series and distance profile)."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser


def plot_csv(csv_path, xcol, xlabel, ylabel, title, fmt, xscale=None, invert_y=False):
    df = pl.read_csv(csv_path)
    x = df[xcol].to_numpy()
    ana = df["analytical"].to_numpy()

    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    ax.plot(x, ana, "k-", label="analytical", zorder=10)
    ax.plot(x, df["stehfest"], "s", color="C0", markevery=5, label="Stehfest")
    ax.plot(x, df["talbot"], "^", color="C1", markevery=5, label="Talbot")
    ax.plot(x, df["dehoog"], "o", color="C2", markevery=5, label="De Hoog")

    if xscale:
        ax.set_xscale(xscale)
    if invert_y:
        ax.invert_yaxis()
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in [("stehfest", "s", "C0"), ("talbot", "^", "C1"), ("dehoog", "o", "C2")]:
        rel_err = np.abs(df[col].to_numpy() - ana) / ana_abs
        ax_err.semilogy(x, rel_err, marker=marker, color=color, markevery=5)

    ax_err.set_xlabel(xlabel)
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv_path.with_suffix(f".{fmt}")
    fig.savefig(out, dpi=600)
    plt.close(fig)
    print(f"  {out}")


def main() -> None:
    args = create_argument_parser(__doc__, __file__).parse_args()

    # Time series
    csv_time = args.dir / "groundwater_theis_well_time.csv"
    if csv_time.exists():
        plot_csv(
            csv_time, "t", "Time [days]", "Drawdown [ft]",
            r"Theis well — $s(r,t) = \frac{114.6\,Q}{T}\,W(u)$" "\n" r"$r = 100$ ft",
            args.format, xscale="log", invert_y=True,
        )

    # Distance profile
    csv_dist = args.dir / "groundwater_theis_well_distance.csv"
    if csv_dist.exists():
        plot_csv(
            csv_dist, "r", "Distance from pumped well [ft]", "Drawdown [ft]",
            r"Theis well — $s(r,t) = \frac{114.6\,Q}{T}\,W(u)$" "\n" r"$t = 2$ days (48 h)",
            args.format, invert_y=True,
        )


if __name__ == "__main__":
    main()
