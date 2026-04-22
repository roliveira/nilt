#!/usr/bin/env python3
"""Plot sphere diffusion results."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser


def main() -> None:
    args = create_argument_parser(__doc__, __file__).parse_args()

    csv_path = args.dir / "transport_sphere_diffusion.csv"
    if not csv_path.exists():
        print(f"  {csv_path} not found, skipping")
        return

    df = pl.read_csv(csv_path)
    t = df["t"].to_numpy()
    ana = df["analytical"].to_numpy()

    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
        sharex=True, constrained_layout=True,
    )

    ax.plot(t, ana, "k-", label="analytical", zorder=10)
    ax.plot(t, df["stehfest"], "s", color="C0", markevery=5, label="Stehfest")
    ax.plot(t, df["talbot"], "^", color="C1", markevery=5, label="Talbot")
    ax.plot(t, df["dehoog"], "o", color="C2", markevery=5, label="De Hoog")
    ax.set_ylabel("Average concentration [mol/m³]")
    ax.set_title(r"$\bar{C}_{\mathrm{avg}} = \frac{C_0}{s}\,\frac{3}{q}\left[\coth q - \frac{1}{q}\right]$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ana_abs = np.maximum(np.abs(ana), 1e-30)
    for col, marker, color in [("stehfest", "s", "C0"), ("talbot", "^", "C1"), ("dehoog", "o", "C2")]:
        rel_err = np.abs(df[col].to_numpy() - ana) / ana_abs
        ax_err.semilogy(t, rel_err, marker=marker, color=color, markevery=5)

    ax_err.set_xlabel("Time [s]")
    ax_err.set_ylabel("Relative error")
    ax_err.grid(True, alpha=0.3, which="both")

    out = csv_path.with_suffix(f".{args.format}")
    fig.savefig(out, dpi=600)
    plt.close(fig)
    print(f"  {out}")


if __name__ == "__main__":
    main()
