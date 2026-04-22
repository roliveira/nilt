#!/usr/bin/env python3
"""Plot cylinder diffusion results (time series and 2D cross-section)."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser


def main() -> None:
    args = create_argument_parser(__doc__, __file__).parse_args()

    csv_path = args.dir / "transport_cylinder_diffusion.csv"
    if csv_path.exists():
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
        ax.set_title(r"$\bar{C}_{\mathrm{avg}} = \frac{C_0}{s}\,\frac{2}{q}\,\frac{I_1(q)}{I_0(q)}$")
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

    spatial = args.dir / "transport_cylinder_diffusion_spatial.csv"
    if spatial.exists():
        df = pl.read_csv(spatial)
        r_data = df["r"].to_numpy()
        C_data = df["talbot"].to_numpy()

        a = r_data[-1]
        pad = 1.5
        n = 301
        x = np.linspace(-a * pad, a * pad, n)
        y = np.linspace(-a * pad, a * pad, n)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)

        Z = np.where(R <= a, np.interp(R.ravel(), r_data, C_data).reshape(n, n), 0.0)

        fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
        cf = ax.contourf(X * 1e3, Y * 1e3, Z, levels=30)
        cb = fig.colorbar(cf, ax=ax)
        cb.set_label("Concentration [mol/m³]")
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(a * 1e3 * np.cos(theta), a * 1e3 * np.sin(theta), "w--")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(r"$C_{\mathrm{avg}} = \frac{2}{q}\,\frac{I_1(q)}{I_0(q)}$ at $t = 500$ s")
        ax.set_aspect("equal")

        out = spatial.with_name("transport_cylinder_diffusion_spatial").with_suffix(f".{args.format}")
        fig.savefig(out, dpi=600)
        plt.close(fig)
        print(f"  {out}")


if __name__ == "__main__":
    main()
