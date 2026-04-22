#!/usr/bin/env python3
"""Plot advection-diffusion plume results (time series and 2D spatial)."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser


def main() -> None:
    args = create_argument_parser(__doc__, __file__).parse_args()

    csv_path = args.dir / "transport_advection_plume_2d.csv"
    if csv_path.exists():
        df = pl.read_csv(csv_path)
        t = df["t"].to_numpy()
        ana = df["analytical"].to_numpy()

        fig, (ax, ax_err) = plt.subplots(
            2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]},
            sharex=True, constrained_layout=True,
        )

        ax.plot(t, ana, "k-", label="analytical", zorder=10)
        ax.plot(t, df["talbot"], "^", color="C1", markevery=5, label="Talbot")
        ax.plot(t, df["dehoog"], "o", color="C2", markevery=5, label="De Hoog")
        ax.set_ylabel("Concentration [kg/m³]")
        ax.set_title(r"$C(x,y,t) = \frac{M}{4\pi Dt}\,e^{-[(x-vt)^2+y^2]/4Dt}$")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ana_abs = np.maximum(np.abs(ana), 1e-30)
        for col, marker, color in [("talbot", "^", "C1"), ("dehoog", "o", "C2")]:
            rel_err = np.abs(df[col].to_numpy() - ana) / ana_abs
            ax_err.semilogy(t, rel_err, marker=marker, color=color, markevery=5)

        ax_err.set_xlabel("Time [s]")
        ax_err.set_ylabel("Relative error")
        ax_err.grid(True, alpha=0.3, which="both")

        out = csv_path.with_suffix(f".{args.format}")
        fig.savefig(out, dpi=600)
        plt.close(fig)
        print(f"  {out}")

    spatial = args.dir / "transport_advection_plume_2d_spatial.csv"
    if spatial.exists():
        ds = pl.read_csv(spatial, infer_schema_length=10000)
        xs = np.sort(ds["x"].unique().to_numpy())
        ys = np.sort(ds["y"].unique().to_numpy())
        nx, ny = len(xs), len(ys)
        Z = ds.sort(["y", "x"])["talbot"].to_numpy().reshape(ny, nx)
        X, Y = np.meshgrid(xs, ys)

        fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
        cf = ax.contourf(X, Y, Z, levels=30)
        cb = fig.colorbar(cf, ax=ax)
        cb.set_label("Concentration [kg/m³]")
        ax.plot(0, 0, "wx", label="Release point")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(r"$C = \frac{M}{4\pi Dt}\,e^{-[(x-vt)^2+y^2]/4Dt}$ at $t = 5$ s")
        ax.set_aspect("equal")
        ax.legend(loc="upper left")

        out = spatial.with_name("transport_advection_plume_2d_spatial").with_suffix(f".{args.format}")
        fig.savefig(out, dpi=600)
        plt.close(fig)
        print(f"  {out}")


if __name__ == "__main__":
    main()
