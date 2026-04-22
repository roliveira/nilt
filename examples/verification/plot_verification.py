#!/usr/bin/env python3
"""Plot verification results for numerical inverse Laplace transform methods."""

import pathlib
import sys

import polars as pl
import matplotlib.pyplot as plt

METHOD_NAMES = {"Stehfest": "Stehfest", "Talbot": "Talbot", "DeHoog": "De Hoog"}
METHOD_STYLES = {
    "Stehfest": ("s", "C0"),
    "Talbot":   ("^", "C1"),
    "DeHoog":   ("o", "C2"),
}
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

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "utils"))
from plot_args import create_argument_parser as _base_parser


def create_argument_parser():
    parser = _base_parser(__doc__, __file__)
    parser.add_argument(
        "--functions",
        type=int,
        nargs="+",
        choices=range(1, 11),
        default=range(1, 11),
        metavar="N",
        help="Test function numbers to plot (1-10, default: all)",
    )
    return parser


def main() -> None:
    args = create_argument_parser().parse_args()

    for func in args.functions:
        fig, (ax_val, ax_err) = plt.subplots(
            2, 1, figsize=(7, 5), gridspec_kw={"height_ratios": [3, 1]}, sharex=True,
            constrained_layout=True,
        )

        plotted_analytical = False
        for method_key, label in METHOD_NAMES.items():
            csv_path = args.dir / f"verification_func{func}_{method_key}.csv"
            if not csv_path.exists():
                continue

            df = pl.read_csv(csv_path)
            marker, color = METHOD_STYLES[method_key]

            if not plotted_analytical:
                ax_val.plot(df["t"], df["fta"], "k-", label="analytical", zorder=10)
                plotted_analytical = True

            ax_val.plot(df["t"], df["ftn"], marker=marker, color=color, label=label)
            ax_err.semilogy(df["t"], df["err"], marker=marker, color=color, label=label)

        ax_val.set_title(FUNC_LABELS[func])
        ax_val.set_ylabel(r"$f(t)$")
        ax_val.legend(loc="best")
        ax_val.grid(True, alpha=0.3)

        ax_err.set_xlabel(r"$t$")
        ax_err.set_ylabel("Relative error")
        ax_err.grid(True, alpha=0.3, which="both")

        out = args.dir / f"verification_func{func}.{args.format}"
        fig.savefig(out, dpi=600)
        plt.close(fig)
        print(f"  {out}")


if __name__ == "__main__":
    main()
