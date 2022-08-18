#!/usr/bin/env python
import sys
import pathlib
import argparse

import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input",     required=True)
    parser.add_argument("-o", "--output",    required=False, default="benchmark")
    parser.add_argument("-m", "--methods",   required=True,  nargs="+", choices=["talbot", "stehfest", "dehoog"])
    
    parser.add_argument("-e", "--extension", required=False, default=".csv")
    parser.add_argument("-p", "--prefix",    required=False, default="benchmark")
    parser.add_argument("-f", "--functions", required=False, choices=range(1, 6), default=range(1, 6))

    args = parser.parse_args(sys.argv[1:])

    dfs = {}

    for imethod in args.methods:
        dfs[imethod] = {}

        for ifunc in args.functions:
            fname = "_".join([args.prefix, imethod, "func" + str(ifunc)]) + args.extension
            fname = pathlib.Path(args.input).joinpath(fname)

            dfs[imethod][ifunc] = pd.read_csv(fname)

    for ifunc in args.functions:
        fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

        for i, imethod in enumerate(dfs.keys()):
            if i == 0:
                _ = ax[0].plot(dfs[imethod][ifunc]["t"], dfs[imethod][ifunc]["fta"], "-k", label="analytical")

            _ = ax[0].plot(dfs[imethod][ifunc]["t"], dfs[imethod][ifunc]["ftn"], "-o", label=imethod)
            _ = ax[1].plot(dfs[imethod][ifunc]["t"], dfs[imethod][ifunc]["err"], "-o", label=imethod)

        _ = ax[0].set_title("func" + str(ifunc))
        _ = ax[0].legend(loc="upper right")
        
        _ = ax[0].set_ylabel("f(t)")
        _ = ax[1].set_ylabel("error")
        _ = ax[1].set_xlabel("t")

        _ = fig.savefig(args.output + "_func" + str(ifunc) + ".pdf")
        _ = fig.savefig(args.output + "_func" + str(ifunc) + ".png", dpi=300)
