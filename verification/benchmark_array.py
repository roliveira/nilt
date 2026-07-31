"""Array-size timing benchmark - measures wall-clock time for whole-array
inversions as len(t) is varied.  Algorithm parameters are kept at their
defaults (Stehfest N=18, Talbot N=50, DeHoog M=40).

Test function: F(s) = 1/(s+1),  f(t) = exp(-t).

Output: py_benchmark_array.csv with columns
  method, nt, time_us  (microseconds per whole-array inversion)
"""

import time

import numpy as np
import nilt

WARMUP  = 20
REPEATS = 100

def Fs(s):
    return 1.0 / (s + 1.0)


def time_array(algo, t):
    for _ in range(WARMUP):
        algo(Fs, t)

    t0 = time.perf_counter()
    for _ in range(REPEATS):
        algo(Fs, t)
    t1 = time.perf_counter()
    return (t1 - t0) / REPEATS * 1e6


SIZES = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]

out = "py_benchmark_array.csv"
with open(out, "w") as f:
    f.write("method,nt,time_us\n")

    for name, algo in [("Stehfest", nilt.Stehfest()),
                       ("Talbot",   nilt.Talbot()),
                       ("DeHoog",   nilt.DeHoog())]:
        for nt in SIZES:
            t = np.linspace(0.1, 5.0, nt) if nt > 1 else np.array([1.0])
            us = time_array(algo, t)
            f.write(f"{name},{nt},{us:.6f}\n")
            print(f"{name:9s} nt={nt:6d}  {us:.1f} us")

print(out)
