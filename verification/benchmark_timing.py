"""Timing benchmark - measures wall-clock time per inversion as the main
tuning parameter is varied for each algorithm.

Uses func4: F(s) = 1/(s+1), f(t) = exp(-t) as the test function (cheap to
evaluate so the timings reflect the algorithm cost, not the user function).

Output: benchmark_timing_python.csv with columns
  method, param, time_us  (microseconds per single inversion at t=1)
"""

import os
import time

import numpy as np
import nilt

WARMUP  = 50
REPEATS = 500
T_EVAL  = 1.0

# Benchmark function: F(s) = 1/(s+1)
def Fs_real(s):
    return 1.0 / (s + 1.0)

def Fs_cplx(s):
    return 1.0 / (s + 1.0)


def time_inversion(algo, Fs):
    """Return average time in microseconds for a single inversion."""
    for _ in range(WARMUP):
        nilt.invert(algo, Fs, T_EVAL)

    t0 = time.perf_counter()
    for _ in range(REPEATS):
        nilt.invert(algo, Fs, T_EVAL)
    t1 = time.perf_counter()

    return (t1 - t0) / REPEATS * 1e6  # seconds -> microseconds


out_dir = os.path.join(os.path.dirname(__file__), "build")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "benchmark_timing_python.csv")

with open(out, "w") as f:
    f.write("method,param,time_us\n")

    # Stehfest: vary N (must be even)
    for N in [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]:
        algo = nilt.Stehfest()
        algo.N = N
        us = time_inversion(algo, Fs_real)
        f.write(f"Stehfest,{N},{us:.6f}\n")
        print(f"Stehfest  N={N:3d}  {us:.1f} us")

    # Talbot: vary n
    for n in [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200]:
        algo = nilt.Talbot()
        algo.n = n
        us = time_inversion(algo, Fs_cplx)
        f.write(f"Talbot,{n},{us:.6f}\n")
        print(f"Talbot    n={n:3d}  {us:.1f} us")

    # DeHoog: vary M
    for M in [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200]:
        algo = nilt.DeHoog()
        algo.M = M
        us = time_inversion(algo, Fs_cplx)
        f.write(f"DeHoog,{M},{us:.6f}\n")
        print(f"DeHoog    M={M:3d}  {us:.1f} us")

print(out)
