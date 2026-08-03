"""Timing benchmark - measures wall-clock time per inversion as the main
tuning parameter is varied for each algorithm.

Uses func4: F(s) = 1/(s+1), f(t) = exp(-t) as the test function (cheap to
evaluate so the timings reflect the algorithm cost, not the user function).

Output: py_benchmark_timing.csv with columns method, param, time_us.
"""

import csv
import time

import nilt

WARMUP  = 50
REPEATS = 500
T_EVAL  = 1.0

# Benchmark function: F(s) = 1/(s+1)
def Fs_real(s):
    return 1.0 / (s + 1.0)

def Fs_cplx(s):
    return 1.0 / (s + 1.0)


STEHFEST_NS = [4, 6, 8, 10, 12, 14, 16, 18, 20]
TALBOT_NS   = [1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
               64, 68, 72, 76, 80, 84, 88, 92, 96, 100]
DEHOOG_MS   = [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200]


def time_inversion(algo, Fs):
    """Return average time in microseconds for a single inversion."""
    for _ in range(WARMUP):
        algo(Fs, T_EVAL)

    t0 = time.perf_counter()
    for _ in range(REPEATS):
        algo(Fs, T_EVAL)
    t1 = time.perf_counter()

    return (t1 - t0) / REPEATS * 1e6  # seconds -> microseconds


out = "py_benchmark_timing.csv"

with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["method", "param", "time_us"])

    for N in STEHFEST_NS:
        algo = nilt.Stehfest(); algo.N = N
        us = time_inversion(algo, Fs_real)
        w.writerow(["Stehfest", N, f"{us:.6f}"])
        print(f"Stehfest  N={N:3d}  {us:.1f} us")

    for N in TALBOT_NS:
        algo = nilt.Talbot(); algo.N = N
        us = time_inversion(algo, Fs_cplx)
        w.writerow(["Talbot", N, f"{us:.6f}"])
        print(f"Talbot    N={N:3d}  {us:.1f} us")

    for M in DEHOOG_MS:
        algo = nilt.DeHoog(); algo.M = M
        us = time_inversion(algo, Fs_cplx)
        w.writerow(["DeHoog", M, f"{us:.6f}"])
        print(f"DeHoog    M={M:3d}  {us:.1f} us")

print(out)

