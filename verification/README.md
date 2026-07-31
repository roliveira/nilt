# Verification - Standard Test Functions

Evaluates all three inversion algorithms against 10 Laplace transform pairs
with known analytical solutions.

## Test Functions

|  # | f(t)              | F(s)                             | Source               |
|----|-------------------|----------------------------------|----------------------|
|  1 | $1/\sqrt{\pi t}$  | $1/\sqrt{s}$                     | Stehfest (1970)      |
|  2 | $-\gamma - \ln t$ | $\ln(s)/s$                       | Stehfest (1970)      |
|  3 | $t^3/6$           | $s^{-4}$                         | Stehfest (1970)      |
|  4 | $e^{-t}$          | $1/(s+1)$                        | Stehfest (1970)      |
|  5 | $\sin\sqrt{2t}$   | $\sqrt{\pi/(2s^3)}\,e^{-1/(2s)}$ | Stehfest (1970)      |
|  6 | $t$               | $1/s^2$                          | Abate & Whitt (2006) |
|  7 | $t\,e^{-t}$       | $1/(s+1)^2$                      | Abate & Whitt (2006) |
|  8 | $\sin t$          | $1/(s^2+1)$                      | Abate & Whitt (2006) |
|  9 | $\cos t$          | $s/(s^2+1)$                      | Abate & Whitt (2006) |
| 10 | $e^{-t}\sin t$    | $1/((s+1)^2+1)$                  | Abate & Whitt (2006) |

## Output

Produces 30 CSV files (10 functions x 3 methods), each with columns
`t, fta, ftn, err` for t = 1, 2, ..., 10.

## Running

```bash
# Build (from the repo root)
cmake -B build -DNILT_BUILD_VERIFICATION=ON
cmake --build build

# Run from build/ directory (all output goes to cwd)
cd build
./verification
python ../verification/verification.py

# Timing benchmarks
./benchmark_timing         # scalar-t sweep over N/M       -> cpp_benchmark_timing.csv
python ../verification/benchmark_timing.py                  # ->  py_benchmark_timing.csv

./benchmark_array          # array-length sweep at defaults -> cpp_benchmark_array.csv
python ../verification/benchmark_array.py                   # ->  py_benchmark_array.csv

python ../verification/plot_benchmark.py  # -> benchmark_timing.png, benchmark_array.png
```

The array benchmark uses the batched Python entry point (`algo(Fs, t_array)`),
which invokes `Fs` exactly once per call - see the main README for details.

## References

- Stehfest, H. (1970). *Algorithm 368: Numerical inversion of Laplace
  transforms*. Commun. ACM 13(1), 47-49.
- Abate, J. & Whitt, W. (2006). *A unified framework for numerically inverting
  Laplace transforms*. INFORMS J. Comput. 18(4), 408-421.
