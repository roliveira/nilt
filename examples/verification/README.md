# Verification - Standard Test Functions

Evaluates all three inversion algorithms against 10 Laplace transform pairs
with known analytical solutions.

## Test Functions

| # | f(t) | F(s) | Source |
|---|------|------|--------|
| 1 | $1/\sqrt{\pi t}$ | $1/\sqrt{s}$ | Stehfest (1970) |
| 2 | $-\gamma - \ln t$ | $\ln(s)/s$ | Stehfest (1970) |
| 3 | $t^3/6$ | $s^{-4}$ | Stehfest (1970) |
| 4 | $e^{-t}$ | $1/(s+1)$ | Standard |
| 5 | $\sin\sqrt{2t}$ | $\sqrt{\pi/(2s^3)}\,e^{-1/(2s)}$ | Stehfest (1970) |
| 6 | $t$ | $1/s^2$ | Abate & Whitt (2006) |
| 7 | $t\,e^{-t}$ | $1/(s+1)^2$ | Abate & Whitt (2006) |
| 8 | $\sin t$ | $1/(s^2+1)$ | Abate & Whitt (2006) |
| 9 | $\cos t$ | $s/(s^2+1)$ | Abate & Whitt (2006) |
| 10 | $e^{-t}\sin t$ | $1/((s+1)^2+1)$ | Abate & Whitt (2006) |

## Output

Produces 30 CSV files (10 functions x 3 methods), each with columns
`t, fta, ftn, err` for t = 1, 2, ..., 10.

## Running

```bash
# from the examples/verification/build directory
./verification

# plot results (reads/writes from build/ by default)
python ../plot_verification.py
python ../plot_benchmark.py
```

## References

- Stehfest, H. (1970). *Algorithm 368: Numerical inversion of Laplace
  transforms*. Commun. ACM 13(1), 47-49.
- Abate, J. & Whitt, W. (2006). *A unified framework for numerically inverting
  Laplace transforms*. INFORMS J. Comput. 18(4), 408-421.
