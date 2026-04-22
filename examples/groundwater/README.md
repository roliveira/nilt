# Groundwater - Theis Well Function

## Drawdown from a pumping well (`theis_well.cpp`)

Models drawdown in a confined aquifer pumped at constant rate $Q$ from a
fully-penetrating well, using the parameters from the original Theis (1935)
publication. Two outputs are produced:

1. **Time series** - drawdown vs time at a fixed observation distance ($r = 100$ ft).
2. **Distance profile** - drawdown vs distance at a fixed snapshot time ($t = 2$ days).

**Laplace domain (consistent ft/day units):**

$$\bar{s}(r,p) = \frac{Q}{2\pi T\,p}\,K_0\!\left(r\sqrt{\frac{pS}{T}}\right)$$

where $T$ is transmissibility, $S$ is specific yield, and $K_0$ is the
modified Bessel function of the second kind, order zero.

**Analytical (Theis, 1935, eq. 5 - US customary units):**

$$s = \frac{114.6\,Q}{T}\,W(u), \qquad u = \frac{1.87\,r^2 S}{Tt}$$

where $W(u) = E_1(u)$ is the exponential integral (well function).

**Parameters:** $Q = 525$ gal/min, $T = 90{,}000$ gal/day/ft, $S = 0.2225$.

**Reference:** Theis, C.V. (1935). Trans. Am. Geophys. Union, 16(2), 519-524.

## 2. Pumping-injection well dipole - **asymmetric** (`well_dipole.cpp`)

Two wells at different locations with different rates ($Q_p \neq Q_i$): one pumping,
one injecting. The resulting drawdown field has no radial symmetry.

**Laplace domain (superposition):**

$$\bar{s}(x,y,p) = \frac{Q_p}{2\pi Tp}\,K_0\!\left(r_p\sqrt{\frac{p}{\alpha}}\right) - \frac{Q_i}{2\pi Tp}\,K_0\!\left(r_i\sqrt{\frac{p}{\alpha}}\right)$$

**Analytical:**

$$s(x,y,t) = \frac{Q_p}{4\pi T}\,E_1(u_p) - \frac{Q_i}{4\pi T}\,E_1(u_i), \qquad u = \frac{r^2 S}{4Tt}$$

**Parameters:** $T = 10^{-3}$ m$^2$/s, $S = 10^{-4}$.
Pumping well: $Q_p = 0.010$ m$^3$/s at ($-15$, $-5$) m.
Injection well: $Q_i = 0.007$ m$^3$/s at ($+20$, $+8$) m.

**Reference:** Bear (1979). *Hydraulics of Groundwater*, §8.3.

## Running

```bash
# from the examples/groundwater/build directory
./groundwater_theis_well
./groundwater_well_dipole

# plot results (reads/writes from build/ by default)
python ../plot_theis_well.py
python ../plot_well_dipole.py
```
