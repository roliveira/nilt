"""Axisymmetric diffusion from a long circular cylinder into an infinite medium.

A long solid cylinder of radius a, initially at concentration C0, is placed
in an infinite medium at concentration 0 at t=0. We track the average
concentration inside the cylinder.

PDE (cylindrical coordinates, axisymmetric, no z-dependence):
  dC/dt = D * (1/r) * d/dr(r * dC/dr),  0 < r < a

Laplace domain (average remaining concentration):
  C_avg_bar(s) = C0/s * (1 - (2/q) * I1(q) / I0(q)),   q = a*sqrt(s/D)

where I0, I1 are modified Bessel functions of the first kind.

Analytical (series solution):
  C_avg(t) = C0 * 4 * sum_{n=1}^{inf} exp(-D*alpha_n^2*t/a^2) / (a^2 * alpha_n^2)

where alpha_n are the positive roots of J0(alpha) = 0 (Bessel zeros).

For simplicity we use the first ~40 Bessel zeros for the series.

Parameters:
  D = 5e-10 m^2/s (slower diffusion, e.g. large molecule in gel)
  a = 5e-4 m (0.5 mm radius)
  C0 = 1.0 mol/m^3

Reference:
  Crank, J. (1975). The Mathematics of Diffusion, 2nd ed., S5.3,
  Oxford University Press.
"""

import os
import numpy as np
from scipy.special import iv, j0, jn_zeros
import nilt

# Parameters
D  = 5.0e-10   # diffusivity [m^2/s]
a  = 5.0e-4    # cylinder radius [m]
C0 = 1.0       # initial concentration [mol/m^3]

def I0(z):
    return iv(0, z)

def I1(z):
    return iv(1, z)

def Fs(s):
    q = a * np.sqrt(s / D)
    return C0 / s * (1.0 - (2.0 / q) * I1(q) / I0(q))

# First 40 positive zeros of J0
j0_zeros = jn_zeros(0, 40)

def analytical(t):
    result = np.zeros_like(t, dtype=float)
    for an in j0_zeros:
        result += np.exp(-D * an**2 * t / a**2) / an**2
    return C0 * 4.0 * result

stehfest = nilt.Stehfest()
talbot   = nilt.Talbot()
dehoog   = nilt.DeHoog()

# Stehfest needs real-valued F(s)
def Fs_real(s):
    q = a * np.sqrt(s / D)
    return float(np.real(C0 / s * (1.0 - (2.0 / q) * I1(q) / I0(q))))

t = np.geomspace(10.0, 5000.0, 80)

ana = analytical(t)
st  = nilt.invert(stehfest, Fs_real, t)
ta  = nilt.invert(talbot,   Fs, t)
dh  = nilt.invert(dehoog,   Fs, t)

out = os.path.join(os.path.dirname(__file__), "build", "py_transport_cylinder_diffusion.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
np.savetxt(out, np.column_stack([t, ana, st, ta, dh]),
           delimiter=",", header="t,analytical,stehfest,talbot,dehoog", comments="")
print(out)
