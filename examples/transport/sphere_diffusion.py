"""Diffusion from a sphere into an infinite medium.

A solid sphere of radius a, initially at concentration C0, is placed in
an infinite medium at concentration 0 at t=0. We track the average
concentration inside the sphere.

PDE (radial diffusion):
  dC/dt = D * (1/r^2) * d/dr(r^2 * dC/dr),  0 < r < a

Laplace domain (average remaining concentration):
  C_avg_bar(s) = (C0/s) * (1 - (3/q)*(coth(q) - 1/q)),  q = a*sqrt(s/D)

Analytical (series solution):
  C_avg(t) = C0 * (6/pi^2) * sum_{n=1}^{inf} (1/n^2) * exp(-D*n^2*pi^2*t/a^2)

Parameters:
  D = 1e-9 m^2/s (typical ionic diffusion in water)
  a = 1e-3 m (1 mm radius)
  C0 = 1.0 mol/m^3

Reference:
  Crank, J. (1975). The Mathematics of Diffusion, 2nd ed., S6.3,
  Oxford University Press.
"""

import os
import numpy as np
import nilt

# Parameters
D  = 1.0e-9   # diffusivity [m^2/s]
a  = 1.0e-3   # sphere radius [m]
C0 = 1.0      # initial concentration [mol/m^3]

def Fs(s):
    q = a * np.sqrt(s / D)
    coth_q = np.cosh(q) / np.sinh(q)
    return C0 / s * (1.0 - (3.0 / q) * (coth_q - 1.0 / q))

def analytical(t):
    result = np.zeros_like(t, dtype=float)
    for n in range(1, 501):
        n2 = n * n
        result += np.exp(-D * n2 * np.pi**2 * t / a**2) / n2
    return C0 * 6.0 / np.pi**2 * result

stehfest = nilt.Stehfest()
talbot   = nilt.Talbot()
dehoog   = nilt.DeHoog()

t = np.geomspace(10.0, 5000.0, 80)

ana = analytical(t)
st  = nilt.invert(stehfest, Fs, t)
ta  = nilt.invert(talbot,   Fs, t)
dh  = nilt.invert(dehoog,   Fs, t)

out = os.path.join(os.path.dirname(__file__), "build", "py_transport_sphere_diffusion.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
np.savetxt(out, np.column_stack([t, ana, st, ta, dh]),
           delimiter=",", header="t,analytical,stehfest,talbot,dehoog", comments="")
print(out)
