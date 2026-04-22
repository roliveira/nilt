"""2D advection-diffusion: instantaneous point release in a uniform flow.

An infinite 2D domain with uniform velocity v in the x-direction and
isotropic diffusion coefficient D. A mass M [kg/m] is released
instantaneously at the origin at t=0. The resulting plume is elongated
downstream - there is NO radial symmetry.

PDE:
  dC/dt + v * dC/dx = D * (d^2C/dx^2 + d^2C/dy^2)

Laplace domain (via Galilean substitution C_bar = exp(vx/(2D)) * psi):
  C_bar(x,y,s) = M / (2*pi*D) * exp(v*x/(2*D))
                  * K0(r * sqrt(v^2/(4D^2) + s/D))

where r = sqrt(x^2 + y^2) and K0 is the modified Bessel function of the
second kind, order zero.

Analytical:
  C(x,y,t) = M / (4*pi*D*t) * exp(-((x - v*t)^2 + y^2) / (4*D*t))

Parameters:
  D = 0.1 m^2/s  (effective isotropic dispersion)
  v = 0.5 m/s   (uniform flow velocity in x)
  M = 1.0 kg/m  (mass released)

Reference:
  Bear, J. (1972). Dynamics of Fluids in Porous Media, S10.6,
  American Elsevier.
"""

import os
import numpy as np
from scipy.special import kv
import nilt

# Parameters
D = 0.1    # dispersion [m^2/s]
v = 0.5    # flow velocity in x [m/s]
M = 1.0    # released mass [kg/m]

def K0(z):
    return kv(0, z)

# Observation point for time series
xp, yp = 3.0, 0.5
rp = np.hypot(xp, yp)

def Fs_pt(s):
    kappa = np.sqrt(v**2 / (4 * D**2) + s / D)
    return M / (2 * np.pi * D) * np.exp(v * xp / (2 * D)) * K0(rp * kappa)

def analytical_pt(t):
    return M / (4 * np.pi * D * t) * np.exp(-((xp - v * t)**2 + yp**2) / (4 * D * t))

talbot = nilt.Talbot()
dehoog = nilt.DeHoog()

t = np.geomspace(0.5, 30.0, 80)
ana = analytical_pt(t)
ta  = nilt.invert(talbot, Fs_pt, t)
dh  = nilt.invert(dehoog, Fs_pt, t)

out_dir = os.path.join(os.path.dirname(__file__), "build")
os.makedirs(out_dir, exist_ok=True)

out = os.path.join(out_dir, "py_transport_advection_plume_2d.csv")
np.savetxt(out, np.column_stack([t, ana, ta, dh]),
           delimiter=",", header="t,analytical,talbot,dehoog", comments="")
print(out)

# 2D spatial field at t=5s
t_snap = 5.0
x_lo, x_hi = -3.0, 8.0
y_lo, y_hi = -3.0, 3.0
nx, ny = 140, 80
xs = np.linspace(x_lo, x_hi, nx)
ys = np.linspace(y_lo, y_hi, ny)

data = []
for j in range(ny):
    for i in range(nx):
        xx, yy = xs[i], ys[j]
        rr = np.hypot(xx, yy)
        anal = M / (4 * np.pi * D * t_snap) * np.exp(-((xx - v * t_snap)**2 + yy**2) / (4 * D * t_snap))
        if rr < 1e-12:
            nilt_val = anal
        else:
            def Fs_xy(s, _xx=xx, _rr=rr):
                kappa = np.sqrt(v**2 / (4 * D**2) + s / D)
                return M / (2 * np.pi * D) * np.exp(v * _xx / (2 * D)) * K0(_rr * kappa)
            nilt_val = nilt.invert(talbot, Fs_xy, t_snap)
        data.append([xx, yy, anal, nilt_val])

out2 = os.path.join(out_dir, "py_transport_advection_plume_2d_spatial.csv")
np.savetxt(out2, np.array(data),
           delimiter=",", header="x,y,analytical,talbot", comments="")
print(out2)
