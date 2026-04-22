"""Pumping-injection well dipole - asymmetric 2D groundwater flow.

A confined aquifer with two wells at different locations:
  - Pumping well at (x_p, y_p) extracting at rate Q_p
  - Injection well at (x_i, y_i) injecting at rate Q_i

The net drawdown is the superposition of individual Theis solutions.
When Q_p != Q_i or wells are not symmetric about the origin, the
drawdown field has NO radial symmetry.

Laplace domain (superposition):
  s_bar(x,y,p) = Q_p / (2*pi*T*p) * K0(r_p * sqrt(p/alpha))
               - Q_i / (2*pi*T*p) * K0(r_i * sqrt(p/alpha))

where r_p, r_i are distances to the pumping and injection wells,
and the negative sign on Q_i represents head recovery from injection.

Analytical (superposition of Theis):
  s(x,y,t) = Q_p/(4*pi*T) * E1(r_p^2*S/(4Tt))
           - Q_i/(4*pi*T) * E1(r_i^2*S/(4Tt))

Parameters:
  T = 1e-3 m^2/s, S = 1e-4
  Well 1 (pumping):   Q_p = 0.010 m^3/s  at (-15, -5) m
  Well 2 (injection): Q_i = 0.007 m^3/s  at (+20, +8) m

Reference:
  Bear, J. (1979). Hydraulics of Groundwater, S8.3, McGraw-Hill.
"""

import os
import numpy as np
from scipy.special import kv, exp1
import nilt

# Parameters
T_aq  = 1.0e-3   # transmissivity [m^2/s]
S     = 1.0e-4   # storativity
alpha = T_aq / S

Qp, xp, yp = 0.010, -15.0, -5.0    # pumping well
Qi, xi, yi  = 0.007,  20.0,  8.0    # injection well

def K0(z):
    return kv(0, z)

# Observation point
xo, yo = 5.0, 0.0
rp_obs = np.hypot(xo - xp, yo - yp)
ri_obs = np.hypot(xo - xi, yo - yi)

def Fs_obs(s):
    sa = np.sqrt(s / alpha)
    return (Qp / (2 * np.pi * T_aq * s) * K0(rp_obs * sa)
          - Qi / (2 * np.pi * T_aq * s) * K0(ri_obs * sa))

def analytical_obs(t):
    up = rp_obs**2 * S / (4.0 * T_aq * t)
    ui = ri_obs**2 * S / (4.0 * T_aq * t)
    return Qp / (4 * np.pi * T_aq) * exp1(up) - Qi / (4 * np.pi * T_aq) * exp1(ui)

talbot = nilt.Talbot()
dehoog = nilt.DeHoog()

t = np.geomspace(10.0, 100000.0, 80)
ana = analytical_obs(t)
ta  = nilt.invert(talbot, Fs_obs, t)
dh  = nilt.invert(dehoog, Fs_obs, t)

out_dir = os.path.join(os.path.dirname(__file__), "build")
os.makedirs(out_dir, exist_ok=True)

out = os.path.join(out_dir, "py_groundwater_well_dipole.csv")
np.savetxt(out, np.column_stack([t, ana, ta, dh]),
           delimiter=",", header="t,analytical,talbot,dehoog", comments="")
print(out)

# 2D spatial drawdown at t=7200s
t_snap = 7200.0
half   = 60.0
nx, ny = 120, 120
xs = np.linspace(-half, half, nx)
ys = np.linspace(-half, half, ny)

data = []
for j in range(ny):
    for i in range(nx):
        xx, yy = xs[i], ys[j]
        dp = np.hypot(xx - xp, yy - yp)
        di = np.hypot(xx - xi, yy - yi)
        up = dp**2 * S / (4 * T_aq * t_snap)
        ui = di**2 * S / (4 * T_aq * t_snap)
        anal = Qp / (4 * np.pi * T_aq) * exp1(up) - Qi / (4 * np.pi * T_aq) * exp1(ui)

        def Fs_xy(s, _dp=dp, _di=di):
            sa = np.sqrt(s / alpha)
            return (Qp / (2 * np.pi * T_aq * s) * K0(_dp * sa)
                  - Qi / (2 * np.pi * T_aq * s) * K0(_di * sa))
        nilt_val = nilt.invert(talbot, Fs_xy, t_snap)
        data.append([xx, yy, anal, nilt_val])

out2 = os.path.join(out_dir, "py_groundwater_well_dipole_spatial.csv")
np.savetxt(out2, np.array(data),
           delimiter=",", header="x,y,analytical,talbot", comments="")
print(out2)
