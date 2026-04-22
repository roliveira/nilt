"""Theis well function - drawdown from a pumping well.

An infinite confined aquifer, initially at head h0, is pumped at constant
rate Q from a fully-penetrating well at r = 0 starting at t = 0.
We compute the drawdown s(r,t) = h0 - h(r,t).

Two outputs are produced:
  1. Drawdown vs time at a fixed observation distance (time series).
  2. Drawdown vs distance at a fixed snapshot time (distance profile).

Equation (Theis, 1935 - eq. 5, US customary units):
  s = (114.6 Q / T) * W(u),   u = 1.87 r^2 S / (T t)

where W(u) = E_1(u) is the exponential integral (well function).

Laplace domain (consistent ft/day units):
  s_bar(r,p) = Q / (2*pi*T*p) * K0(r * sqrt(p*S/T))

Parameters (from Theis, 1935, Fig. 1):
  Q = 525 gal/min    (pumping rate)
  T = 90,000 gal/day/ft (transmissibility)
  S = 0.2225          (specific yield)
  r = 100 ft          (observation distance, time series)
  t = 2 days          (snapshot time, distance profile)

Reference:
  Theis, C.V. (1935). The relation between the lowering of the
  piezometric surface and the rate and duration of discharge of a well
  using ground-water storage. Trans. Am. Geophys. Union, 16(2), 519-524.
"""

import os
import numpy as np
from scipy.special import exp1, kv
import nilt

# Parameters in US customary units (Theis, 1935)
Q_gpm = 525.0          # pumping rate [gal/min]
T_gpd = 90_000.0       # transmissibility [gal/day/ft]
S     = 0.2225          # specific yield [-]

# Unit conversions to consistent ft/day system
_GPM_TO_FT3DAY = 1440.0 * 0.133681
_GPD_TO_FT2DAY = 0.133681

Q     = Q_gpm * _GPM_TO_FT3DAY   # ft^3/day
T     = T_gpd * _GPD_TO_FT2DAY   # ft^2/day
alpha = T / S                     # ft^2/day


def K0(z):
    return kv(0, z)


def analytical(r_ft, t_days):
    u = 1.87 * r_ft**2 * S / (T_gpd * t_days)
    return 114.6 * Q_gpm / T_gpd * exp1(u)


stehfest = nilt.Stehfest()
talbot   = nilt.Talbot()
dehoog   = nilt.DeHoog()

out_dir = os.path.join(os.path.dirname(__file__), "build")
os.makedirs(out_dir, exist_ok=True)

r_obs = 100.0  # observation distance [ft]


def Fs_time(s):
    return Q / (2.0 * np.pi * T * s) * K0(r_obs * np.sqrt(s * S / T))


t = np.geomspace(0.01, 100.0, 80)

ana = analytical(r_obs, t)
st  = nilt.invert(stehfest, Fs_time, t)
ta  = nilt.invert(talbot,   Fs_time, t)
dh  = nilt.invert(dehoog,   Fs_time, t)

out = os.path.join(out_dir, "py_groundwater_theis_well_time.csv")
np.savetxt(out, np.column_stack([t, ana, st, ta, dh]),
           delimiter=",", header="t,analytical,stehfest,talbot,dehoog", comments="")
print(out)

t_snap = 2.0  # snapshot time [days] (48 h)

r_arr = np.linspace(1.0, 1000.0, 300)
ana_r = analytical(r_arr, t_snap)

st_r = np.empty_like(r_arr)
ta_r = np.empty_like(r_arr)
dh_r = np.empty_like(r_arr)
for i, ri in enumerate(r_arr):
    def Fs_r(s, _r=ri):
        return Q / (2.0 * np.pi * T * s) * K0(_r * np.sqrt(s * S / T))
    st_r[i] = nilt.invert(stehfest, Fs_r, t_snap)
    ta_r[i] = nilt.invert(talbot,   Fs_r, t_snap)
    dh_r[i] = nilt.invert(dehoog,   Fs_r, t_snap)

out2 = os.path.join(out_dir, "py_groundwater_theis_well_distance.csv")
np.savetxt(out2, np.column_stack([r_arr, ana_r, st_r, ta_r, dh_r]),
           delimiter=",", header="r,analytical,stehfest,talbot,dehoog", comments="")
print(out2)
