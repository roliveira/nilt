"""Verification - standard test functions for numerical inverse Laplace transforms.

Evaluates all three inversion algorithms (Stehfest, Talbot, De Hoog) against
10 Laplace transform pairs with known analytical solutions. The first five
are from Stehfest (1970) and the remaining five from Abate & Whitt (2006).

Output: one CSV per (function, method) combination, containing columns
  t, fta (analytical), ftn (numerical), err (relative error).

References:
  Stehfest, H. (1970). Commun. ACM 13(1), 47-49.
  Abate, J. & Whitt, W. (2006). INFORMS J. Comput. 18(4), 408-421.
"""

import os
import numpy as np
import nilt

pi = np.pi

# Test functions: (name, F(s), f(t))
# Stehfest 1970
def Fs1(s): return 1.0 / np.sqrt(s)
def ft1(t): return 1.0 / np.sqrt(pi * t)

def Fs2(s): return np.log(s) / s
def ft2(t): return -0.57722 - np.log(t)

def Fs3(s): return s ** (-4.0)
def ft3(t): return t ** 3 / 6.0

def Fs4(s): return 1.0 / (s + 1.0)
def ft4(t): return np.exp(-t)

def Fs5(s): return np.sqrt(pi / (2.0 * s ** 3)) * np.exp(-1.0 / (2.0 * s))
def ft5(t): return np.sin(np.sqrt(2.0 * t))

# Abate & Whitt 2006
def Fs6(s): return 1.0 / (s * s)
def ft6(t): return t

def Fs7(s): return 1.0 / ((s + 1.0) ** 2)
def ft7(t): return t * np.exp(-t)

def Fs8(s): return 1.0 / (s * s + 1.0)
def ft8(t): return np.sin(t)

def Fs9(s): return s / (s * s + 1.0)
def ft9(t): return np.cos(t)

def Fs10(s): return 1.0 / ((s + 1.0) ** 2 + 1.0)
def ft10(t): return np.exp(-t) * np.sin(t)

functions = [
    ("func1",  Fs1,  ft1),
    ("func2",  Fs2,  ft2),
    ("func3",  Fs3,  ft3),
    ("func4",  Fs4,  ft4),
    ("func5",  Fs5,  ft5),
    ("func6",  Fs6,  ft6),
    ("func7",  Fs7,  ft7),
    ("func8",  Fs8,  ft8),
    ("func9",  Fs9,  ft9),
    ("func10", Fs10, ft10),
]

methods = [
    ("Stehfest", nilt.Stehfest()),
    ("Talbot",   nilt.Talbot()),
    ("DeHoog",   nilt.DeHoog()),
]

t_values = np.arange(1.0, 11.0)

out_dir = os.path.join(os.path.dirname(__file__), "build")
os.makedirs(out_dir, exist_ok=True)

for fname, Fs, ft in functions:
    for mname, algo in methods:
        rows = []
        for t in t_values:
            fta = ft(t)
            if isinstance(fta, complex):
                fta = fta.real
            ftn = nilt.invert(algo, Fs, t)
            err = abs(ftn - fta) / max(abs(fta), 1e-30)
            rows.append([t, fta, ftn, err])

        out = os.path.join(out_dir, f"py_{fname}_{mname}.csv")
        np.savetxt(out, np.array(rows),
                   delimiter=",", header="t,fta,ftn,err", comments="")
        print(out)
