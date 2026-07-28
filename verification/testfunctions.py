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

import numpy as np


pi = np.pi


# Stehfest (1970) test functions

# f1: f(t) = 1/sqrt(pi*t),        F(s) = 1/sqrt(s)
def Fs1(s): return 1.0 / np.sqrt(s)
def ft1(t): return 1.0 / np.sqrt(pi * t)

# f2: f(t) = -gamma - ln(t),       F(s) = ln(s)/s
def Fs2(s): return np.log(s) / s
def ft2(t): return -0.57722 - np.log(t)

# f3: f(t) = t^3/6,                F(s) = 1/s^4
def Fs3(s): return s ** (-4.0)
def ft3(t): return t ** 3 / 6.0

# f4: f(t) = exp(-t),              F(s) = 1/(s+1)
def Fs4(s): return 1.0 / (s + 1.0)
def ft4(t): return np.exp(-t)

# f5: f(t) = sin(sqrt(2t)),        F(s) = sqrt(pi/(2s^3)) * exp(-1/(2s))
def Fs5(s): return np.sqrt(pi / 2.0) * (s ** -1.5) * np.exp(-1.0 / (2.0 * s))
def ft5(t): return np.sin(np.sqrt(2.0 * t))

# Abate & Whitt (2006) test functions

# f6: f(t) = t,                    F(s) = 1/s^2
def Fs6(s): return 1.0 / (s * s)
def ft6(t): return t

# f7: f(t) = t*exp(-t),            F(s) = 1/(s+1)^2
def Fs7(s): return 1.0 / ((s + 1.0) ** 2)
def ft7(t): return t * np.exp(-t)

# f8: f(t) = sin(t),               F(s) = 1/(s^2+1)
def Fs8(s): return 1.0 / (s * s + 1.0)
def ft8(t): return np.sin(t)

# f9: f(t) = cos(t),               F(s) = s/(s^2+1)
def Fs9(s): return s / (s * s + 1.0)
def ft9(t): return np.cos(t)

# f10: f(t) = exp(-t)*sin(t),      F(s) = 1/((s+1)^2+1)
def Fs10(s): return 1.0 / ((s + 1.0) ** 2 + 1.0)
def ft10(t): return np.exp(-t) * np.sin(t)


T_VALUES = list(range(1, 11))

TEST_FUNCTIONS = {
    "f1" : {"ft": ft1,  "Fs": Fs1 },
    "f2" : {"ft": ft2,  "Fs": Fs2 },
    "f3" : {"ft": ft3,  "Fs": Fs3 },
    "f4" : {"ft": ft4,  "Fs": Fs4 },
    "f5" : {"ft": ft5,  "Fs": Fs5 },
    "f6" : {"ft": ft6,  "Fs": Fs6 },
    "f7" : {"ft": ft7,  "Fs": Fs7 },
    "f8" : {"ft": ft8,  "Fs": Fs8 },
    "f9" : {"ft": ft9,  "Fs": Fs9 },
    "f10": {"ft": ft10, "Fs": Fs10}
}