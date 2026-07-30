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
import nilt

from testfunctions import TEST_FUNCTIONS


methods = [
    ("Stehfest", nilt.Stehfest()),
    ("Talbot",   nilt.Talbot()),
    ("DeHoog",   nilt.DeHoog()),
]

t_values = np.arange(1.0, 11.0)

for fname, funcs in TEST_FUNCTIONS.items():
    Fs = funcs["Fs"]
    ft = funcs["ft"]
    for mname, algo in methods:
        rows = []
        for t in t_values:
            fta = ft(t)
            if isinstance(fta, complex):
                fta = fta.real
            ftn = nilt.invert(algo, Fs, t)
            err = abs(ftn - fta) / max(abs(fta), 1e-30)
            rows.append([t, fta, ftn, err])

        out = f"py_{fname}_{mname}.csv"
        np.savetxt(out, np.array(rows),
                   delimiter=",", header="t,fta,ftn,err", comments="")
        print(out)
