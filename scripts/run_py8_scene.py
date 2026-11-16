"""Py8 Q5: per-date scene statistics over the baseline coral mask (medians,
fixed units). Prints aggregates only; no pixel values are written.

Columns: brightness (sum of the 8 bands), Blue, Green1, Red, raw NIR, the
glint correction (raw minus glint-corrected) for Blue/Green1/Red, the share
of raw Blue removed by the glint correction, the share of NIR pixels that
compat mode floors at 1e-6, and 2*Blue/(Red+Green1).
"""
import warnings

import numpy as np

from coralremote import io
from coralremote.features import prepare_spectra
from coralremote.preprocess import raw_stack
from coralremote.units import COMPAT, FIXED

warnings.simplefilter("ignore")
B, G1, R, NIR = 1, 2, 5, 7
g = io.load_glint_stack()
r = raw_stack()
m = io.load_baseline_mask()[0]
print("date      bright  B      G1     R      NIRraw | glint B  G1     R      | glint%B | NIRfloor% | 2B/(R+G)")
for d in g.dates:
    G = prepare_spectra(g.get(d), FIXED)[m]
    Rw = r.get(d)[m]
    ok = np.isfinite(G).all(1) & np.isfinite(Rw).all(1)
    G, Rw = G[ok], Rw[ok]
    med = np.median(G, 0)
    gl = np.median(Rw - G, 0)
    floor = 100 * (prepare_spectra(g.get(d)[m][ok], COMPAT)[:, NIR] < COMPAT.floor).mean()
    print(f"{d}  {np.median(G.sum(1)):.4f}  " + "  ".join(f"{med[i]:.4f}" for i in (B, G1, R))
          + f"  {np.median(Rw[:, NIR]):.4f} | " + "  ".join(f"{gl[i]:.4f}" for i in (B, G1, R))
          + f" | {100 * np.median((Rw[:, B] - G[:, B]) / Rw[:, B]):6.1f}  | {floor:7.1f}   | "
          + f"{np.nanmedian(np.where(G[:, R] + G[:, G1] > 0, 2 * G[:, B] / (G[:, R] + G[:, G1]), np.nan)):.3f}")
