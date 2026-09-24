"""E(L_rel=1) - E(L_rel=max) for two quasiholes at nu=2/5, versus system size N.

Reads 2_qhs_2_5_filling_factor_coulomb_dispersion/N_particles.csv
(columns: L_rel as 'p//q', energy, MCMC error) and plots the gap vs N.
"""
import glob
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "2_qhs_2_5_filling_factor_coulomb_dispersion")
OUT = os.path.dirname(os.path.abspath(__file__))


def load(path):
    L, E, err = [], [], []
    for line in open(path):
        if not line.strip():
            continue
        l, e, s = [x.strip() for x in line.split(",")]
        p, q = l.split("//")
        L.append(int(p) / int(q))
        E.append(float(e))
        err.append(float(s))
    return np.array(L), np.array(E), np.array(err)


rows = []
for path in glob.glob(os.path.join(DATA, "*_particles.csv")):
    N = int(re.match(r"(\d+)_particles", os.path.basename(path)).group(1))
    L, E, err = load(path)
    i1, imax = np.argmin(np.abs(L - 1)), np.argmax(L)
    gap = E[i1] - E[imax]
    gap_err = np.hypot(err[i1], err[imax])
    rows.append((N, L[imax], gap, gap_err))
rows.sort()
N, Lmax, gap, gap_err = map(np.array, zip(*rows))

with open(os.path.join(OUT, "qh_2_5_gap_vs_N.csv"), "w") as f:
    f.write("N,L_max,E(L=1)-E(L_max),error\n")
    for r in rows:
        f.write(f"{r[0]},{r[1]:g},{r[2]:.6e},{r[3]:.6e}\n")
for r in rows:
    print(f"N={r[0]:3d}  L_max={r[1]:g}  dE={r[2]:+.6f} +/- {r[3]:.6f}")

fig, ax = plt.subplots(figsize=(5, 3.6))
ax.errorbar(N, gap, yerr=gap_err, fmt="o", ms=7, capsize=4, lw=1.5, color="#2a6fdb")
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$E(L_{\rm rel}=1)-E(L_{\rm rel}^{\max})$  $[e^2/\epsilon\ell]$")
ax.set_title(r"2 QHs at $\nu=2/5$")
ax.grid(alpha=0.25)
ax.set_xlim(N.min() - 5, N.max() + 5)
ax.set_xticks(N)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "qh_2_5_gap_vs_N.png"), dpi=200)
