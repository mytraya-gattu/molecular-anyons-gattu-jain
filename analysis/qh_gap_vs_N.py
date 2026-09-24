"""E(L_rel=1) - E(L_rel=max) for two quasiholes versus system size N.

For each dataset in DATASETS, reads <dir>/N_particles.csv (columns: L_rel as
'p//q', energy, MCMC error), computes the gap dE(N) = E(L=1) - E(L_max) and
plots it vs N.  When there are >= 3 sizes, it also finds a, b minimising the
(error-weighted) spread of dE(N) - a/N - b/sqrt(N).  That is a weighted
least-squares fit of dE(N) = c + a/N + b/sqrt(N), where c is the N -> inf gap.

Usage: python3 analysis/qh_gap_vs_N.py
"""
import glob
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.dirname(os.path.abspath(__file__))
DATASETS = {
    "2_5": ("2_qhs_2_5_filling_factor_coulomb_dispersion", r"2 QHs at $\nu=2/5$"),
    "3_7": ("2_qhs_3_7_filling_factor_coulomb_dispersion", r"2 QHs at $\nu=3/7$"),
}
BLUE, ORANGE = "#2a6fdb", "#e0701a"
YLABEL = r"$E(L_{\rm rel}=1)-E(L_{\rm rel}^{\max})$  $[e^2/\epsilon\ell]$"


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


def gaps(directory):
    rows = []
    for path in glob.glob(os.path.join(ROOT, directory, "*_particles.csv")):
        N = int(re.match(r"(\d+)_particles", os.path.basename(path)).group(1))
        L, E, err = load(path)
        i1, imax = np.argmin(np.abs(L - 1)), np.argmax(L)
        rows.append((N, L[imax], E[i1] - E[imax], np.hypot(err[i1], err[imax])))
    rows.sort()
    return [np.array(c) for c in zip(*rows)]


def fit(N, y, s):
    """Weighted LSQ of y = c + a/N + b/sqrt(N). Returns (c, a, b), covariance, chi2."""
    X = np.column_stack([np.ones_like(N, dtype=float), 1 / N, 1 / np.sqrt(N)])
    W = 1 / s**2
    A = X.T @ (X * W[:, None])
    cov = np.linalg.inv(A)
    p = cov @ (X.T @ (W * y))
    chi2 = np.sum(W * (y - X @ p) ** 2)
    return p, cov, chi2


for tag, (directory, title) in DATASETS.items():
    N, Lmax, gap, gap_err = gaps(directory)
    print(f"== {title}")
    for row in zip(N, Lmax, gap, gap_err):
        print("  N=%3d  L_max=%2g  dE=%+.6f +/- %.6f" % row)

    with open(os.path.join(OUT, f"qh_{tag}_gap_vs_N.csv"), "w") as f:
        f.write("N,L_max,E(L=1)-E(L_max),error\n")
        for row in zip(N, Lmax, gap, gap_err):
            f.write("%d,%g,%.6e,%.6e\n" % row)

    fit_it = len(N) >= 3
    fig, axes = plt.subplots(1, 2 if fit_it else 1, figsize=(9.6 if fit_it else 5, 3.8),
                             squeeze=False)
    ax = axes[0, 0]
    ax.errorbar(N, gap, yerr=gap_err, fmt="o", ms=7, capsize=4, lw=1.5, color=BLUE)
    ax.set_ylabel(YLABEL)
    ax.set_title(title + " (raw)" if fit_it else title)

    if fit_it:
        (c, a, b), cov, chi2 = fit(N, gap, gap_err)
        da, db, dc = np.sqrt(cov[1, 1]), np.sqrt(cov[2, 2]), np.sqrt(cov[0, 0])
        rho = cov[1, 2] / (da * db)
        corrected = gap - a / N - b / np.sqrt(N)
        print(f"  a = {a:+.5f} +/- {da:.5f}")
        print(f"  b = {b:+.5f} +/- {db:.5f}   (corr(a,b) = {rho:+.3f})")
        print(f"  c = dE(N->inf) = {c:+.6f} +/- {dc:.6f}")
        print(f"  chi2 = {chi2:.3f} for {len(N) - 3} dof")
        print(f"  spread (std) raw = {gap.std():.6f}, corrected = {corrected.std():.6f}")
        for n, y in zip(N, corrected):
            print(f"    N={n:3d}  dE - a/N - b/sqrt(N) = {y:+.6f}")

        with open(os.path.join(OUT, f"qh_{tag}_gap_fit.txt"), "w") as f:
            f.write("Fit dE(N) = c + a/N + b/sqrt(N), weighted by MCMC errors\n")
            f.write(f"a = {a:+.6e} +/- {da:.6e}\n")
            f.write(f"b = {b:+.6e} +/- {db:.6e}\n")
            f.write(f"corr(a,b) = {rho:+.4f}\n")
            f.write(f"c = {c:+.6e} +/- {dc:.6e}\n")
            f.write(f"chi2 = {chi2:.4f}, dof = {len(N) - 3}\n")
            f.write("N,dE - a/N - b/sqrt(N)\n")
            for n, y in zip(N, corrected):
                f.write(f"{n},{y:.6e}\n")

        Nf = np.linspace(N.min(), N.max(), 200)
        ax.plot(Nf, c + a / Nf + b / np.sqrt(Nf), color=BLUE, lw=1, alpha=0.5)

        ax2 = axes[0, 1]
        ax2.errorbar(N, corrected, yerr=gap_err, fmt="s", ms=6, capsize=4, lw=1.5,
                     color=ORANGE)
        ax2.axhline(c, color="0.4", lw=1, ls="--")
        ax2.set_ylabel(r"$\Delta E(N) - a/N - b/\sqrt{N}$  $[e^2/\epsilon\ell]$")
        ax2.set_title(rf"corrected: $a={a:+.3f}$, $b={b:+.3f}$", fontsize=10.5)
        ax2.text(0.03, 0.04, rf"$N\to\infty$: ${c:+.4f}\pm{dc:.4f}$",
                 transform=ax2.transAxes, fontsize=8.5, color="0.3")

    for axi in axes.ravel():
        axi.set_xlabel(r"$N$")
        axi.grid(alpha=0.25)
        axi.set_xticks(N)
        axi.set_xlim(N.min() - 5, N.max() + 5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"qh_{tag}_gap_vs_N.png"), dpi=200)
