"""Remove only a/sqrt(N) from the two-QH gap E(L_rel=1) - E(L_rel=max).

Finds the a minimising the (error-weighted) spread of dE(N) - a/sqrt(N), i.e.
a weighted least-squares fit of dE(N) = c + a/sqrt(N), for each dataset in
qh_gap_vs_N.DATASETS.  With two sizes (nu=2/5) the fit is exact.

Usage: python3 analysis/qh_gap_sqrtN_fit.py
"""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qh_gap_vs_N import BLUE, DATASETS, ORANGE, OUT, YLABEL, gaps


def fit(N, y, s):
    """Weighted LSQ of y = c + a/sqrt(N). Returns (c, a), covariance, chi2."""
    X = np.column_stack([np.ones_like(N, dtype=float), 1 / np.sqrt(N)])
    W = 1 / s**2
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    p = cov @ (X.T @ (W * y))
    return p, cov, np.sum(W * (y - X @ p) ** 2)


def main():
    for tag, (directory, title) in DATASETS.items():
        N, _, gap, gap_err = gaps(directory)
        (c, a), cov, chi2 = fit(N, gap, gap_err)
        dc, da = np.sqrt(np.diag(cov))
        corrected = gap - a / np.sqrt(N)
        print(f"== {title}")
        print(f"  a = {a:+.6f} +/- {da:.6f}")
        print(f"  c = dE(N->inf) = {c:+.6f} +/- {dc:.6f}")
        print(f"  chi2 = {chi2:.3f} for {len(N) - 2} dof")
        print(f"  spread (std) raw = {gap.std():.6f}, corrected = {corrected.std():.6f}")
        for n, y, e in zip(N, corrected, gap_err):
            print(f"    N={n:3d}  dE - a/sqrt(N) = {y:+.6f} +/- {e:.6f}")

        with open(os.path.join(OUT, f"qh_{tag}_gap_sqrtN_fit.txt"), "w") as f:
            f.write("Fit dE(N) = c + a/sqrt(N), weighted by MCMC errors\n")
            f.write(f"a = {a:+.6e} +/- {da:.6e}\n")
            f.write(f"c = {c:+.6e} +/- {dc:.6e}\n")
            f.write(f"chi2 = {chi2:.4f}, dof = {len(N) - 2}\n")
            f.write("N,dE - a/sqrt(N),error\n")
            for n, y, e in zip(N, corrected, gap_err):
                f.write(f"{n},{y:.6e},{e:.6e}\n")

        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
        ax.errorbar(N, gap, yerr=gap_err, fmt="o", ms=7, capsize=4, lw=1.5, color=BLUE)
        Nf = np.linspace(N.min(), N.max(), 200)
        ax.plot(Nf, c + a / np.sqrt(Nf), color=BLUE, lw=1, alpha=0.5)
        ax.set_ylabel(YLABEL)
        ax.set_title(title + " (raw)")

        ax2.errorbar(N, corrected, yerr=gap_err, fmt="s", ms=6, capsize=4, lw=1.5,
                     color=ORANGE)
        ax2.axhline(c, color="0.4", lw=1, ls="--")
        ax2.set_ylabel(r"$\Delta E(N) - a/\sqrt{N}$  $[e^2/\epsilon\ell]$")
        ax2.set_title(rf"corrected: $a={a:+.4f}\pm{da:.4f}$", fontsize=10.5)
        ax2.text(0.03, 0.95, rf"$N\to\infty$: ${c:+.5f}\pm{dc:.5f}$", va="top",
                 transform=ax2.transAxes, fontsize=8.5, color="0.3")

        for axi in (ax, ax2):
            axi.set_xlabel(r"$N$")
            axi.grid(alpha=0.25)
            axi.set_xticks(N)
            axi.set_xlim(N.min() - 5, N.max() + 5)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"qh_{tag}_gap_sqrtN_fit.png"), dpi=200)


if __name__ == "__main__":
    main()
