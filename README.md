# molecular-anyons-gattu-jain

Numerical datasets for the molecular-anyon study, organised by observable and filling fraction. Each subdirectory ships with a README covering file formats, units, and analysis tips.

## Directory Guide

**Dispersion spectra**
- `2_qhs_5_2_dispersion/`: Mean pair-density dispersions for topological and non-topological quasiholes at ν = 5/2.
- `*_qhs_*_filling_factor_coulomb_dispersion/`: Coulomb-dispersion curves for quasihole molecules at fillings such as 2/5 and 3/7 (`N_particles.csv` with relative angular momentum, energy, and MCMC error).
- `*_qps_*_filling_factor_coulomb_dispersion/`: Coulomb-dispersion curves for quasiparticle molecules (same layout as the quasihole directories).

**Molecular density profiles**
- `qpA_qpB_r_s_filling_factor/`: Average radial and pair-density profiles for QP\_A at the north pole and QP\_B at the south pole at filling `r/s`.
- `qhA_qhB_r_s_filling_factor/`: Analogous density profiles for quasihole molecules.
- `densities/`: Fig. 1c radial-density datasets; filenames encode QP/QH count and filling.

**Interaction inputs**
- `lda_interactions/`: Interaction tables indexed by `well_width_and_density.csv`, giving `r` and `V(r)` pairs for each quantum-well width and electron density.

**Other**
- `LICENSE`: Usage terms for the released data.

See each directory’s README before consuming the data—those notes document column definitions, sampling details, and suggestions for recomputing energies with alternative interactions.

## Citing

If you use this work, please cite the paper:

```bibtex
@article{scl5-8pv6,
  title = {Molecular Anyons in the Fractional Quantum Hall Effect},
  author = {Gattu, Mytraya and Jain, J. K.},
  journal = {Phys. Rev. Lett.},
  volume = {135},
  issue = {23},
  pages = {236601},
  numpages = {7},
  year = {2025},
  month = {Dec},
  publisher = {American Physical Society},
  doi = {10.1103/scl5-8pv6},
  url = {https://link.aps.org/doi/10.1103/scl5-8pv6}
}
```
