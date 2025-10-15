# qh2_qh1_3_7_filling_factor

This folder contains Monte Carlo averages for a QH\_2 molecule anchored at the north pole and a QH\_1 molecule at the south pole at filling 3/7.

**File naming**
- `N_particles_density.csv`: `N` is the total number of particles in the configuration.
- `N_particles_pair_density.csv`: same `N` as above.

**Columns**
- `N_particles_density.csv`: distance from the north pole, density normalised to the bulk filling.
- `N_particles_pair_density.csv`: particle separation, pair-distribution probability.

These profiles can be combined with arbitrary interaction potentials to evaluate molecular anyon binding energies. To reduce storage, only the mean profiles across all MCMC chains are stored.
