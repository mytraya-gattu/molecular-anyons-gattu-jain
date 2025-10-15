qhs_5_2_dispersion

The files in this directory store mean pair density functions for quasihole states at filling factor 5/2. They are generated separately for topological and non-topological sectors and follow these naming patterns:

```
topological_qhs_dispersion_<N>_particles_<Lnum>_<Ldenom>.csv
non_topological_qhs_dispersion_<N>_particles_<Lnum>_<Ldenom>.csv
```

Here:
- `N` is the total number of particles in the system.
- `Lnum` is the numerator of the relative angular momentum.
- `Ldenom` is the denominator of the relative angular momentum.
- `L=Lnum/Ldenom` is the relative angular momentum.

Each CSV file provides the mean pair density as a function of interparticle separation for the specified `L`. By convolving this density with an arbitrary interaction potential, the energy of the corresponding state can be obtained.
