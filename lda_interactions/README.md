This directory stores the interaction potentials used in the LDA analysis.

- `well_width_and_density.csv` lists the quantum well widths (Å) and planar electron densities (cm⁻²) for each parameter set. Row `i` corresponds to the potential stored in `i.csv`.
- `i.csv` tabulates the interaction as two columns:
  - `r`: particle separation in units of the magnetic length `ℓ_B`.
  - `V(r)`: interaction energy in units of `e² / (ε ℓ_B)`.
