# Fig. 1c Density Profiles

This directory contains the radial density profiles used in Fig. 1c. Each file encodes the density of either a quasiparticle (QP) or quasihole (QH) molecule at a specified filling fraction.

**File naming**
- `qp_q_r_s.csv`: QP\_{q} molecule evaluated at filling `r/s`.
- `qh_q_r_s.csv`: QH\_{q} molecule evaluated at filling `r/s`.

**Data format**
- Column 1 (`r`): Radial distance from the molecule center.
- Column 2 (`delta_nu_over_nu`): Relative density change `Δν(r)/ν` with respect to the background FQH state.

Values are stored as comma-separated columns with no header.
