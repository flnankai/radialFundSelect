# radialFundSelect 0.3.0

- Adds named IN-BH, SS-BH and LS-BH selection functions.
- Recomputes the paper's simulation contrasts, three rolling portfolios and four figures.
- Adds original MAT and French CSV/ZIP input alongside the aligned CSV route.

# radialFundSelect 0.2.0

- The method name remains IN-BH (`method = "IN"`); smoothing is an
  implementation parameter, not a new method name.
- IN defaults to relative smoothing strength 0.5 in `select_funds`,
  `fit_three`, and `simulate_study`. Set `smoothing = 0` for the legacy fit.
- Lambda is fixed from PQR-standardized unrestricted OLS residual norms.
  The positive-lambda fit uses the general radial variance expression,
  permits individual zero residual radii, and retains strict score checks.
- Simulation metadata records the smoothing strength, implementation,
  pilot rule, and studentizer. Intermediate simulation data remain in memory.
- The original Python comparison scripts explicitly use the legacy control.
- Regularization is distinguished from exact inverse-norm optimality and
  from an unproved FDR guarantee for the pilot-dependent rule.
