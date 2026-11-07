# Monte Carlo scenarios

Section 5 uses 200 funds and three observable factors. The funds share the same factor observations in a panel. Twenty funds have positive alphas, twenty have negative alphas, and the remaining 160 have zero alphas. The signal increment shifts the positive alphas. The residual correlation between funds separated by k positions is rho to the power k.

The table lists the 78 distinct settings returned by `scenario_inventory()` and identifies the figures in which each setting appears. Figure 1 varies signal strength at T = 60. Figure 2 repeats the signal experiment at T = 120. Figure 3 varies residual correlation at T = 60. Each distribution contributes one setting common to Figures 1 and 3.

| Scenario ID | Error distribution | Months | Residual correlation | Positive-alpha increment | Figure |
|---|---|---:|---:|---:|---|
| `normal_T60_rho0.50_delta0.06` | Gaussian | 60 | 0.50 | 0.06 | 1 |
| `normal_T60_rho0.50_delta0.12` | Gaussian | 60 | 0.50 | 0.12 | 1 |
| `normal_T60_rho0.50_delta0.18` | Gaussian | 60 | 0.50 | 0.18 | 1 |
| `normal_T60_rho0.50_delta0.24` | Gaussian | 60 | 0.50 | 0.24 | 1 |
| `normal_T60_rho0.50_delta0.30` | Gaussian | 60 | 0.50 | 0.30 | 1 |
| `normal_T60_rho0.50_delta0.36` | Gaussian | 60 | 0.50 | 0.36 | 1 |
| `normal_T60_rho0.50_delta0.42` | Gaussian | 60 | 0.50 | 0.42 | 1 |
| `normal_T60_rho0.50_delta0.48` | Gaussian | 60 | 0.50 | 0.48 | 1 and 3 |
| `normal_T60_rho0.50_delta0.54` | Gaussian | 60 | 0.50 | 0.54 | 1 |
| `normal_T120_rho0.50_delta0.04` | Gaussian | 120 | 0.50 | 0.04 | 2 |
| `normal_T120_rho0.50_delta0.08` | Gaussian | 120 | 0.50 | 0.08 | 2 |
| `normal_T120_rho0.50_delta0.12` | Gaussian | 120 | 0.50 | 0.12 | 2 |
| `normal_T120_rho0.50_delta0.16` | Gaussian | 120 | 0.50 | 0.16 | 2 |
| `normal_T120_rho0.50_delta0.20` | Gaussian | 120 | 0.50 | 0.20 | 2 |
| `normal_T120_rho0.50_delta0.24` | Gaussian | 120 | 0.50 | 0.24 | 2 |
| `normal_T120_rho0.50_delta0.28` | Gaussian | 120 | 0.50 | 0.28 | 2 |
| `normal_T120_rho0.50_delta0.32` | Gaussian | 120 | 0.50 | 0.32 | 2 |
| `normal_T120_rho0.50_delta0.36` | Gaussian | 120 | 0.50 | 0.36 | 2 |
| `normal_T60_rho0.10_delta0.48` | Gaussian | 60 | 0.10 | 0.48 | 3 |
| `normal_T60_rho0.20_delta0.48` | Gaussian | 60 | 0.20 | 0.48 | 3 |
| `normal_T60_rho0.30_delta0.48` | Gaussian | 60 | 0.30 | 0.48 | 3 |
| `normal_T60_rho0.40_delta0.48` | Gaussian | 60 | 0.40 | 0.48 | 3 |
| `normal_T60_rho0.60_delta0.48` | Gaussian | 60 | 0.60 | 0.48 | 3 |
| `normal_T60_rho0.70_delta0.48` | Gaussian | 60 | 0.70 | 0.48 | 3 |
| `normal_T60_rho0.80_delta0.48` | Gaussian | 60 | 0.80 | 0.48 | 3 |
| `normal_T60_rho0.90_delta0.48` | Gaussian | 60 | 0.90 | 0.48 | 3 |
| `t3_T60_rho0.50_delta0.06` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.06 | 1 |
| `t3_T60_rho0.50_delta0.12` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.12 | 1 |
| `t3_T60_rho0.50_delta0.18` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.18 | 1 |
| `t3_T60_rho0.50_delta0.24` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.24 | 1 |
| `t3_T60_rho0.50_delta0.30` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.30 | 1 |
| `t3_T60_rho0.50_delta0.36` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.36 | 1 |
| `t3_T60_rho0.50_delta0.42` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.42 | 1 |
| `t3_T60_rho0.50_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.48 | 1 and 3 |
| `t3_T60_rho0.50_delta0.54` | Covariance-standardized Student t with 3 df | 60 | 0.50 | 0.54 | 1 |
| `t3_T120_rho0.50_delta0.04` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.04 | 2 |
| `t3_T120_rho0.50_delta0.08` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.08 | 2 |
| `t3_T120_rho0.50_delta0.12` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.12 | 2 |
| `t3_T120_rho0.50_delta0.16` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.16 | 2 |
| `t3_T120_rho0.50_delta0.20` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.20 | 2 |
| `t3_T120_rho0.50_delta0.24` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.24 | 2 |
| `t3_T120_rho0.50_delta0.28` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.28 | 2 |
| `t3_T120_rho0.50_delta0.32` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.32 | 2 |
| `t3_T120_rho0.50_delta0.36` | Covariance-standardized Student t with 3 df | 120 | 0.50 | 0.36 | 2 |
| `t3_T60_rho0.10_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.10 | 0.48 | 3 |
| `t3_T60_rho0.20_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.20 | 0.48 | 3 |
| `t3_T60_rho0.30_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.30 | 0.48 | 3 |
| `t3_T60_rho0.40_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.40 | 0.48 | 3 |
| `t3_T60_rho0.60_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.60 | 0.48 | 3 |
| `t3_T60_rho0.70_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.70 | 0.48 | 3 |
| `t3_T60_rho0.80_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.80 | 0.48 | 3 |
| `t3_T60_rho0.90_delta0.48` | Covariance-standardized Student t with 3 df | 60 | 0.90 | 0.48 | 3 |
| `mixture_T60_rho0.50_delta0.06` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.06 | 1 |
| `mixture_T60_rho0.50_delta0.12` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.12 | 1 |
| `mixture_T60_rho0.50_delta0.18` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.18 | 1 |
| `mixture_T60_rho0.50_delta0.24` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.24 | 1 |
| `mixture_T60_rho0.50_delta0.30` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.30 | 1 |
| `mixture_T60_rho0.50_delta0.36` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.36 | 1 |
| `mixture_T60_rho0.50_delta0.42` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.42 | 1 |
| `mixture_T60_rho0.50_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.48 | 1 and 3 |
| `mixture_T60_rho0.50_delta0.54` | Symmetric two-scale normal mixture | 60 | 0.50 | 0.54 | 1 |
| `mixture_T120_rho0.50_delta0.04` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.04 | 2 |
| `mixture_T120_rho0.50_delta0.08` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.08 | 2 |
| `mixture_T120_rho0.50_delta0.12` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.12 | 2 |
| `mixture_T120_rho0.50_delta0.16` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.16 | 2 |
| `mixture_T120_rho0.50_delta0.20` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.20 | 2 |
| `mixture_T120_rho0.50_delta0.24` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.24 | 2 |
| `mixture_T120_rho0.50_delta0.28` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.28 | 2 |
| `mixture_T120_rho0.50_delta0.32` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.32 | 2 |
| `mixture_T120_rho0.50_delta0.36` | Symmetric two-scale normal mixture | 120 | 0.50 | 0.36 | 2 |
| `mixture_T60_rho0.10_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.10 | 0.48 | 3 |
| `mixture_T60_rho0.20_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.20 | 0.48 | 3 |
| `mixture_T60_rho0.30_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.30 | 0.48 | 3 |
| `mixture_T60_rho0.40_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.40 | 0.48 | 3 |
| `mixture_T60_rho0.60_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.60 | 0.48 | 3 |
| `mixture_T60_rho0.70_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.70 | 0.48 | 3 |
| `mixture_T60_rho0.80_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.80 | 0.48 | 3 |
| `mixture_T60_rho0.90_delta0.48` | Symmetric two-scale normal mixture | 60 | 0.90 | 0.48 | 3 |

Both the main three-method experiment and the smoothing comparison use the same settings and L'Ecuyer-CMRG random streams. The paper uses 2,000 repetitions, seed 20261105, and BH at 0.10. For each comparison, a failed fit makes zero rejections and remains in the denominator. `simulate_panel()` generates a chosen setting and repetition. `reproduce_simulations()` recomputes the setting-level and paired results used in Section 5.
