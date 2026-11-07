# radialFundSelect

**radialFundSelect** implements the three positive-alpha fund selection procedures in Long Feng's *Optimal radial weighting for robust mutual fund selection*. Apply the methods to your own return matrix or rerun the Monte Carlo experiments and rolling fund portfolios for the November 2026 manuscript.

## Installation

```r
install.packages("remotes")  # once
remotes::install_github("flnankai/radialFundSelect", ref = "v0.3.0")
library(radialFundSelect)
```

The package requires R 4.1 or later. The estimators, simulations and PDF plots run in R. Reading the authors' original MAT file also requires `install.packages("R.matlab")`. The aligned CSV input uses base R.

## Call each method on your own data

Supply a numeric matrix of **excess returns**, one date per row and one fund per column. Supply pricing factors in the same date order, without an intercept column. Use a complete estimation window. With no pricing factors, set `factors = NULL`. If your input contains raw monthly fund returns, first subtract the corresponding monthly risk-free return in the same units.

```r
panel <- simulate_panel(T = 60, N = 40, distribution = "t3", seed = 10)
Y <- panel$Y
F <- panel$factors

inverse <- select_in_bh(Y, F, q = 0.05)
signs   <- select_ss_bh(Y, F, q = 0.05)
ols     <- select_ls_bh(Y, F, q = 0.05)

which(inverse$selected)          # selected columns
inverse$pvalue                   # one-sided p-values in column order
inverse$fit_converged            # check before interpreting the fitted values
inverse$fit$alpha_first_order    # fund-specific alpha estimates
```

IN-BH fixes pairwise-quantile scales and uses smoothing strength 0.5. Set `smoothing = 0` for the unregularized comparison. SS-BH jointly estimates location and diagonal scale. LS-BH uses the OLS intercept with residual variance divisor T. Each method uses the requested BH level. The shared interfaces `select_funds()` and `fit_three()` are also available. See `?select_funds` for the estimators, numerical checks and failure policy.

## Reproduce the paper

| Paper output | Fresh computation | Returned results or PDF |
| --- | --- | --- |
| Section 5, Figures 1 and 2 | Signal-strength grid, T = 60 and 120 | `simulation_facts`, `signal_T60.pdf`, `signal_T120.pdf` |
| Section 5, Figure 3 | Residual-correlation grid, T = 60 | `study$summary`, `study$paired_differences`, `correlation.pdf` |
| Section 5, effect of smoothing | Pair unregularized and smoothed fits on the same simulated panels | `smoothing`, `smoothing_facts` |
| Section 6, Figure 4 | Four-factor selection on 69-month histories, followed by 3-, 6- or 9-month holding periods | `empirical$summary`, baseline and stressed wealth paths, `empirical_T69.pdf` |

### Simulations and smoothing comparison

A short run exercises every setting and figure. Its Monte Carlo accuracy depends on the number of repetitions.

```r
quick <- reproduce_simulations(reps = 2, workers = 1, output_dir = tempdir())
quick$study$summary
```

To regenerate the paper's estimates, use 2,000 repetitions at each of the 78 distinct settings and seed 20261105, embedded in the reproduction function. The figures include pointwise Monte Carlo intervals. The paired smoothing comparison uses the same draws. Workers process panels in memory. Choose a worker count that fits your machine, up to 25.

```r
sim <- reproduce_simulations(reps = 2000, workers = 4, output_dir = "paper-results")
sim$simulation_facts
sim$smoothing_facts
sim$study$paired_differences
```

The full simulation workflow took about 7 minutes 25 seconds with 25 workers and R 4.5.2 on the authors' Windows machine. Runtime depends on the CPU and the number of workers.

To rerun an individual case, select its row from the design inventory. Use `simulate_panel()` for a single draw at your own signal increment or residual correlation. The full experiment uses Gaussian, Student t with three degrees of freedom, and normal-mixture errors. For each distribution the signal designs vary nine increments at fixed correlation 0.5. The correlation designs vary nine correlations at signal increment 0.48. One setting per distribution contributes to both displays. [All settings, identifiers and commands](inst/reproduce/simulation-scenarios.md)

```r
cases <- scenario_inventory()
print(cases)
example_case <- subset(cases, distribution == "t3" & T == 60 & abs(delta - 0.48) < 1e-9 & abs(rho - 0.5) < 1e-9)
one <- simulate_study(scenarios = example_case, reps = 20, workers = 1)
one$summary
```

### Empirical results and Figure 4

Obtain `CaseStudy/mutual_fund.mat` from the Gao, Peng and Yan replication archive and the historical French monthly RF ZIP. The package reads the MAT file directly and aligns the risk-free series by calendar month. The inputs are distributed by their respective providers. [Data sources, download links and CSV alternative](inst/replication/empirical-data.md)

```r
install.packages("R.matlab")  # once, for MAT input
fund_data <- "path/to/CaseStudy/mutual_fund.mat"
risk_free <- "path/to/F-F_Research_Data_Factors_CSV.zip"
emp <- reproduce_empirical(fund_data, risk_free, figure = "empirical_T69.pdf")
emp$summary
```

The model uses fund excess returns and the archived market, size, value and momentum factors. Each window uses its own historically complete fund universe, BH at 5%, equal initial allocations and buy-and-hold weights. Raw returns compound during holding periods. Empty selections earn the risk-free rate. The summary includes compound growth, maximum drawdown, a common April 2017 endpoint and a full-loss sensitivity check for missing held positions. Reported rates are decimals.

### Run both sections together

```r
paper <- reproduce_paper(data = fund_data, rf_file = risk_free,
                         output_dir = "paper-results", workers = 4)
paper$figures
paper$empirical$summary
paper$simulations$smoothing_facts$overall_paired
```

The four PDFs are written to `output_dir`. The returned R object contains the numerical summaries. Use `citation("radialFundSelect")` for the working-paper citation.

## License

The R package is distributed under GPL-3. Third-party data remain subject to their providers' terms. [License text](LICENSE.md)
