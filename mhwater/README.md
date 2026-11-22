# mhwater

STAT\*4000 (Fall 2026) project package: Metropolis-Hastings sampling, with an
application to Ontario PWQMN 2025 total phosphorus data.

Authors: Félix Mainville, Matthew Kai Kai

## Install

From the package folder (open `mhwater.Rproj` in RStudio):

```r
devtools::document()   # rebuild help files + NAMESPACE after editing R/
devtools::install()    # install the package
devtools::test()       # run the tests
devtools::check()      # full R CMD check
```

## What is in it

| Function | Purpose |
|---|---|
| `read_pwqmn()` | Read the raw PWQMN CSV (Latin-1, all columns as text) |
| `clean_tp()` | Keep one analyte, drop duplicate rows, handle `<DL` values |
| `describe_tp()` | Descriptive summary + empirical exceedance proportion |
| `fit_candidates()` | MLE fits of lognormal, gamma, Weibull (or others) |
| `compare_fits()` | AIC/BIC, KS/CvM/AD, exceedance probabilities, quantiles |
| `exceed_prob()` | P(X > threshold) under a fitted distribution |
| `make_target()` | Turn a fit into an MH target (density, log-density, true P(X > 30)) |
| `mh_sampler()` | Metropolis-Hastings sampler with a pluggable proposal; returns chain, acceptance rate, proposal name |
| `normal_rw()`, `lognormal_rw()`, `gamma_independence()` | Proposals for `mh_sampler()`: normal random walk, multiplicative (log-scale) random walk, gamma independence sampler |
| `plot_fit_diagnostics()` | 4-panel plot: exceedance curve, Q-Q, CDF, P-P |
| `plot_log_density()` | Fitted densities on the log10 scale |
| `censoring_check()`, `site_median_check()` | Sensitivity checks |
| `tp_clean` | Cleaned dataset (3,690 samples, 449 sites) |

Running the sampler on the lognormal target (`tgt` from the quick example below):

```r
chain <- mh_sampler(tgt$log_density, 10000, x0 = 10, proposal = normal_rw(20))
mean(chain$x[-(1:1000)] > 30)
```

## Reproducing the environmental analysis

```r
source(system.file("analysis", "environmental_fit.R", package = "mhwater"))
```

The cleaned dataset is rebuilt from the raw CSV with `data-raw/tp_clean.R`.

## Quick example

```r
library(mhwater)
fits <- fit_candidates(tp_clean$value)
compare_fits(fits, threshold = 30)
tgt <- make_target(fits$lnorm)
tgt$log_density(10)   # pass this to the MH sampler
tgt$p_exceed          # exact answer the MH estimate should reproduce
```
