# =============================================================================
# Environmental application (Objective 2): reproduce the distribution
# comparison for PWQMN 2025 total phosphorus using the mhwater package.
#
# Run from any folder after installing the package (devtools::install()).
# Writes tables to ./output and figures to ./figures.
# =============================================================================
library(mhwater)

set.seed(4000)
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

PWQO_TP <- 30   # Ontario PWQO for total phosphorus, rivers & streams (ug/L)

x <- tp_clean$value

# 1. Descriptive summary -------------------------------------------------------
desc <- describe_tp(x, threshold = PWQO_TP)
print(desc, width = Inf)

# 2. Fit and compare candidate distributions -----------------------------------
fits <- fit_candidates(x, families = c("lnorm", "gamma", "weibull"))
comparison <- compare_fits(fits, threshold = PWQO_TP)
print(comparison, width = Inf)
utils::write.csv(comparison, "output/fit_comparison.csv", row.names = FALSE)

# 3. Figures --------------------------------------------------------------------
plot_fit_diagnostics(fits, threshold = PWQO_TP,
                     file = "figures/tp_fit_diagnostics.png")
ggplot2::ggsave("figures/tp_fitted_densities_log.png",
                plot_log_density(fits, threshold = PWQO_TP),
                width = 7, height = 4.5, dpi = 300)

# 4. Sensitivity checks ---------------------------------------------------------
print(censoring_check(tp_clean, "lnorm"))   # DL/2 vs censored MLE
print(site_median_check(tp_clean))           # one value per site

# 5. Target for the Metropolis-Hastings sampler ---------------------------------
best <- comparison$distribution[1]
tgt  <- make_target(fits[[best]], threshold = PWQO_TP)
cat("Best fit by AIC:", best, "\n")
print(tgt$params)
cat("Exact P(TP >", PWQO_TP, ") under the fit:", round(tgt$p_exceed, 4), "\n")
saveRDS(tgt, "output/mh_target.rds")
