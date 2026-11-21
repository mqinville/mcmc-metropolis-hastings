# =============================================================================
# STAT*4000 MCMC Project - Environmental Application (Objective 2)
# Step 1: Clean PWQMN 2025 total phosphorus data
# Step 2: Fit and compare candidate distributions (lognormal, gamma, Weibull)
#
# Run from the "MCMC Research Project" folder (setwd() there, or open it as
# the working directory). Expects: data/PWQMN_2025.csv
#
# Outputs
#   data/tp_clean.csv          cleaned TP observations (one row per sample)
#   output/cleaning_log.csv    row counts after each cleaning step
#   output/fit_comparison.csv  AIC / BIC / GOF statistics + exceedance probs
#   output/tp_fits.rds         fitted parameters -> target for the MH sampler
#   figures/*.png              histogram + 4-panel fit diagnostics
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fitdistrplus)
})

set.seed(4000)
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

PWQO_TP <- 30  # µg/L. Ontario PWQO (interim) for TP in rivers & streams = 0.03 mg/L

# -----------------------------------------------------------------------------
# 1. IMPORT
# -----------------------------------------------------------------------------
# The file is Latin-1 encoded (the "µ" in µg/L breaks under UTF-8).
# Read every column as character first so "<0.2" style results survive.
raw <- read_csv("data/PWQMN_2025.csv",
                locale = locale(encoding = "latin1"),
                col_types = cols(.default = col_character()))

log_step <- function(df, step) tibble(step = step, n_rows = nrow(df))
clean_log <- log_step(raw, "Raw file (all analytes)")

# -----------------------------------------------------------------------------
# 2. KEEP TOTAL PHOSPHORUS ONLY
# -----------------------------------------------------------------------------
# "Phosphorus; total" (µg/L, method E3558) is the analyte named in the
# proposal. "Phosphorus; phosphate" is a different analyte (mg/L, ~12% BDL).
tp <- raw |>
  filter(Analyte == "Phosphorus; total") |>
  rename(site = `Collection Site`, date = `Collection Date`,
         time = `Collection Time`, result_raw = Result,
         result_call = `Result Call`, det_limit = `Detection Limit`,
         lab_sample = `Lab Sample ID`, remark = `Remark Code`) |>
  mutate(site = sub("^'", "", site))   # site IDs carry a leading apostrophe
clean_log <- bind_rows(clean_log, log_step(tp, "Total phosphorus rows"))

# Every TP row is in µg/L (checked without typing "µ", which is locale-fragile)
stopifnot(length(unique(tp$Units)) == 1, grepl("g/L$", tp$Units[1]),
          !grepl("^mg", tp$Units[1]))

# -----------------------------------------------------------------------------
# 3. REMOVE DUPLICATE RECORDS
# -----------------------------------------------------------------------------
# 74 site/date/time combinations appear twice with an identical result and
# lab sample ID: these are repeated rows, not field replicates.
tp <- tp |> distinct(site, date, time, result_raw, .keep_all = TRUE)
clean_log <- bind_rows(clean_log, log_step(tp, "After removing duplicate rows"))

# Safety check: no remaining sample has two different results
stopifnot(!any(duplicated(tp[, c("site", "date", "time")])))

# -----------------------------------------------------------------------------
# 4. PARSE RESULTS AND HANDLE NON-DETECTS
# -----------------------------------------------------------------------------
# Below-detection-limit (BDL) values are stored as text, e.g. "<0.2".
# Only 20 of ~3,700 samples (~0.5%) are BDL, so we substitute DL/2 for the
# main analysis and check that a censored-likelihood fit gives the same answer.
tp <- tp |>
  mutate(
    censored   = result_call == "BDL",
    det_limit  = as.numeric(det_limit),
    result_num = suppressWarnings(as.numeric(result_raw)),
    tp_ugL     = if_else(censored, det_limit / 2, result_num),
    date       = as.Date(date),
    month      = factor(format(date, "%b"), levels = month.abb),
    remark_flag = !is.na(remark)   # lab remark codes kept but flagged
  )

stopifnot(!any(is.na(tp$tp_ugL)), all(tp$tp_ugL > 0))

tp_clean <- tp |>
  dplyr::select(site, date, time, month, tp_ugL, censored, det_limit,
         remark, remark_flag, lab_sample)
clean_log <- bind_rows(clean_log, log_step(tp_clean, "Final clean TP samples"))

write_csv(tp_clean, "data/tp_clean.csv")
write_csv(clean_log, "output/cleaning_log.csv")
print(clean_log)

x <- tp_clean$tp_ugL

# -----------------------------------------------------------------------------
# 5. DESCRIPTIVE SUMMARY
# -----------------------------------------------------------------------------
desc <- tibble(
  n = length(x), n_sites = n_distinct(tp_clean$site),
  n_BDL = sum(tp_clean$censored), n_remark = sum(tp_clean$remark_flag),
  mean = mean(x), sd = sd(x), median = median(x),
  q90 = quantile(x, .90), q99 = quantile(x, .99), max = max(x),
  skewness = mean((x - mean(x))^3) / sd(x)^3,
  emp_P_exceed_PWQO = mean(x > PWQO_TP)
)
print(desc, width = Inf)
write_csv(desc, "output/tp_descriptive.csv")

# Raw and log-scale histograms: the raw scale is extremely right-skewed,
# log(TP) looks roughly symmetric, which already hints at a lognormal.
p_hist <- ggplot(tp_clean, aes(tp_ugL)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = PWQO_TP, linetype = "dashed", colour = "red") +
  scale_x_log10() +
  labs(title = "PWQMN 2025 total phosphorus (log scale)",
       subtitle = "Dashed line: PWQO 30 µg/L",
       x = "Total phosphorus (µg/L, log10 axis)", y = "Count") +
  theme_minimal()
ggsave("figures/tp_histogram_log.png", p_hist, width = 7, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# 6. FIT CANDIDATE DISTRIBUTIONS (maximum likelihood)
# -----------------------------------------------------------------------------
# All three are positive, right-skewed families (proposal: gamma + lognormal;
# Weibull added as a common third environmental candidate).
fit_ln <- fitdist(x, "lnorm")
# Gamma needs a sensible start (default moment start is poor for heavy tails).
fit_ga <- fitdist(x, "gamma", start = list(shape = 0.5, rate = 0.5 / mean(x)))
fit_wb <- fitdist(x, "weibull")
fits   <- list(lognormal = fit_ln, gamma = fit_ga, weibull = fit_wb)

lapply(fits, summary)

# -----------------------------------------------------------------------------
# 7. COMPARE FITS
# -----------------------------------------------------------------------------
gof <- gofstat(fits, fitnames = names(fits))

# Model-implied exceedance probability P(X > 30) under each fit
p_exceed <- c(
  lognormal = plnorm(PWQO_TP, fit_ln$estimate["meanlog"], fit_ln$estimate["sdlog"],
                     lower.tail = FALSE),
  gamma     = pgamma(PWQO_TP, fit_ga$estimate["shape"], fit_ga$estimate["rate"],
                     lower.tail = FALSE),
  weibull   = pweibull(PWQO_TP, fit_wb$estimate["shape"], fit_wb$estimate["scale"],
                       lower.tail = FALSE)
)

# Model-implied upper quantiles vs empirical (tail fit matters for exceedance)
q_model <- sapply(fits, function(f) quantile(f, probs = c(.5, .9, .99))$quantiles)

comparison <- tibble(
  distribution = names(fits),
  param_1 = sapply(fits, \(f) sprintf("%s = %.4g", names(f$estimate)[1], f$estimate[1])),
  param_2 = sapply(fits, \(f) sprintf("%s = %.4g", names(f$estimate)[2], f$estimate[2])),
  logLik  = sapply(fits, \(f) f$loglik),
  AIC     = gof$aic, BIC = gof$bic,
  KS      = gof$ks, CvM = gof$cvm, AD = gof$ad,
  P_exceed_30 = p_exceed,
  q50 = unlist(q_model[1, ]), q90 = unlist(q_model[2, ]), q99 = unlist(q_model[3, ])
) |>
  arrange(AIC) |>
  mutate(delta_AIC = AIC - min(AIC))

empirical_row <- tibble(distribution = "EMPIRICAL",
                        P_exceed_30 = mean(x > PWQO_TP),
                        q50 = quantile(x, .5), q90 = quantile(x, .9),
                        q99 = quantile(x, .99))
comparison_out <- bind_rows(comparison, empirical_row)
print(comparison_out, width = Inf)
write_csv(comparison_out, "output/fit_comparison.csv")

# NOTE on GOF stats: AD = Inf for gamma/Weibull because their fitted CDF
# rounds to exactly 1 at the largest observations (10,300 ug/L) -> log(0).
# That is itself evidence their tails are far too light.

# Diagnostic plots. Log axes are needed, otherwise everything is squashed at 0.
png("figures/tp_fit_diagnostics.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))
leg <- c("Lognormal", "Gamma", "Weibull")
# Panel 1: exceedance curve P(X > x), log-log. This is what matters for the
# PWQO exceedance question, and shows tail behaviour clearly.
xs <- sort(x); emp_surv <- 1 - (seq_along(xs) - 0.5) / length(xs)
plot(xs, emp_surv, log = "xy", pch = 16, cex = 0.3, col = "grey40",
     main = "Exceedance curve P(X > x)", xlab = "TP (\u00b5g/L)", ylab = "P(X > x)")
xg2 <- exp(seq(log(min(xs)), log(max(xs)), length.out = 400))
lines(xg2, plnorm(xg2, fit_ln$estimate[1], fit_ln$estimate[2], lower.tail = FALSE), col = 2, lwd = 2)
lines(xg2, pgamma(xg2, fit_ga$estimate[1], fit_ga$estimate[2], lower.tail = FALSE), col = 3, lwd = 2)
lines(xg2, pweibull(xg2, fit_wb$estimate[1], fit_wb$estimate[2], lower.tail = FALSE), col = 4, lwd = 2)
abline(v = PWQO_TP, lty = 2)
legend("bottomleft", c("Empirical", leg), col = c("grey40", 2:4),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), bty = "n")
qqcomp(fits, legendtext = leg, fitcol = 2:4, xlogscale = TRUE, ylogscale = TRUE,
       main = "Q-Q plot (log-log)")
cdfcomp(fits, legendtext = leg, fitcol = 2:4, xlogscale = TRUE, main = "Empirical vs fitted CDF",
        xlab = "TP (\u00b5g/L)")
ppcomp(fits, legendtext = leg, fitcol = 2:4, main = "P-P plot")
dev.off()

# Fitted densities on the log10 scale (whole range visible).
# If Y = log10(X) then f_Y(y) = f_X(10^y) * 10^y * log(10).
ygrid <- seq(min(log10(x)), max(log10(x)), length.out = 500)
xg    <- 10^ygrid
dens_log <- tibble(
  y = ygrid,
  Lognormal = dlnorm(xg, fit_ln$estimate[1], fit_ln$estimate[2]),
  Gamma     = dgamma(xg, fit_ga$estimate[1], fit_ga$estimate[2]),
  Weibull   = dweibull(xg, fit_wb$estimate[1], fit_wb$estimate[2])
) |>
  pivot_longer(-y, names_to = "distribution", values_to = "fx") |>
  mutate(fy = fx * 10^y * log(10))

p_dens <- ggplot() +
  geom_histogram(aes(x = log10(x), y = after_stat(density)), bins = 60,
                 fill = "grey80", colour = "white") +
  geom_line(data = dens_log, aes(y, fy, colour = distribution), linewidth = 0.9) +
  geom_vline(xintercept = log10(PWQO_TP), linetype = "dashed") +
  labs(title = "Fitted candidate densities, log10 scale",
       subtitle = "Dashed line: PWQO 30 \u00b5g/L",
       x = "log10(TP, \u00b5g/L)", y = "Density", colour = NULL) +
  scale_colour_manual(values = c(Lognormal = palette()[2], Gamma = palette()[3],
                                 Weibull = palette()[4])) +
  theme_minimal() + theme(legend.position = "top")
ggsave("figures/tp_fitted_densities_log.png", p_dens, width = 7, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. SENSITIVITY CHECKS
# -----------------------------------------------------------------------------
# (a) Treat the 20 BDL values as left-censored instead of DL/2 substitution.
cens_df <- data.frame(
  left  = ifelse(tp_clean$censored, NA, tp_clean$tp_ugL),
  right = ifelse(tp_clean$censored, tp_clean$det_limit, tp_clean$tp_ugL)
)
fit_ln_cens <- fitdistcens(cens_df, "lnorm")
sens_cens <- rbind(DL_half = fit_ln$estimate, censored_MLE = fit_ln_cens$estimate)
print(round(sens_cens, 4))

# (b) Samples from the same site are not independent (~8 per site).
#     Refit using one value per site (site median) and see whether the ranking
#     of the candidate families changes.
site_med <- tp_clean |> group_by(site) |> summarise(tp = median(tp_ugL)) |> pull(tp)
fits_site <- list(
  lognormal = fitdist(site_med, "lnorm"),
  gamma     = fitdist(site_med, "gamma",
                      start = list(shape = 0.5, rate = 0.5 / mean(site_med))),
  weibull   = fitdist(site_med, "weibull")
)
gof_site <- gofstat(fits_site, fitnames = names(fits_site))
sens_site <- tibble(distribution = names(fits_site), AIC = gof_site$aic,
                    BIC = gof_site$bic, AD = gof_site$ad) |> arrange(AIC)
print(sens_site)

sensitivity <- list(censoring = sens_cens, site_median = sens_site)

# -----------------------------------------------------------------------------
# 9. SAVE TARGET FOR THE METROPOLIS-HASTINGS SAMPLER
# -----------------------------------------------------------------------------
# The MH code can load this and use the best fitted density as its target:
#   tgt <- readRDS("output/tp_fits.rds")
#   dtarget <- function(x) dlnorm(x, tgt$best$meanlog, tgt$best$sdlog)
# The analytic P(X > 30) under the fit is the "truth" the MH estimate
# should reproduce.
best_name <- comparison$distribution[1]
saveRDS(list(
  best          = as.list(fits[[best_name]]$estimate),
  best_family   = best_name,
  all_estimates = lapply(fits, \(f) as.list(f$estimate)),
  vcov          = lapply(fits, \(f) f$vcov),
  PWQO_ugL      = PWQO_TP,
  P_exceed_true = p_exceed,
  emp_P_exceed  = mean(x > PWQO_TP),
  n             = length(x),
  sensitivity   = sensitivity
), "output/tp_fits.rds")

cat("\nBest-fitting family by AIC:", best_name, "\n")
sessionInfo()
