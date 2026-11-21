#' Sensitivity check: treat below-detection values as left-censored
#'
#' Refits a distribution with [fitdistrplus::fitdistcens()], treating BDL
#' samples as left-censored at their detection limit instead of substituting
#' a value. Returns both sets of estimates side by side so the effect of the
#' DL/2 substitution in [clean_tp()] can be judged.
#'
#' @param data cleaned data, the `data` element of [clean_tp()] output.
#' @param family distribution name. Default `"lnorm"`.
#'
#' @returns A matrix with rows `substituted` and `censored_MLE`.
#' @importFrom fitdistrplus fitdistcens
#' @examples
#' set.seed(2)
#' v <- rlnorm(300, 2.6, 1.4)
#' d <- data.frame(value = pmax(v, 0.1), censored = v < 0.2, det_limit = 0.2)
#' censoring_check(d)
#' @export
censoring_check <- function(data, family = "lnorm") {
  cens_df <- data.frame(
    left  = ifelse(data$censored, NA, data$value),
    right = ifelse(data$censored, data$det_limit, data$value)
  )
  sub  <- fitdistrplus::fitdist(data$value, family)
  cens <- fitdistrplus::fitdistcens(cens_df, family)
  rbind(substituted = sub$estimate, censored_MLE = cens$estimate)
}

#' Sensitivity check: one value per site
#'
#' Repeated samples from the same site are not independent. This refits the
#' candidate families to the per-site median so the ranking can be compared
#' with the pooled fit.
#'
#' @param data cleaned data with columns `site` and `value`.
#' @param families candidate families, passed to [fit_candidates()].
#'
#' @returns A tibble with AIC, BIC and AD for each family, sorted by AIC.
#' @importFrom dplyr group_by summarise pull arrange
#' @examples
#' set.seed(3)
#' d <- data.frame(site = rep(1:60, each = 5), value = rlnorm(300, 2.6, 1.4))
#' site_median_check(d)
#' @export
site_median_check <- function(data, families = c("lnorm", "gamma", "weibull")) {
  med <- data |>
    dplyr::group_by(.data$site) |>
    dplyr::summarise(m = stats::median(.data$value)) |>
    dplyr::pull(.data$m)
  fits <- fit_candidates(med, families)
  gof  <- fitdistrplus::gofstat(fits, fitnames = names(fits))
  tibble::tibble(distribution = names(fits), n_sites = length(med),
                 AIC = unname(gof$aic), BIC = unname(gof$bic),
                 AD = unname(gof$ad)) |>
    dplyr::arrange(.data$AIC)
}
