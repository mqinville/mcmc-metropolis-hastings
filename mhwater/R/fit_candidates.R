#' Fit candidate distributions by maximum likelihood
#'
#' Fits each candidate family to a vector of positive measurements using
#' [fitdistrplus::fitdist()]. The gamma fit is given a start of
#' `shape = 0.5, rate = 0.5 / mean(x)` because the default moment-based start
#' fails for very heavy-tailed data.
#'
#' @param x numeric vector of positive values.
#' @param families character vector of distribution names understood by
#'   `fitdist()` (`"lnorm"`, `"gamma"`, `"weibull"`, `"exp"`, ...).
#'   Default: lognormal, gamma and Weibull.
#'
#' @returns A named list of `fitdist` objects (names = `families`).
#' @importFrom fitdistrplus fitdist
#' @examples
#' set.seed(1)
#' x <- rlnorm(500, meanlog = 2.6, sdlog = 1.4)
#' fits <- fit_candidates(x)
#' fits$lnorm$estimate
#' @export
fit_candidates <- function(x, families = c("lnorm", "gamma", "weibull")) {
  if (any(x <= 0, na.rm = TRUE)) stop("All values of x must be positive.")
  x <- x[!is.na(x)]
  fits <- lapply(families, function(fam) {
    if (fam == "gamma") {
      fitdistrplus::fitdist(x, "gamma",
                            start = list(shape = 0.5, rate = 0.5 / mean(x)))
    } else {
      fitdistrplus::fitdist(x, fam)
    }
  })
  stats::setNames(fits, families)
}

#' Exceedance probability under a fitted distribution
#'
#' Computes \eqn{P(X > threshold)} using the fitted parameters, i.e. the
#' "true" exceedance probability that a Metropolis-Hastings estimate should
#' reproduce when this fit is used as the target.
#'
#' @param fit a `fitdist` object.
#' @param threshold numeric threshold. Default 30 (µg/L).
#'
#' @returns A single number.
#' @examples
#' fit <- fitdistrplus::fitdist(rlnorm(200, 2.6, 1.4), "lnorm")
#' exceed_prob(fit, 30)
#' @export
exceed_prob <- function(fit, threshold = 30) {
  pfun <- match.fun(paste0("p", fit$distname))
  do.call(pfun, c(list(q = threshold, lower.tail = FALSE),
                  as.list(fit$estimate)))
}

#' Compare fitted candidate distributions
#'
#' Builds one comparison table with log-likelihood, AIC, BIC, the
#' Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling statistics,
#' model-implied exceedance probability and quantiles, sorted by AIC.
#' An extra row `"empirical"` gives the observed exceedance proportion and
#' sample quantiles for reference.
#'
#' Note: the Anderson-Darling statistic is `Inf` when a fitted CDF rounds to
#' exactly 1 at the largest observation. That itself indicates a tail that
#' is far too light.
#'
#' @param fits named list of `fitdist` objects, e.g. from [fit_candidates()].
#' @param threshold guideline value for exceedance. Default 30.
#' @param probs quantile levels to report. Default 0.5, 0.9, 0.99.
#'
#' @returns A tibble with one row per distribution plus an `"empirical"` row.
#' @importFrom fitdistrplus gofstat
#' @importFrom dplyr arrange mutate bind_rows
#' @examples
#' set.seed(1)
#' x <- rlnorm(500, 2.6, 1.4)
#' compare_fits(fit_candidates(x))
#' @export
compare_fits <- function(fits, threshold = 30, probs = c(0.5, 0.9, 0.99)) {
  x   <- fits[[1]]$data
  gof <- fitdistrplus::gofstat(fits, fitnames = names(fits))

  qmat <- t(sapply(fits, function(f) {
    unlist(stats::quantile(f, probs = probs)$quantiles)
  }))
  colnames(qmat) <- paste0("q", probs * 100)

  par_str <- function(f, i) {
    sprintf("%s = %.4g", names(f$estimate)[i], f$estimate[i])
  }

  out <- tibble::tibble(
    distribution = names(fits),
    param_1  = vapply(fits, par_str, "", i = 1),
    param_2  = vapply(fits, function(f) {
      if (length(f$estimate) > 1) par_str(f, 2) else NA_character_
    }, ""),
    logLik   = vapply(fits, function(f) f$loglik, 0),
    AIC      = unname(gof$aic), BIC = unname(gof$bic),
    KS       = unname(gof$ks), CvM = unname(gof$cvm), AD = unname(gof$ad),
    P_exceed = vapply(fits, exceed_prob, 0, threshold = threshold)
  )
  out <- cbind(out, qmat) |>
    tibble::as_tibble() |>
    dplyr::arrange(.data$AIC) |>
    dplyr::mutate(delta_AIC = .data$AIC - min(.data$AIC))

  emp <- tibble::as_tibble(as.list(stats::quantile(x, probs)))
  names(emp) <- colnames(qmat)
  emp <- cbind(tibble::tibble(distribution = "empirical",
                              P_exceed = mean(x > threshold)), emp)
  dplyr::bind_rows(out, tibble::as_tibble(emp))
}

#' Build a Metropolis-Hastings target from a fitted distribution
#'
#' Wraps a `fitdist` object into the pieces a Metropolis-Hastings sampler
#' needs: the target density, its log, and the exact exceedance probability
#' to compare the MCMC estimate against.
#'
#' @param fit a `fitdist` object (usually the best fit by AIC).
#' @param threshold guideline value for exceedance. Default 30.
#'
#' @returns A list with elements `family`, `params` (named numeric),
#'   `density` (function of x, the raw-scale target to pass to
#'   [mh_sampler()]), `log_density` (function of x),
#'   `threshold` and `p_exceed`.
#' @examples
#' fit <- fitdistrplus::fitdist(rlnorm(200, 2.6, 1.4), "lnorm")
#' tgt <- make_target(fit)
#' tgt$density(10)
#' tgt$p_exceed
#' @export
make_target <- function(fit, threshold = 30) {
  dfun   <- match.fun(paste0("d", fit$distname))
  params <- fit$estimate
  list(
    family      = fit$distname,
    params      = params,
    density     = function(x) do.call(dfun, c(list(x = x), as.list(params))),
    log_density = function(x) do.call(dfun, c(list(x = x, log = TRUE),
                                              as.list(params))),
    threshold   = threshold,
    p_exceed    = exceed_prob(fit, threshold)
  )
}
