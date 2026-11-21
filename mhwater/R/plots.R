#' Four-panel diagnostic plot for candidate fits
#'
#' Draws (1) the empirical vs fitted exceedance curve \eqn{P(X > x)} on
#' log-log axes, (2) a log-log Q-Q plot, (3) empirical vs fitted CDFs and
#' (4) a P-P plot. Log axes are needed because environmental concentrations
#' span several orders of magnitude.
#'
#' @param fits named list of `fitdist` objects from [fit_candidates()].
#' @param threshold guideline value drawn as a dashed line. Default 30.
#' @param file optional path to a `.png` file. If `NULL` (default) the plot
#'   is drawn on the current device.
#'
#' @returns `NULL`, invisibly. Called for its side effect.
#' @importFrom graphics par plot lines abline legend
#' @importFrom grDevices png dev.off palette
#' @importFrom fitdistrplus qqcomp cdfcomp ppcomp
#' @examples
#' set.seed(1)
#' fits <- fit_candidates(rlnorm(300, 2.6, 1.4))
#' plot_fit_diagnostics(fits)
#' @export
plot_fit_diagnostics <- function(fits, threshold = 30, file = NULL) {
  if (!is.null(file)) {
    grDevices::png(file, width = 10, height = 8, units = "in", res = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  op <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(op), add = TRUE)

  x    <- fits[[1]]$data
  cols <- seq_along(fits) + 1
  leg  <- names(fits)

  # Panel 1: exceedance curve
  xs <- sort(x)
  emp_surv <- 1 - (seq_along(xs) - 0.5) / length(xs)
  graphics::plot(xs, emp_surv, log = "xy", pch = 16, cex = 0.3, col = "grey40",
                 main = "Exceedance curve P(X > x)", xlab = "x",
                 ylab = "P(X > x)")
  xg <- exp(seq(log(min(xs)), log(max(xs)), length.out = 400))
  for (i in seq_along(fits)) {
    graphics::lines(xg, vapply(xg, function(q) exceed_prob(fits[[i]], q), 0),
                    col = cols[i], lwd = 2)
  }
  graphics::abline(v = threshold, lty = 2)
  graphics::legend("bottomleft", c("Empirical", leg), col = c("grey40", cols),
                   pch = c(16, rep(NA, length(fits))),
                   lty = c(NA, rep(1, length(fits))), bty = "n")

  # Panels 2-4: fitdistrplus comparison plots
  fitdistrplus::qqcomp(fits, legendtext = leg, fitcol = cols,
                       xlogscale = TRUE, ylogscale = TRUE,
                       main = "Q-Q plot (log-log)")
  fitdistrplus::cdfcomp(fits, legendtext = leg, fitcol = cols,
                        xlogscale = TRUE, main = "Empirical vs fitted CDF")
  fitdistrplus::ppcomp(fits, legendtext = leg, fitcol = cols,
                       main = "P-P plot")
  invisible(NULL)
}

#' Fitted densities on the log10 scale
#'
#' Histogram of \eqn{\log_{10}(x)} with each fitted density transformed to
#' that scale using the change of variables
#' \eqn{f_Y(y) = f_X(10^y)\,10^y \ln 10}, so that the curves remain valid
#' densities on the log axis.
#'
#' @param fits named list of `fitdist` objects from [fit_candidates()].
#' @param threshold guideline value drawn as a dashed line. Default 30.
#' @param bins number of histogram bins. Default 60.
#'
#' @returns A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline labs
#'   theme_minimal theme scale_colour_manual after_stat
#' @importFrom tidyr pivot_longer
#' @examples
#' set.seed(1)
#' fits <- fit_candidates(rlnorm(300, 2.6, 1.4))
#' plot_log_density(fits)
#' @export
plot_log_density <- function(fits, threshold = 30, bins = 60) {
  x  <- fits[[1]]$data
  yg <- seq(min(log10(x)), max(log10(x)), length.out = 500)
  dens <- lapply(names(fits), function(nm) {
    f    <- fits[[nm]]
    dfun <- match.fun(paste0("d", f$distname))
    fx   <- do.call(dfun, c(list(x = 10^yg), as.list(f$estimate)))
    tibble::tibble(y = yg, distribution = nm, fy = fx * 10^yg * log(10))
  })
  dens <- dplyr::bind_rows(dens)
  cols <- stats::setNames(grDevices::palette()[seq_along(fits) + 1], names(fits))

  ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = tibble::tibble(lx = log10(x)),
      ggplot2::aes(x = .data$lx, y = ggplot2::after_stat(density)),
      bins = bins, fill = "grey80", colour = "white") +
    ggplot2::geom_line(data = dens,
                       ggplot2::aes(.data$y, .data$fy, colour = .data$distribution),
                       linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = log10(threshold), linetype = "dashed") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "Fitted candidate densities, log10 scale",
                  subtitle = paste("Dashed line: threshold =", threshold),
                  x = "log10(x)", y = "Density", colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}
