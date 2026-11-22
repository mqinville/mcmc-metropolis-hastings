#' Metropolis-Hastings sampler
#'
#' Draws a Markov chain whose stationary distribution is the target density
#' \eqn{f}. At each step a candidate \eqn{Y} is drawn from the proposal
#' \eqn{g(\cdot \mid X_t)} and accepted with probability
#' \eqn{\min(1, \alpha)}, where (Rizzo 2019, Section 11.2.1, eq. 11.3)
#' \deqn{\alpha = \frac{f(Y)\, g(X_t \mid Y)}{f(X_t)\, g(Y \mid X_t)}.}
#' When the candidate is rejected the chain stays at \eqn{X_t}.
#'
#' The ratio is computed on the log scale to avoid underflow for tiny
#' densities:
#' `log_alpha = (log f(Y) - log f(X_t)) + (log g(X_t | Y) - log g(Y | X_t))`,
#' and \eqn{Y} is accepted if `log(U) <= log_alpha` with
#' \eqn{U \sim Uniform(0, 1)}. For a symmetric proposal such as
#' [normal_rw()], \eqn{g(X_t \mid Y) = g(Y \mid X_t)}, so the correction term
#' cancels and this is the random walk Metropolis sampler (Rizzo 2019,
#' Sections 11.2.2-11.2.3). Asymmetric proposals ([lognormal_rw()]) and
#' independence proposals ([gamma_independence()]) rely on the correction.
#'
#' Candidates where `log_target` is `-Inf`, `NA` or `NaN` (e.g. outside the
#' support) are always rejected.
#'
#' @param log_target function of one number returning the log target density
#'   (e.g. `make_target(fit)$log_density`). It may return `-Inf` outside the
#'   support.
#' @param n_iter total chain length, **including** the starting value (so
#'   `x[1] == x0`). A single whole number, at least 2.
#' @param x0 starting value; `log_target(x0)` must be finite.
#' @param proposal the proposal distribution: a list with functions
#'   `draw(x)` (one candidate given the current state `x`) and
#'   `log_q(to, from)` (log density of proposing `to` from `from`), as
#'   returned by [normal_rw()], [lognormal_rw()] or [gamma_independence()].
#'
#' @returns A list with elements `x` (numeric vector of length `n_iter`, the
#'   chain), `accept_rate` (proportion of the `n_iter - 1` candidates that
#'   were accepted) and `proposal` (the proposal's `name`, or
#'   `"custom proposal"` if it has none).
#' @references Rizzo, M. L. (2019). *Statistical Computing with R* (2nd ed.),
#'   Chapter 11. Chapman & Hall/CRC.
#' @importFrom stats runif
#' @examples
#' set.seed(1)
#' chain <- mh_sampler(function(x) dnorm(x, log = TRUE),
#'                     n_iter = 5000, x0 = 0, proposal = normal_rw(2.4))
#' chain$accept_rate
#' mean(chain$x[-(1:500)])        # close to 0
#'
#' fits <- fit_candidates(tp_clean$value)
#' tgt  <- make_target(fits$lnorm, threshold = 30)
#' chain <- mh_sampler(tgt$log_density, n_iter = 10000, x0 = 10,
#'                     proposal = normal_rw(20))
#' mean(chain$x[-(1:1000)] > tgt$threshold)   # compare with tgt$p_exceed
#'
#' # Multiplicative random walk: suits the heavy right tail of the lognormal
#' chain <- mh_sampler(tgt$log_density, n_iter = 10000, x0 = 10,
#'                     proposal = lognormal_rw(1.5))
#' chain$accept_rate
#' mean(chain$x[-(1:1000)] > tgt$threshold)
#' @export
mh_sampler <- function(log_target, n_iter, x0, proposal) {
  if (!is.numeric(n_iter) || length(n_iter) != 1 || !is.finite(n_iter) || n_iter < 2 || n_iter != round(n_iter)) {
    print("n_iter must be a single whole number, at least 2.")
    return(NULL)
  }

  if (!is.list(proposal) || !is.function(proposal$draw) || !is.function(proposal$log_q)) {
    print("proposal must be a list with functions draw(x) and log_q(to, from), ", "e.g. normal_rw(1).")
    return(NULL)
  }

  current_logp <- log_target(x0)
  if (!is.finite(current_logp)) {
    print("log_target(x0) is not finite: choose a starting value inside the ", "support of the target.")
    return(NULL)
  }

  x <- numeric(n_iter)
  x[1] <- x0
  n_accept <- 0

  for (i in 2:n_iter) {
    x_t <- x[i - 1]
    y <- proposal$draw(x_t)

    # A NA/NaN target value counts as zero density: the candidate is rejected.
    y_logp <- log_target(y)
    if (is.na(y_logp)) {
      y_logp <- -Inf
    }

    # Log of the MH ratio (eq. 11.3). No need for min(1, ratio):
    # log(u) <= 0 always.
    log_alpha <- (y_logp - current_logp) + proposal$log_q(x_t, y) - proposal$log_q(y, x_t)
    if (is.na(log_alpha)) {
      log_alpha <- -Inf
    }
    
    if (log(stats::runif(1)) <= log_alpha) {
      x[i] <- y
      current_logp <- y_logp
      n_accept <- n_accept + 1
    } else {
      x[i] <- x_t
    }
  }

  proposal_name <- proposal$name
  if (is.null(proposal_name)) {
    proposal_name <- "custom proposal"
  }

  list(x = x, accept_rate = n_accept / (n_iter - 1), proposal = proposal_name)
}
