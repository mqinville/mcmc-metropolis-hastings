#' Metropolis-Hastings sampler
#'
#' Draws a Markov chain whose stationary distribution is the target density
#' \eqn{f}. At each step a candidate \eqn{Y} is drawn from the proposal
#' \eqn{g(\cdot \mid X_t)} and accepted when \eqn{U \le r(X_t, Y)}, with
#' \eqn{U \sim Uniform(0, 1)} and (Rizzo 2019, Section 11.2.1, eq. 11.3)
#' \deqn{r(X_t, Y) = \frac{f(Y)\, g(X_t \mid Y)}{f(X_t)\, g(Y \mid X_t)}.}
#' When the candidate is rejected the chain stays at \eqn{X_t}. Because
#' \eqn{U \le 1}, comparing \eqn{U} with \eqn{r} is the same as accepting with
#' probability \eqn{\min(1, r)}.
#'
#' For a symmetric proposal such as [normal_rw()],
#' \eqn{g(X_t \mid Y) = g(Y \mid X_t)}, so the correction cancels and
#' \eqn{r = f(Y) / f(X_t)}: this is the random walk Metropolis sampler
#' (Rizzo 2019, Sections 11.2.2-11.2.3). Asymmetric proposals
#' ([lognormal_rw()]) and independence proposals ([gamma_independence()])
#' rely on the correction.
#'
#' Candidates where `target` is 0, `NA` or `NaN` (e.g. outside the support)
#' are always rejected.
#'
#' @details The ratio is computed on the raw density scale, as in Rizzo
#'   (2019). Densities can underflow (round to 0) in the far tails or for
#'   targets that are products of many small terms. If `target(x0)` is 0 for
#'   that reason, move the starting value toward the bulk of the target.
#'
#' @param target function of one number returning the target density
#'   \eqn{f(x)} (e.g. `make_target(fit)$density`). It may be unnormalised and
#'   should return 0 outside the support.
#' @param n_iter total chain length, **including** the starting value (so
#'   `x[1] == x0`). A single whole number, at least 2.
#' @param x0 starting value; `target(x0)` must be finite and > 0.
#' @param proposal the proposal distribution: a list with functions
#'   `draw(x)` (one candidate given the current state `x`) and
#'   `q(to, from)` (density of proposing `to` from `from`), as returned by
#'   [normal_rw()], [lognormal_rw()] or [gamma_independence()].
#'
#' @returns A list with elements `x` (numeric vector of length `n_iter`, the
#'   chain), `accept_rate` (proportion of the `n_iter - 1` candidates that
#'   were accepted) and `proposal` (the proposal's `name`, or
#'   `"custom proposal"` if it has none). If an input is invalid, a message
#'   is printed and `NULL` is returned.
#' @references Rizzo, M. L. (2019). *Statistical Computing with R* (2nd ed.),
#'   Chapter 11. Chapman & Hall/CRC.
#' @importFrom stats runif
#' @examples
#' set.seed(1)
#' chain <- mh_sampler(function(x) dnorm(x),
#'                     n_iter = 5000, x0 = 0, proposal = normal_rw(2.4))
#' chain$accept_rate
#' mean(chain$x[-(1:500)])        # close to 0
#'
#' fits <- fit_candidates(tp_clean$value)
#' tgt  <- make_target(fits$lnorm, threshold = 30)
#' chain <- mh_sampler(tgt$density, n_iter = 10000, x0 = 10,
#'                     proposal = normal_rw(20))
#' mean(chain$x[-(1:1000)] > tgt$threshold)   # compare with tgt$p_exceed
#'
#' # Multiplicative random walk: suits the heavy right tail of the lognormal
#' chain <- mh_sampler(tgt$density, n_iter = 10000, x0 = 10,
#'                     proposal = lognormal_rw(1.5))
#' chain$accept_rate
#' mean(chain$x[-(1:1000)] > tgt$threshold)
#' @export
mh_sampler <- function(target, n_iter, x0, proposal) {
  if (!is.numeric(n_iter) || length(n_iter) != 1 || !is.finite(n_iter) ||
      n_iter != round(n_iter)) {
    print("n_iter must be a single whole number, at least 2.")
    return(NULL)
  }
  if (n_iter < 2) {
    print("n_iter must be at least 2.")
    return(NULL)
  }

  if (!is.list(proposal) || !is.function(proposal$draw) ||
      !is.function(proposal$q)) {
    print(paste0("proposal must be a list with functions draw(x) and ",
                 "q(to, from), e.g. normal_rw(1)."))
    return(NULL)
  }

  current_f <- target(x0)
  if (length(current_f) != 1 || !is.finite(current_f) || current_f <= 0) {
    print(paste0("target(x0) is not finite and > 0: choose a starting value ",
                 "inside the support of the target (and away from the far ",
                 "tails)."))
    return(NULL)
  }

  x <- numeric(n_iter)
  x[1] <- x0
  n_accept <- 0

  for (i in 2:n_iter) {
    x_t <- x[i - 1]
    y <- proposal$draw(x_t)

    # A NA/NaN target value counts as zero density: the candidate is rejected.
    f_y <- target(y)
    if (is.na(f_y)) {
      f_y <- 0
    }

    # MH ratio (eq. 11.3). A zero proposal density in the denominator, or an
    # undefined ratio, rejects the candidate.
    q_forward <- proposal$q(y, x_t)
    if (is.na(q_forward) || q_forward == 0) {
      r <- 0
    } else {
      r <- (f_y * proposal$q(x_t, y)) / (current_f * q_forward)
      if (is.na(r)) {
        r <- 0
      }
    }

    # No need for min(1, r): u <= 1 always.
    u <- stats::runif(1)
    if (u <= r) {
      x[i] <- y
      current_f <- f_y
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
