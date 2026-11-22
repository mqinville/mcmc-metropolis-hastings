#' Proposal distributions for [mh_sampler()]
#'
#' Constructors for the proposal distribution \eqn{g(\cdot \mid X_t)} used by
#' [mh_sampler()]. Each returns a list of class `"mh_proposal"` with
#' elements:
#' * `draw(x)`: draws one candidate \eqn{Y} given the current state `x`;
#' * `q(to, from)`: density of proposing `to` when the chain is at `from`,
#'   i.e. \eqn{g(to \mid from)};
#' * `name`: a label such as `"normal_rw(sd = 2)"` for plots and tables.
#'
#' * `normal_rw(sd)`: symmetric random walk \eqn{Y \sim N(X_t, sd^2)}. The
#'   Hastings correction cancels, giving the random walk Metropolis sampler
#'   (Rizzo 2019, Section 11.2.3).
#' * `lognormal_rw(sd)`: multiplicative random walk
#'   \eqn{\log Y \sim N(\log X_t, sd^2)}. Proposals are always positive, so
#'   use it only for targets on \eqn{(0, \infty)}, with `x0 > 0`. It is not
#'   symmetric, so it relies on the Hastings correction.
#' * `gamma_independence(shape, rate)`: independence sampler
#'   \eqn{Y \sim Gamma(shape, rate)} regardless of \eqn{X_t}
#'   (Rizzo 2019, Section 11.2.4). Positive targets only; works best when the
#'   gamma is close to the target with tails at least as heavy.
#'
#' Any list with functions `draw` and `q` (and optionally `name`) can be
#' passed to [mh_sampler()] as a custom proposal.
#'
#' @param sd standard deviation of the step (on the log scale for
#'   `lognormal_rw()`). Must be > 0.
#' @param shape,rate gamma shape and rate. Both must be > 0.
#'
#' @returns A list of class `"mh_proposal"` with elements `draw`, `q`
#'   and `name`.
#' @importFrom stats rnorm dnorm rlnorm dlnorm rgamma dgamma
#' @examples
#' p <- normal_rw(2)
#' p$name
#' p$draw(0)
#' p$q(1, 0)
#'
#' lognormal_rw(0.5)$draw(10)        # always positive
#' gamma_independence(2, 0.5)$draw(10)  # ignores the current state
#' @name proposals
NULL

#' @rdname proposals
#' @export
normal_rw <- function(sd) {
  check_positive(sd, "sd")
  structure(
    list(
      draw  = function(x) stats::rnorm(1, mean = x, sd = sd),
      q     = function(to, from) stats::dnorm(to, mean = from, sd = sd),
      name  = sprintf("normal_rw(sd = %g)", sd)
    ),
    class = "mh_proposal"
  )
}

#' @rdname proposals
#' @export
lognormal_rw <- function(sd) {
  check_positive(sd, "sd")
  structure(
    list(
      draw  = function(x) {
        if (x <= 0) {
          stop("lognormal_rw() needs a positive current state.")
        }
        stats::rlnorm(1, meanlog = log(x), sdlog = sd)
      },
      q     = function(to, from) stats::dlnorm(to, meanlog = log(from),
                                               sdlog = sd),
      name  = sprintf("lognormal_rw(sd = %g)", sd)
    ),
    class = "mh_proposal"
  )
}

#' @rdname proposals
#' @export
gamma_independence <- function(shape, rate) {
  check_positive(shape, "shape")
  check_positive(rate, "rate")
  structure(
    list(
      draw  = function(x) stats::rgamma(1, shape = shape, rate = rate),
      q     = function(to, from) stats::dgamma(to, shape = shape, rate = rate),
      name  = sprintf("gamma_independence(shape = %g, rate = %g)",
                      shape, rate)
    ),
    class = "mh_proposal"
  )
}

# Stop unless `value` is a single finite number > 0.
check_positive <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value) ||
      value <= 0) {
    stop(arg, " must be a single positive number.")
  }
}
