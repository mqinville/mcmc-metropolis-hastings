#' Random-walk Metropolis-Hastings sampler
#'
#' Draws a Markov chain whose stationary distribution is the target density.
#' At each step a proposal `y ~ Normal(current, proposal_sd)` is drawn. The
#' normal proposal is symmetric, so the Hastings ratio reduces to
#' `target(y) / target(current)`. The ratio is computed on the log scale to
#' avoid underflow for tiny densities. When a proposal is rejected the chain
#' stays at its current value.
#'
#' @param log_target function of one number returning the log target density
#'   (e.g. `make_target(fit)$log_density`). It may return `-Inf` outside the
#'   support; such proposals are always rejected.
#' @param n_iter total chain length, **including** the starting value (so
#'   `x[1] == x0`). Must be at least 2.
#' @param x0 starting value; `log_target(x0)` must be finite.
#' @param proposal_sd standard deviation of the normal proposal (step size),
#'   must be > 0.
#'
#' @returns A list with elements `x` (numeric vector of length `n_iter`, the
#'   chain) and `accept_rate` (proportion of the `n_iter - 1` proposals that
#'   were accepted).
#' @importFrom stats rnorm runif
#' @examples
#' set.seed(1)
#' chain <- mh_sampler(function(x) dnorm(x, log = TRUE),
#'                     n_iter = 5000, x0 = 0, proposal_sd = 2.4)
#' chain$accept_rate
#' mean(chain$x[-(1:500)])        # close to 0
#'
#' fits <- fit_candidates(tp_clean$value)
#' tgt  <- make_target(fits$lnorm, threshold = 30)
#' chain <- mh_sampler(tgt$log_density, n_iter = 10000, x0 = 10, proposal_sd = 20)
#' mean(chain$x[-(1:1000)] > tgt$threshold)   # compare with tgt$p_exceed
#' @export
mh_sampler <- function(log_target, n_iter, x0, proposal_sd) {
  if (n_iter < 2) {
    print("n_iter must be at least 2")
    return(NULL)
  } 
  
  if (proposal_sd <= 0) {
    print("proposal_sd must be positive.")
    return(NULL)
  }
  
  current_logp <- log_target(x0)
  if (is.na(current_logp) | current_logp == -Inf) {
    print("log_target(x0) is not finite: choose a different starting value.")
    return(NULL)
  }

  x <- numeric(n_iter)
  x[1] <- x0
  n_accept <- 0

  for (i in 2:n_iter) {
    proposal <- stats::rnorm(1, mean = x[i - 1], sd = proposal_sd)
    proposal_logp <- log_target(proposal)

    # Log acceptance ratio. No need for min(1, ratio): log(u) <= 0 always.
    log_alpha <- proposal_logp - current_logp

    # Proposals outside the support give log_alpha = -Inf and are rejected.
    if (log(stats::runif(1)) < log_alpha) {
      x[i] <- proposal
      current_logp <- proposal_logp
      n_accept <- n_accept + 1
    } else {
      x[i] <- x[i - 1]
    }
  }

  list(x = x, accept_rate = n_accept / (n_iter - 1))
}
