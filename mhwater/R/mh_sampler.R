# Metropolis-Hastings sampler (Objective 1) -- to be written by Felix.
#
# Suggested interface so it plugs straight into the water-quality target:
#
#   tgt <- make_target(fit_candidates(tp_clean$value)$lnorm)
#   chain <- mh_sampler(log_target = tgt$log_density, n_iter = 10000,
#                       x0 = 10, proposal_sd = 1)
#   mean(chain$x[-(1:burn)] > tgt$threshold)   # compare with tgt$p_exceed
#
# Template (uncomment, fill in, then run devtools::document()):
#
# #' Random-walk Metropolis-Hastings sampler
# #'
# #' @param log_target function returning the log target density at x.
# #' @param n_iter number of iterations.
# #' @param x0 starting value.
# #' @param proposal_sd standard deviation of the proposal.
# #' @returns A list with the chain `x` and the acceptance rate.
# #' @export
# mh_sampler <- function(log_target, n_iter, x0, proposal_sd) {
#
# }
