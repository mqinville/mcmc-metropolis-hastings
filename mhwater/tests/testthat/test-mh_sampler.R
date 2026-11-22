# The pre-generalisation sampler (hard-coded normal random walk), kept as a
# reference for the backwards-equivalence test.
old_mh_sampler <- function(log_target, n_iter, x0, proposal_sd) {
  current_logp <- log_target(x0)
  x <- numeric(n_iter)
  x[1] <- x0
  n_accept <- 0
  for (i in 2:n_iter) {
    proposal <- stats::rnorm(1, mean = x[i - 1], sd = proposal_sd)
    proposal_logp <- log_target(proposal)
    log_alpha <- proposal_logp - current_logp
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

test_that("normal_rw reproduces the old hard-coded sampler exactly", {
  log_target <- function(x) dlnorm(x, 2.656, 1.384, log = TRUE)
  set.seed(7)
  old <- old_mh_sampler(log_target, 5000, 10, 20)
  set.seed(7)
  new <- mh_sampler(log_target, 5000, 10, normal_rw(20))
  expect_identical(new$x, old$x)
  expect_identical(new$accept_rate, old$accept_rate)
  expect_identical(new$proposal, "normal_rw(sd = 20)")
})

test_that("mh_sampler recovers a standard normal", {
  set.seed(1)
  chain <- mh_sampler(function(x) dnorm(x, log = TRUE), 20000, 0,
                      normal_rw(2.4))
  x <- chain$x[-(1:1000)]
  expect_length(chain$x, 20000)
  expect_equal(chain$x[1], 0)
  expect_lt(abs(mean(x)), 0.1)
  expect_lt(abs(sd(x) - 1), 0.1)
  expect_gt(chain$accept_rate, 0.3)
  expect_lt(chain$accept_rate, 0.6)
})

test_that("every proposal recovers the Gamma(3, 1) mean", {
  log_target <- function(x) dgamma(x, shape = 3, rate = 1, log = TRUE)
  proposals <- list(normal_rw(2.5), lognormal_rw(0.8),
                    gamma_independence(2.5, 0.8))
  for (p in proposals) {
    set.seed(1)
    chain <- mh_sampler(log_target, 20000, 3, p)
    expect_lt(abs(mean(chain$x[-(1:1000)]) - 3), 0.15)
    expect_true(all(chain$x > 0))
  }
})

test_that("every proposal recovers the Lognormal(0, 0.5) mean", {
  log_target <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5, log = TRUE)
  proposals <- list(normal_rw(1), lognormal_rw(0.8),
                    gamma_independence(3, 2.5))
  for (p in proposals) {
    set.seed(1)
    chain <- mh_sampler(log_target, 20000, 1, p)
    expect_lt(abs(mean(chain$x[-(1:1000)]) - exp(0.125)), 0.1)
  }
})

test_that("mh_sampler reproduces the PWQMN lognormal exceedance probability", {
  set.seed(4000)
  tgt <- make_target(fitdistrplus::fitdist(tp_clean$value, "lnorm"),
                     threshold = 30)
  chain <- mh_sampler(tgt$log_density, 50000, 10, normal_rw(20))
  expect_lt(abs(mean(chain$x[-(1:1000)] > 30) - tgt$p_exceed), 0.05)
  expect_true(all(chain$x > 0))
})

test_that("step size drives the acceptance rate", {
  log_target <- function(x) dlnorm(x, 2.656, 1.384, log = TRUE)
  set.seed(1)
  expect_gt(mh_sampler(log_target, 5000, 10, normal_rw(0.01))$accept_rate,
            0.9)
  set.seed(1)
  expect_lt(mh_sampler(log_target, 5000, 10, normal_rw(5000))$accept_rate,
            0.05)
})

test_that("NA/NaN target values are rejected instead of crashing", {
  # NaN on the negative half-line, NA on (0, 0.1): neither may be accepted.
  log_target <- function(x) {
    if (x < 0) return(NaN)
    if (x < 0.1) return(NA_real_)
    dexp(x, log = TRUE)
  }
  set.seed(1)
  chain <- mh_sampler(log_target, 5000, 1, normal_rw(1))
  expect_length(chain$x, 5000)
  expect_true(all(chain$x >= 0.1))
  expect_gt(chain$accept_rate, 0)
})

test_that("mh_sampler rejects invalid inputs", {
  dn <- function(x) dnorm(x, log = TRUE)
  expect_error(mh_sampler(dn, 1, 0, normal_rw(1)), "at least 2")
  expect_error(mh_sampler(dn, c(10, 20), 0, normal_rw(1)), "single whole")
  expect_error(mh_sampler(function(x) dlnorm(x, log = TRUE), 100, -1,
                          normal_rw(1)), "not finite")
  expect_error(mh_sampler(function(x) NaN, 100, 0, normal_rw(1)),
               "not finite")
  expect_error(mh_sampler(dn, 100, 0, 1), "proposal must be a list")
  expect_error(mh_sampler(dn, 100, 0, list(draw = function(x) x)),
               "proposal must be a list")
})

test_that("proposal constructors validate their parameters", {
  expect_error(normal_rw(0), "sd must be a single positive number")
  expect_error(normal_rw(-1), "sd must be")
  expect_error(normal_rw(c(1, 2)), "sd must be")
  expect_error(lognormal_rw(NA), "sd must be")
  expect_error(gamma_independence(0, 1), "shape must be")
  expect_error(gamma_independence(1, -2), "rate must be")
  expect_error(lognormal_rw(1)$draw(0), "positive current state")
})

test_that("proposals return the documented structure", {
  p <- gamma_independence(2, 0.5)
  expect_s3_class(p, "mh_proposal")
  expect_identical(p$name, "gamma_independence(shape = 2, rate = 0.5)")
  # Independence sampler: log_q ignores the current state.
  expect_identical(p$log_q(3, 1), p$log_q(3, 100))
  # Lognormal random walk is asymmetric.
  q <- lognormal_rw(0.5)
  expect_false(isTRUE(all.equal(q$log_q(2, 1), q$log_q(1, 2))))
})
