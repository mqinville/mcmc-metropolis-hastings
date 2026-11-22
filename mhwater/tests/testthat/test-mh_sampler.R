# Reference chains from the log-scale sampler (commit 3d6ffc2), saved before
# the switch to the raw density scale:
#   set.seed(123); mh_sampler(function(x) dnorm(x, log = TRUE), 2000, 0,
#                             normal_rw(2.4))            # with log_q proposals
#   set.seed(123); mh_sampler(function(x) dlnorm(x, 0, 0.5, log = TRUE), 2000,
#                             1, lognormal_rw(0.8))
test_that("raw-scale sampler reproduces the log-scale reference chains", {
  ref <- readRDS(test_path("fixtures", "log_scale_chains.rds"))

  set.seed(123)
  norm <- mh_sampler(function(x) dnorm(x), 2000, 0, normal_rw(2.4))
  expect_equal(norm, ref$norm)

  set.seed(123)
  lnorm <- mh_sampler(function(x) dlnorm(x, 0, 0.5), 2000, 1,
                      lognormal_rw(0.8))
  expect_equal(lnorm, ref$lnorm)
})

test_that("mh_sampler recovers a standard normal", {
  set.seed(1)
  chain <- mh_sampler(function(x) dnorm(x), 20000, 0, normal_rw(2.4))
  x <- chain$x[-(1:1000)]
  expect_length(chain$x, 20000)
  expect_equal(chain$x[1], 0)
  expect_lt(abs(mean(x)), 0.1)
  expect_lt(abs(sd(x) - 1), 0.1)
  expect_gt(chain$accept_rate, 0.3)
  expect_lt(chain$accept_rate, 0.6)
})

test_that("every proposal recovers the Gamma(3, 1) mean", {
  target <- function(x) dgamma(x, shape = 3, rate = 1)
  proposals <- list(normal_rw(2.5), lognormal_rw(0.8),
                    gamma_independence(2.5, 0.8))
  for (p in proposals) {
    set.seed(1)
    chain <- mh_sampler(target, 20000, 3, p)
    expect_lt(abs(mean(chain$x[-(1:1000)]) - 3), 0.15)
    expect_true(all(chain$x > 0))
  }
})

test_that("every proposal recovers the Lognormal(0, 0.5) mean", {
  target <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5)
  proposals <- list(normal_rw(1), lognormal_rw(0.8),
                    gamma_independence(3, 2.5))
  for (p in proposals) {
    set.seed(1)
    chain <- mh_sampler(target, 20000, 1, p)
    expect_lt(abs(mean(chain$x[-(1:1000)]) - exp(0.125)), 0.1)
  }
})

test_that("mh_sampler reproduces the PWQMN lognormal exceedance probability", {
  set.seed(4000)
  tgt <- make_target(fitdistrplus::fitdist(tp_clean$value, "lnorm"),
                     threshold = 30)
  chain <- mh_sampler(tgt$density, 50000, 10, normal_rw(20))
  expect_lt(abs(mean(chain$x[-(1:1000)] > 30) - tgt$p_exceed), 0.05)
  expect_true(all(chain$x > 0))
})

test_that("step size drives the acceptance rate", {
  target <- function(x) dlnorm(x, 2.656, 1.384)
  set.seed(1)
  expect_gt(mh_sampler(target, 5000, 10, normal_rw(0.01))$accept_rate, 0.9)
  set.seed(1)
  expect_lt(mh_sampler(target, 5000, 10, normal_rw(5000))$accept_rate, 0.05)
})

test_that("NA/NaN target values are rejected instead of crashing", {
  # NaN on the negative half-line, NA on (0, 0.1): neither may be accepted.
  target <- function(x) {
    if (x < 0) return(NaN)
    if (x < 0.1) return(NA_real_)
    dexp(x)
  }
  set.seed(1)
  chain <- mh_sampler(target, 5000, 1, normal_rw(1))
  expect_length(chain$x, 5000)
  expect_true(all(chain$x >= 0.1))
  expect_gt(chain$accept_rate, 0)
})

test_that("mh_sampler prints a message and returns NULL on invalid inputs", {
  # Invalid input prints `message` and returns NULL instead of a chain.
  expect_invalid <- function(call, message) {
    expect_output(result <- call, message)
    expect_null(result)
  }
  dn <- function(x) dnorm(x)
  expect_invalid(mh_sampler(dn, 1, 0, normal_rw(1)), "at least 2")
  expect_invalid(mh_sampler(dn, c(10, 20), 0, normal_rw(1)), "single whole")
  expect_invalid(mh_sampler(function(x) dlnorm(x), 100, -1, normal_rw(1)),
                 "target\\(x0\\)")
  expect_invalid(mh_sampler(dn, 100, 1e4, normal_rw(1)), "target\\(x0\\)")
  expect_invalid(mh_sampler(function(x) NaN, 100, 0, normal_rw(1)),
                 "target\\(x0\\)")
  expect_invalid(mh_sampler(dn, 100, 0, 1), "proposal must be a list")
  expect_invalid(mh_sampler(dn, 100, 0, list(draw = function(x) x)),
                 "proposal must be a list")
  # Old log_q-style proposals are no longer accepted.
  expect_invalid(mh_sampler(dn, 100, 0, list(draw = function(x) x,
                                             log_q = function(to, from) 0)),
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
  # Independence sampler: q ignores the current state.
  expect_identical(p$q(3, 1), p$q(3, 100))
  expect_equal(p$q(3, 1), dgamma(3, 2, 0.5))
  # Lognormal random walk is asymmetric.
  q <- lognormal_rw(0.5)
  expect_false(isTRUE(all.equal(q$q(2, 1), q$q(1, 2))))
})
