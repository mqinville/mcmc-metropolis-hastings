test_that("mh_sampler recovers a standard normal", {
  set.seed(1)
  chain <- mh_sampler(function(x) dnorm(x, log = TRUE), 20000, 0, 2.4)
  x <- chain$x[-(1:1000)]
  expect_length(chain$x, 20000)
  expect_equal(chain$x[1], 0)
  expect_lt(abs(mean(x)), 0.1)
  expect_lt(abs(sd(x) - 1), 0.1)
  expect_gt(chain$accept_rate, 0.3)
  expect_lt(chain$accept_rate, 0.6)
})

test_that("mh_sampler reproduces the PWQMN lognormal exceedance probability", {
  set.seed(4000)
  tgt <- make_target(fitdistrplus::fitdist(tp_clean$value, "lnorm"),
                     threshold = 30)
  chain <- mh_sampler(tgt$log_density, 50000, 10, 20)
  expect_lt(abs(mean(chain$x[-(1:1000)] > 30) - tgt$p_exceed), 0.05)
  expect_true(all(chain$x > 0))
})

test_that("step size drives the acceptance rate", {
  log_target <- function(x) dlnorm(x, 2.656, 1.384, log = TRUE)
  set.seed(1)
  expect_gt(mh_sampler(log_target, 5000, 10, 0.01)$accept_rate, 0.9)
  set.seed(1)
  expect_lt(mh_sampler(log_target, 5000, 10, 5000)$accept_rate, 0.05)
})

test_that("mh_sampler rejects invalid inputs", {
  expect_error(mh_sampler(function(x) dlnorm(x, log = TRUE), 100, -1, 1),
               "not finite")
  expect_error(mh_sampler(function(x) dnorm(x, log = TRUE), 100, 0, 0),
               "positive")
  expect_error(mh_sampler(function(x) dnorm(x, log = TRUE), 1, 0, 1),
               "at least 2")
})
