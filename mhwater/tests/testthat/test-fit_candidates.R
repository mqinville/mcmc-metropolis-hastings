test_that("fit_candidates recovers lognormal parameters", {
  set.seed(42)
  x <- rlnorm(5000, meanlog = 2.6, sdlog = 1.4)
  fits <- fit_candidates(x)
  expect_named(fits, c("lnorm", "gamma", "weibull"))
  expect_equal(unname(fits$lnorm$estimate), c(2.6, 1.4), tolerance = 0.05)
  expect_error(fit_candidates(c(1, 0, 2)), "positive")
})

test_that("compare_fits ranks the true family first", {
  set.seed(42)
  cmp <- compare_fits(fit_candidates(rlnorm(2000, 2.6, 1.4)))
  expect_equal(cmp$distribution[1], "lnorm")
  expect_equal(cmp$delta_AIC[1], 0)
  expect_equal(utils::tail(cmp$distribution, 1), "empirical")
})

test_that("exceed_prob and make_target match closed-form values", {
  set.seed(1)
  fit <- fitdistrplus::fitdist(rlnorm(500, 2.6, 1.4), "lnorm")
  p <- plnorm(30, fit$estimate[1], fit$estimate[2], lower.tail = FALSE)
  expect_equal(exceed_prob(fit, 30), unname(p))
  tgt <- make_target(fit, 30)
  expect_equal(tgt$p_exceed, unname(p))
  expect_equal(tgt$density(10),
               dlnorm(10, fit$estimate[1], fit$estimate[2]))
  expect_equal(tgt$log_density(10), log(tgt$density(10)))
})
