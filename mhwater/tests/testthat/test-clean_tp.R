make_raw <- function() {
  data.frame(
    `Collection Site` = c("'001", "'001", "'002", "'003", "'004"),
    Analyte = c(rep("Phosphorus; total", 4), "Nitrate"),
    `Collection Date` = c("2025-05-01", "2025-05-01", "2025-06-01",
                          "2025-07-01", "2025-07-01"),
    `Collection Time` = "10:00:00",
    `Lab Sample ID` = c("1", "1", "2", "3", "4"),
    Result = c("12.5", "12.5", "<0.2", "40", "5"),
    Units = "ug/L",
    `Result Call` = c("Detected", "Detected", "BDL", "Detected", "Detected"),
    `Detection Limit` = "0.2",
    `Remark Code` = c(NA, NA, NA, "UHT", NA),
    check.names = FALSE
  )
}

test_that("clean_tp keeps the analyte, drops duplicates and parses BDL", {
  out <- clean_tp(make_raw())
  d <- out$data
  expect_equal(nrow(d), 3)                       # nitrate + duplicate removed
  expect_false(any(grepl("^'", d$site)))         # apostrophe stripped
  expect_equal(d$value[d$censored], 0.1)         # "<0.2" -> DL/2
  expect_equal(sum(d$remark_flag), 1)
  expect_equal(out$log$n_rows, c(5, 4, 3, 3))
})

test_that("clean_tp respects bdl_factor and errors on a missing analyte", {
  d <- clean_tp(make_raw(), bdl_factor = 1)$data
  expect_equal(d$value[d$censored], 0.2)
  expect_error(clean_tp(make_raw(), analyte = "Gold"), "No rows")
})

test_that("shipped tp_clean data has the expected shape", {
  expect_equal(nrow(tp_clean), 3690)
  expect_true(all(tp_clean$value > 0))
  expect_equal(sum(tp_clean$censored), 20)
})
