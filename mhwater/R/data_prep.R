#' Read the raw PWQMN water-quality CSV
#'
#' Reads the Ontario Provincial (Stream) Water Quality Monitoring Network
#' (PWQMN) yearly CSV file. The file is Latin-1 encoded (the "µ" in µg/L
#' breaks under UTF-8), and every column is read as character so that
#' below-detection results such as `"<0.2"` are not lost.
#'
#' @param path path to the PWQMN CSV file (e.g. `"PWQMN_2025.csv"`).
#'
#' @returns A tibble with one row per measurement and all columns as character.
#' @importFrom readr read_csv locale cols col_character
#' @examples
#' \dontrun{
#' raw <- read_pwqmn("data/PWQMN_2025.csv")
#' }
#' @export
read_pwqmn <- function(path) {
  readr::read_csv(path,
                  locale = readr::locale(encoding = "latin1"),
                  col_types = readr::cols(.default = readr::col_character()))
}

#' Clean PWQMN total phosphorus measurements
#'
#' Extracts one analyte (total phosphorus by default) from the raw PWQMN data
#' and cleans it:
#'
#' 1. keeps only rows for `analyte`;
#' 2. strips the leading apostrophe from site IDs;
#' 3. removes duplicated rows (same site, date, time and result);
#' 4. parses below-detection-limit (BDL) results such as `"<0.2"` and
#'    replaces them with `bdl_factor * detection limit` (default DL/2);
#' 5. flags samples carrying a laboratory remark code.
#'
#' @param raw a data frame as returned by [read_pwqmn()].
#' @param analyte the `Analyte` value to keep. Default `"Phosphorus; total"`.
#' @param bdl_factor multiplier applied to the detection limit for BDL
#'   results. Default `0.5` (DL/2 substitution).
#'
#' @returns A list with two elements:
#'   \describe{
#'     \item{data}{a tibble with columns `site`, `date`, `time`, `month`,
#'       `value` (numeric result), `censored` (logical, `TRUE` if BDL),
#'       `det_limit`, `units`, `remark`, `remark_flag` and `lab_sample`.}
#'     \item{log}{a tibble of row counts after each cleaning step.}
#'   }
#' @importFrom dplyr filter rename mutate distinct select if_else n_distinct
#' @importFrom tibble tibble
#' @examples
#' raw <- data.frame(
#'   `Collection Site` = c("'001", "'001", "'002"),
#'   Analyte = "Phosphorus; total",
#'   `Collection Date` = c("2025-05-01", "2025-05-01", "2025-06-01"),
#'   `Collection Time` = "10:00:00",
#'   `Lab Sample ID` = c("1", "1", "2"),
#'   Result = c("12.5", "12.5", "<0.2"),
#'   Units = "ug/L",
#'   `Result Call` = c("Detected", "Detected", "BDL"),
#'   `Detection Limit` = "0.2",
#'   `Remark Code` = NA,
#'   check.names = FALSE
#' )
#' clean_tp(raw)
#' @export
clean_tp <- function(raw, analyte = "Phosphorus; total", bdl_factor = 0.5) {
  step <- function(df, label) tibble::tibble(step = label, n_rows = nrow(df))
  log <- step(raw, "Raw file (all analytes)")

  d <- raw |>
    dplyr::filter(.data$Analyte == analyte) |>
    dplyr::rename(site = "Collection Site", date = "Collection Date",
                  time = "Collection Time", result_raw = "Result",
                  result_call = "Result Call", det_limit = "Detection Limit",
                  lab_sample = "Lab Sample ID", remark = "Remark Code",
                  units = "Units") |>
    dplyr::mutate(site = sub("^'", "", .data$site))
  if (nrow(d) == 0) stop("No rows found for analyte '", analyte, "'.")
  log <- rbind(log, step(d, paste0("Rows for '", analyte, "'")))

  if (length(unique(d$units)) != 1)
    warning("More than one unit found for this analyte: ",
            paste(unique(d$units), collapse = ", "))

  # Exact repeated rows (same site/date/time/result) are data-entry duplicates
  d <- dplyr::distinct(d, .data$site, .data$date, .data$time,
                       .data$result_raw, .keep_all = TRUE)
  log <- rbind(log, step(d, "After removing duplicate rows"))
  if (any(duplicated(d[, c("site", "date", "time")])))
    warning("Some samples have two different results at the same site/date/time.")

  d <- d |>
    dplyr::mutate(
      censored    = .data$result_call == "BDL" | grepl("^<", .data$result_raw),
      det_limit   = suppressWarnings(as.numeric(.data$det_limit)),
      result_num  = suppressWarnings(as.numeric(.data$result_raw)),
      value       = dplyr::if_else(.data$censored,
                                   bdl_factor * .data$det_limit,
                                   .data$result_num),
      date        = as.Date(.data$date),
      month       = factor(format(.data$date, "%b"), levels = month.abb),
      remark_flag = !is.na(.data$remark)
    ) |>
    dplyr::filter(!is.na(.data$value), .data$value > 0) |>
    dplyr::select("site", "date", "time", "month", "value", "censored",
                  "det_limit", "units", "remark", "remark_flag", "lab_sample")
  log <- rbind(log, step(d, "Final clean samples"))

  list(data = d, log = log)
}

#' Descriptive summary of a positive measurement vector
#'
#' @param x numeric vector of positive measurements.
#' @param threshold guideline value used for the empirical exceedance
#'   proportion. Default 30 (the Ontario PWQO for total phosphorus in
#'   rivers and streams, in µg/L).
#'
#' @returns A one-row tibble: n, mean, sd, median, 90th and 99th percentiles,
#'   max, skewness and the empirical proportion exceeding `threshold`.
#' @importFrom stats sd median quantile
#' @examples
#' describe_tp(rlnorm(100, 2.6, 1.4))
#' @export
describe_tp <- function(x, threshold = 30) {
  tibble::tibble(
    n = length(x), mean = mean(x), sd = stats::sd(x),
    median = stats::median(x),
    q90 = unname(stats::quantile(x, 0.90)),
    q99 = unname(stats::quantile(x, 0.99)),
    max = max(x),
    skewness = mean((x - mean(x))^3) / stats::sd(x)^3,
    emp_P_exceed = mean(x > threshold)
  )
}
