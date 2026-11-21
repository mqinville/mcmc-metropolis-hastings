#' Cleaned PWQMN 2025 total phosphorus measurements
#'
#' Total phosphorus samples from the Ontario Provincial (Stream) Water
#' Quality Monitoring Network, 2025, after cleaning with [clean_tp()]:
#' duplicate rows removed and the 20 below-detection results replaced by
#' half the detection limit. Built by `data-raw/tp_clean.R`.
#'
#' @format A tibble with 3,690 rows (449 sites) and 11 columns:
#' \describe{
#'   \item{site}{PWQMN station ID}
#'   \item{date}{collection date}
#'   \item{time}{collection time (character, HH:MM:SS)}
#'   \item{month}{collection month (factor, Jan-Dec)}
#'   \item{value}{total phosphorus concentration, µg/L}
#'   \item{censored}{`TRUE` if the result was below the detection limit}
#'   \item{det_limit}{laboratory detection limit, µg/L}
#'   \item{units}{measurement units as given in the source file}
#'   \item{remark}{laboratory remark code(s), `NA` if none}
#'   \item{remark_flag}{`TRUE` if a remark code is present}
#'   \item{lab_sample}{laboratory sample ID}
#' }
#' @source Ontario Ministry of the Environment, Conservation and Parks,
#'   PWQMN 2025 data:
#'   <https://files.ontario.ca/moe_mapping/downloads/2Water/PWQMN/PWQMN_2025.csv>
"tp_clean"
