## Builds the `tp_clean` dataset shipped with the package.
## Run from the package root:  source("data-raw/tp_clean.R")
## Raw file: https://files.ontario.ca/moe_mapping/downloads/2Water/PWQMN/PWQMN_2025.csv

devtools::load_all()

raw      <- read_pwqmn("data-raw/PWQMN_2025.csv")
cleaned  <- clean_tp(raw)
tp_clean <- cleaned$data
print(cleaned$log)

usethis::use_data(tp_clean, overwrite = TRUE)
