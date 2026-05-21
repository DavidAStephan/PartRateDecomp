# Download the latest ABS source workbooks and extract the series we need.
#
# Public functions:
#   fetch_lf_table01(dest_dir)   -> path to the saved xlsx
#   fetch_gdp_table01(dest_dir)  -> path to the saved xlsx
#   read_lf_table01(path)        -> tibble with monthly date, unr, prt, civpop
#   read_gdp_table01(path)       -> tibble with quarterly date, gdp_cvm_sa
#
# The ABS publishes Time Series Workbooks with a fixed shape: the "Data1" sheet
# has 10 header rows (one of which is the Series ID row), then dated rows of
# numeric data starting at row 11. We pull series by Series ID so the code is
# robust to ABS reordering or renaming columns.

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(tibble)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(httr)
})

ABS_LF_URL <- paste0(
  "https://www.abs.gov.au/statistics/labour/employment-and-unemployment/",
  "labour-force-australia/latest-release/62020001.xlsx"
)

ABS_GDP_URL <- paste0(
  "https://www.abs.gov.au/statistics/economy/national-accounts/",
  "australian-national-accounts-national-income-expenditure-and-product/",
  "latest-release/5206001_key_aggregates.xlsx"
)

# Series IDs we use (these are stable ABS identifiers).
SERIES <- list(
  unr_sa     = "A84423050A",  # Unemployment rate, Persons, SA, %
  prt_sa     = "A84423051C",  # Participation rate, Persons, SA, %
  civpop_15  = "A84423091W",  # Civilian population 15+, Persons, Original, '000
  gdp_cvm_sa = "A2304402X"    # GDP Chain volume measures, SA, $m
)

.download <- function(url, dest) {
  resp <- httr::GET(
    url,
    httr::user_agent("PartRateDecomp/1.0 (https://github.com/DavidAStephan/PartRateDecomp)"),
    httr::write_disk(dest, overwrite = TRUE),
    httr::timeout(120)
  )
  if (httr::status_code(resp) >= 400) {
    stop("Download failed: ", url, " (HTTP ", httr::status_code(resp), ")")
  }
  invisible(dest)
}

fetch_lf_table01 <- function(dest_dir) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(dest_dir, "lf_table01.xlsx")
  .download(ABS_LF_URL, out)
  out
}

fetch_gdp_table01 <- function(dest_dir) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(dest_dir, "gdp_table01.xlsx")
  .download(ABS_GDP_URL, out)
  out
}

# Extract a single ABS series by Series ID from a Time Series Workbook sheet.
# Header rows 1-10 of "Data1" describe each column; the Series ID is on row 10.
# Returns a tibble(date, value).
.read_abs_series <- function(path, sheet, series_id) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                            .name_repair = "minimal")
  # Series IDs live in row 10 (header offset).
  sid_row <- as.character(raw[10, ])
  col_idx <- which(sid_row == series_id)
  if (length(col_idx) == 0) {
    stop("Series ID '", series_id, "' not found in ", path, " / ", sheet)
  }
  if (length(col_idx) > 1) col_idx <- col_idx[1]

  # Date is column 1; data starts at row 11.
  dates_raw <- raw[[1]][11:nrow(raw)]
  vals_raw <- raw[[col_idx]][11:nrow(raw)]

  # When read with col_names=FALSE, readxl returns column 1 as character (since
  # row 1 contains a string). Excel stores dates as serial numbers from
  # 1899-12-30, so coerce numerically first; fall back to date parsing.
  dates_num <- suppressWarnings(as.numeric(dates_raw))
  dates <- suppressWarnings(as.Date(dates_num, origin = "1899-12-30"))
  if (all(is.na(dates))) {
    dates <- suppressWarnings(as.Date(as.POSIXct(dates_raw, tz = "UTC")))
  }
  vals <- suppressWarnings(as.numeric(vals_raw))

  tibble::tibble(date = dates, value = vals) |>
    dplyr::filter(!is.na(date))
}

read_lf_table01 <- function(path) {
  unr    <- .read_abs_series(path, "Data1", SERIES$unr_sa)    |> dplyr::rename(unr = value)
  prt    <- .read_abs_series(path, "Data1", SERIES$prt_sa)    |> dplyr::rename(prt = value)
  civpop <- .read_abs_series(path, "Data1", SERIES$civpop_15) |> dplyr::rename(civpop = value)
  unr |>
    dplyr::full_join(prt, by = "date") |>
    dplyr::full_join(civpop, by = "date") |>
    dplyr::arrange(date)
}

read_gdp_table01 <- function(path) {
  .read_abs_series(path, "Data1", SERIES$gdp_cvm_sa) |>
    dplyr::rename(gdp_cvm_sa = value) |>
    dplyr::filter(!is.na(gdp_cvm_sa))
}
