# Splice the bundled historical (pre-1978) series with the latest ABS data and
# return a single quarterly tibble ready for the Stan model.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readxl)
  library(lubridate)
})

# Map a date to its quarter using the ABS / EViews convention where the
# quarter is labelled by its end month (Mar/Jun/Sep/Dec, first day).
.qtr_end <- function(d) {
  m <- ((lubridate::month(d) - 1) %/% 3) * 3 + 3
  lubridate::make_date(lubridate::year(d), m, 1)
}

# Aggregate monthly LF data to quarterly averages and align dates to the first
# day of the quarter (consistent with the bundled historical series).
aggregate_lf_to_quarterly <- function(monthly) {
  monthly |>
    dplyr::mutate(qtr = .qtr_end(date)) |>
    dplyr::group_by(qtr) |>
    dplyr::summarise(
      unr    = if (sum(!is.na(unr))    == 3) mean(unr)    else NA_real_,
      prt    = if (sum(!is.na(prt))    == 3) mean(prt)    else NA_real_,
      civpop = if (sum(!is.na(civpop)) == 3) mean(civpop) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::rename(date = qtr) |>
    dplyr::arrange(date)
}

# Read the bundled historical (1964Q1-2019Q2) series, used to splice in the
# pre-1978 period that the ABS LF table does not cover.
read_historical_bundle <- function(path = "data/lf_hist.xlsx") {
  raw <- readxl::read_excel(path, sheet = "eviews", col_names = TRUE,
                            .name_repair = "minimal")
  names(raw)[1] <- "date"
  raw |>
    dplyr::select(date, rgdppc, unr, prt) |>
    dplyr::mutate(date = as.Date(date)) |>
    dplyr::filter(!is.na(date))
}

# Combine: history before 1978Q1 + freshly fetched ABS data thereafter.
# Returns a tibble(date, rgdppc, lrgdppc, unr, prt) covering 1964Q1 to the
# latest available quarter. Assumes the ABS reader functions are already in
# scope (sourced by the caller -- the run.R orchestrator does this).
build_dataset <- function(hist_path = "data/lf_hist.xlsx",
                          lf_xlsx, gdp_xlsx,
                          splice_date = as.Date("1978-06-01")) {
  hist <- read_historical_bundle(hist_path)

  lf_monthly <- read_lf_table01(lf_xlsx)
  lf_qtr <- aggregate_lf_to_quarterly(lf_monthly)
  gdp <- read_gdp_table01(gdp_xlsx) |>
    dplyr::mutate(date = .qtr_end(date))

  modern <- lf_qtr |>
    dplyr::inner_join(gdp, by = "date") |>
    dplyr::filter(!is.na(unr), !is.na(prt), !is.na(civpop), !is.na(gdp_cvm_sa)) |>
    dplyr::mutate(rgdppc_raw = gdp_cvm_sa / civpop)

  # Match the level of the bundled rgdppc series at the splice quarter so the
  # combined series has no level break. The bundled rgdppc is defined as
  # GDP($m) / civpop('000) up to rebasing differences in CVM and pop vintages.
  hist_pre <- hist |> dplyr::filter(date < splice_date)
  hist_at_splice <- hist |> dplyr::filter(date == splice_date) |> dplyr::pull(rgdppc)
  modern_at_splice <- modern |> dplyr::filter(date == splice_date) |> dplyr::pull(rgdppc_raw)

  if (length(hist_at_splice) == 1 && length(modern_at_splice) == 1 &&
      is.finite(hist_at_splice) && is.finite(modern_at_splice) && modern_at_splice > 0) {
    scale <- hist_at_splice / modern_at_splice
  } else {
    scale <- 1
  }
  modern <- modern |>
    dplyr::mutate(rgdppc = rgdppc_raw * scale) |>
    dplyr::select(date, rgdppc, unr, prt)

  combined <- dplyr::bind_rows(hist_pre, modern) |>
    dplyr::arrange(date) |>
    dplyr::distinct(date, .keep_all = TRUE)

  combined |>
    dplyr::filter(!is.na(rgdppc), !is.na(unr), !is.na(prt)) |>
    dplyr::mutate(lrgdppc = log(rgdppc) * 100)
}
