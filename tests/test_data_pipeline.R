# Run with: Rscript -e 'testthat::test_dir("tests")'

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(tibble)
})

# testthat may run tests with a different working directory; resolve project
# root by walking up until we find DESCRIPTION.
.find_root <- function(start = getwd()) {
  p <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "DESCRIPTION"))) return(p)
    parent <- dirname(p)
    if (parent == p) stop("Could not find project root")
    p <- parent
  }
}
ROOT <- .find_root()
source(file.path(ROOT, "R/abs_fetch.R"))
source(file.path(ROOT, "R/data_prep.R"))
abs_path <- function(rel) file.path(ROOT, rel)

test_that("aggregate_lf_to_quarterly averages complete quarters and NA-fills partial ones", {
  monthly <- tibble::tibble(
    date = seq(as.Date("2024-01-01"), as.Date("2024-08-01"), by = "month"),
    unr  = c(4.0, 4.2, 4.4,  4.1, 4.3, 4.5,  4.6, 4.7),
    prt  = c(66.0, 66.1, 66.2,  66.3, 66.4, 66.5,  66.6, 66.7),
    civpop = rep(21000, 8)
  )
  out <- aggregate_lf_to_quarterly(monthly)
  # Full quarters (Q1, Q2) and a partial Q3 placeholder.
  expect_equal(out$date, as.Date(c("2024-03-01", "2024-06-01", "2024-09-01")))
  expect_equal(out$unr, c(mean(c(4.0, 4.2, 4.4)), mean(c(4.1, 4.3, 4.5)), NA_real_))
  expect_equal(out$prt, c(mean(c(66.0, 66.1, 66.2)), mean(c(66.3, 66.4, 66.5)), NA_real_))
})

test_that("read_lf_table01 extracts the three series by Series ID", {
  skip_if_not(file.exists(abs_path("data/raw/lf_table01.xlsx")),
              "no cached ABS workbook (CI fetches it fresh)")
  out <- read_lf_table01(abs_path("data/raw/lf_table01.xlsx"))
  expect_true(all(c("date", "unr", "prt", "civpop") %in% names(out)))
  expect_true(min(out$date) <= as.Date("1978-03-01"))
  # The published participation rate lives in [50, 80] %.
  expect_true(all(out$prt > 50 & out$prt < 80, na.rm = TRUE))
  # Unemployment rate has been within [1, 12] %.
  expect_true(all(out$unr >= 1 & out$unr <= 12, na.rm = TRUE))
})

test_that("read_gdp_table01 returns monotone-ish quarterly GDP", {
  skip_if_not(file.exists(abs_path("data/raw/gdp_table01.xlsx")),
              "no cached ABS workbook (CI fetches it fresh)")
  out <- read_gdp_table01(abs_path("data/raw/gdp_table01.xlsx"))
  expect_true(all(c("date", "gdp_cvm_sa") %in% names(out)))
  expect_true(nrow(out) > 200)
  # GDP grew over the full sample.
  expect_gt(tail(out$gdp_cvm_sa, 1), head(out$gdp_cvm_sa, 1))
})

test_that("build_dataset produces a continuous quarterly series at the splice", {
  skip_if_not(file.exists(abs_path("data/raw/lf_table01.xlsx")) &&
              file.exists(abs_path("data/raw/gdp_table01.xlsx")),
              "no cached ABS workbooks")
  ds <- build_dataset(hist_path = abs_path("data/lf_hist.xlsx"),
                      lf_xlsx = abs_path("data/raw/lf_table01.xlsx"),
                      gdp_xlsx = abs_path("data/raw/gdp_table01.xlsx"))
  # Series complete (no NAs in the three key columns).
  expect_true(all(!is.na(ds$rgdppc)))
  expect_true(all(!is.na(ds$unr)))
  expect_true(all(!is.na(ds$prt)))
  # No date gaps (each successive date is exactly one quarter later).
  gaps <- as.integer(diff(ds$date))
  expect_true(all(gaps >= 89 & gaps <= 92))
  # No jump > 5% in log-GDP-per-capita between adjacent quarters at the splice.
  splice_idx <- which(ds$date == as.Date("1978-06-01"))
  if (length(splice_idx) == 1 && splice_idx > 1) {
    step <- abs(ds$lrgdppc[splice_idx] - ds$lrgdppc[splice_idx - 1])
    expect_lt(step, 5)
  }
})
