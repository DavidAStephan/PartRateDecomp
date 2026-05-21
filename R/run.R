# End-to-end pipeline. Run from the repo root:
#   Rscript R/run.R
#
# Steps:
#   1. Download the latest ABS Labour Force (Table 1) and National Accounts
#      (Table 1) workbooks.
#   2. Splice with the bundled historical series (data/lf_hist.xlsx) covering
#      1964Q1-1977Q4.
#   3. Compile and fit the Stan model defined in stan/lf_uc.stan.
#   4. Write outputs/{states,parameters}.csv and the chart in PNG + PDF.
#   5. Refresh docs/index.html with the latest run timestamp.
#
# Environment knobs (for the GitHub Actions runner):
#   PRD_CACHE_DIR  -- where to put the downloaded xlsx files (default: data/raw/)
#   PRD_CHAINS     -- number of MCMC chains (default: 4)
#   PRD_WARMUP     -- warmup iterations (default: 1000)
#   PRD_SAMPLING   -- post-warmup iterations (default: 1000)

suppressPackageStartupMessages({
  library(dplyr)
})

source("R/abs_fetch.R")
source("R/data_prep.R")
source("R/fit.R")
source("R/outputs.R")

main <- function() {
  cache_dir <- Sys.getenv("PRD_CACHE_DIR", "data/raw")
  chains    <- as.integer(Sys.getenv("PRD_CHAINS", "4"))
  warmup    <- as.integer(Sys.getenv("PRD_WARMUP", "1000"))
  sampling  <- as.integer(Sys.getenv("PRD_SAMPLING", "1000"))

  lf_xlsx  <- file.path(cache_dir, "lf_table01.xlsx")
  gdp_xlsx <- file.path(cache_dir, "gdp_table01.xlsx")
  if (!file.exists(lf_xlsx) || !file.exists(gdp_xlsx) ||
      isTRUE(as.logical(Sys.getenv("PRD_FORCE_FETCH", "FALSE")))) {
    message("[1/5] Fetching ABS data...")
    lf_xlsx  <- fetch_lf_table01(cache_dir)
    gdp_xlsx <- fetch_gdp_table01(cache_dir)
  } else {
    message("[1/5] Using cached ABS data in ", cache_dir)
  }

  message("[2/5] Building combined dataset...")
  dataset <- build_dataset(lf_xlsx = lf_xlsx, gdp_xlsx = gdp_xlsx)
  message(sprintf("    Sample: %s to %s (%d quarters)",
                  format(min(dataset$date)), format(max(dataset$date)),
                  nrow(dataset)))

  message("[3/5] Preparing Stan inputs...")
  stan_inputs <- prepare_stan_data(dataset)

  message("[4/5] Fitting model (Stan)...")
  fit_obj <- fit_model(stan_inputs,
                       chains = chains, iter_warmup = warmup, iter_sampling = sampling)
  fit_obj$fit$cmdstan_diagnose()

  message("[5/5] Writing outputs...")
  write_outputs(fit_obj)

  # Update docs page with the run timestamp + last data quarter.
  refresh_docs(last_quarter = max(stan_inputs$est$date))
  message("Done. See outputs/ and docs/.")
}

refresh_docs <- function(last_quarter, docs_dir = "docs") {
  dir.create(docs_dir, showWarnings = FALSE, recursive = TRUE)
  # Copy chart so docs/ is self-contained for GitHub Pages.
  file.copy("outputs/participation_rate_chart.png",
            file.path(docs_dir, "participation_rate_chart.png"),
            overwrite = TRUE)
  file.copy("outputs/participation_rate_chart.pdf",
            file.path(docs_dir, "participation_rate_chart.pdf"),
            overwrite = TRUE)
  file.copy("outputs/states.csv", file.path(docs_dir, "states.csv"),
            overwrite = TRUE)

  index_path <- file.path(docs_dir, "index.html")
  tmpl <- readLines("docs/index.template.html", warn = FALSE)
  tmpl <- gsub("{{LAST_QUARTER}}", format(last_quarter, "%b %Y"), tmpl, fixed = TRUE)
  tmpl <- gsub("{{RUN_TIMESTAMP}}",
               format(Sys.time(), "%Y-%m-%d %H:%M %Z"), tmpl, fixed = TRUE)
  writeLines(tmpl, index_path)
}

if (sys.nframe() == 0) {
  main()
}
