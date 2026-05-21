# Compile and fit the Stan model. Returns a list with the CmdStanMCMC fit
# object, the input data, and the quarterly date vector.

suppressPackageStartupMessages({
  library(dplyr)
  library(cmdstanr)
  library(posterior)
})

# Heuristic prior on initial state mean using HP-filtered values at t=2.
# Lambda=1600 is the standard quarterly value (Hodrick-Prescott 1997).
.hp_filter <- function(x, lambda = 1600) {
  n <- length(x)
  if (n < 4) return(x)
  D <- diag(n)
  D2 <- diff(diff(D))
  trend <- solve(diag(n) + lambda * crossprod(D2), x)
  as.numeric(trend)
}

# Estimation sample: 1966Q3 onward (matches the EViews ssest sample). This is
# the earliest date with non-missing rgdppc, unr, prt in the bundled data.
ESTIMATION_START <- as.Date("1966-09-01")
BREAK_DATE       <- as.Date("1983-12-01")

prepare_stan_data <- function(dataset) {
  est <- dataset |> dplyr::filter(date >= ESTIMATION_START)
  stopifnot(all(!is.na(est$rgdppc)), all(!is.na(est$unr)), all(!is.na(est$prt)))

  # HP-filter the three input series to get initial-state moments.
  y_hp <- .hp_filter(est$lrgdppc)
  u_hp <- .hp_filter(est$unr)
  p_hp <- .hp_filter(est$prt)
  cycle_ini <- est$lrgdppc - y_hp

  m0 <- c(
    y_hp[2],   # ystar_0
    u_hp[2],   # unrstar_0
    p_hp[2],   # prtstar_0
    0,         # cycle_0
    0          # cyclelag_0
  )
  # Prior variances loosely informed by Evans/Moore/Rees' EViews priors.
  P0 <- diag(c(2.0, 0.5, 0.5, 0.25, 0.25))

  break_idx <- max(which(est$date <= BREAK_DATE))

  list(
    stan_data = list(
      T = nrow(est),
      y = est$lrgdppc,
      u = est$unr,
      p = est$prt,
      break_idx = break_idx,
      m0 = m0,
      P0 = P0,
      prior_delta_mean = mean(diff(y_hp)),
      prior_delta_sd = 0.5
    ),
    dates = est$date,
    est = est
  )
}

fit_model <- function(stan_inputs,
                      stan_file = "stan/lf_uc.stan",
                      chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                      seed = 20260521,
                      parallel_chains = max(1, parallel::detectCores() - 1)) {
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = stan_inputs$stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 200,
    adapt_delta = 0.95,
    max_treedepth = 12
  )
  list(fit = fit, dates = stan_inputs$dates, est = stan_inputs$est)
}
