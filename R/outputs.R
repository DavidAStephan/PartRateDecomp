# Post-processing: build the published artefacts -- a CSV of state estimates,
# a CSV of posterior parameter summaries, and a 2-panel chart replicating the
# Evans/Moore/Rees figure (trend participation rate + cyclical participation
# rate with 95% credible bands).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(posterior)
})

.summarise_state <- function(draws_mat) {
  apply(draws_mat, 2, function(v) {
    c(mean = mean(v), sd = sd(v),
      q025 = unname(quantile(v, 0.025)),
      q975 = unname(quantile(v, 0.975)))
  }) |>
    t() |>
    as.data.frame() |>
    tibble::as_tibble()
}

build_state_table <- function(fit_obj) {
  draws <- posterior::as_draws_matrix(fit_obj$fit$draws(
    c("ystar_sm", "unrstar_sm", "prtstar_sm",
      "cycle_sm", "cyclelag_sm", "prt_cycle_sm")
  ))
  dates <- fit_obj$dates

  pick <- function(name) {
    cols <- grep(paste0("^", name, "\\["), colnames(draws), value = TRUE)
    cols <- cols[order(as.integer(gsub("[^0-9]", "", cols)))]
    draws[, cols, drop = FALSE]
  }

  parts <- list(
    ystar    = .summarise_state(pick("ystar_sm")),
    unrstar  = .summarise_state(pick("unrstar_sm")),
    prtstar  = .summarise_state(pick("prtstar_sm")),
    cycle    = .summarise_state(pick("cycle_sm")),
    prt_cyc  = .summarise_state(pick("prt_cycle_sm"))
  )

  bind_rename <- function(df, prefix) {
    df |> dplyr::rename_with(~ paste0(prefix, "_", .x))
  }

  tibble::tibble(date = dates) |>
    dplyr::bind_cols(bind_rename(parts$ystar,   "ystar")) |>
    dplyr::bind_cols(bind_rename(parts$unrstar, "unrstar")) |>
    dplyr::bind_cols(bind_rename(parts$prtstar, "prtstar")) |>
    dplyr::bind_cols(bind_rename(parts$cycle,   "cycle")) |>
    dplyr::bind_cols(bind_rename(parts$prt_cyc, "prt_cycle"))
}

build_parameter_table <- function(fit_obj) {
  pars <- c("delta", "phi1", "phi2", "pacf1", "pacf2",
            "kappa1", "kappa2", "theta1", "theta2",
            paste0("sigma[", rep(1:2, each = 4), ",", rep(1:4, 2), "]"),
            "rho[1]", "rho[2]", "nu")
  s <- fit_obj$fit$summary(variables = pars,
                           mean = mean, sd = sd,
                           q2.5 = ~ quantile(.x, 0.025),
                           q97.5 = ~ quantile(.x, 0.975),
                           rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  tibble::as_tibble(s)
}

make_chart <- function(state_tbl, est_data, png_path, pdf_path,
                       since = as.Date("1980-01-01"),
                       run_date = Sys.Date()) {
  st <- state_tbl |> dplyr::filter(date >= since)
  raw <- est_data |> dplyr::filter(date >= since)

  p_trend <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = st,
                         ggplot2::aes(date, ymin = prtstar_q025, ymax = prtstar_q975),
                         fill = "grey80") +
    ggplot2::geom_line(data = raw, ggplot2::aes(date, prt),
                       colour = "grey40", linewidth = 0.35) +
    ggplot2::geom_line(data = st, ggplot2::aes(date, prtstar_mean),
                       colour = "#1f4e79", linewidth = 0.9) +
    ggplot2::labs(title = "Trend participation rate",
                  subtitle = "Posterior mean with 95% credible interval; thin line = observed",
                  y = "Per cent", x = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  p_cycle <- ggplot2::ggplot(st, ggplot2::aes(date, prt_cycle_mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = prt_cycle_q025, ymax = prt_cycle_q975),
                         fill = "grey80") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_line(colour = "#c75146", linewidth = 0.9) +
    ggplot2::labs(title = "Cyclical participation rate",
                  subtitle = "Posterior mean with 95% credible interval",
                  y = "Percentage points", x = NULL,
                  caption = sprintf("Last updated %s. Source: ABS 6202.0 and 5206.0. Model: Evans/Moore/Rees (2018) replicated in Stan.",
                                    format(run_date, "%d %b %Y"))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  combined <- p_trend / p_cycle
  ggplot2::ggsave(png_path, combined, width = 9, height = 7, dpi = 150)
  ggplot2::ggsave(pdf_path, combined, width = 9, height = 7)
  invisible(list(png = png_path, pdf = pdf_path))
}

write_outputs <- function(fit_obj, out_dir = "outputs") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  state_tbl <- build_state_table(fit_obj)
  param_tbl <- build_parameter_table(fit_obj)
  readr::write_csv(state_tbl, file.path(out_dir, "states.csv"))
  readr::write_csv(param_tbl, file.path(out_dir, "parameters.csv"))
  make_chart(state_tbl, fit_obj$est,
             file.path(out_dir, "participation_rate_chart.png"),
             file.path(out_dir, "participation_rate_chart.pdf"))
  list(state_tbl = state_tbl, param_tbl = param_tbl)
}
