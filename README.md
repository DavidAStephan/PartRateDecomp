# PartRateDecomp

Bayesian unobserved-components decomposition of the Australian
labour-force participation rate. The model is implemented in
[Stan](https://mc-stan.org/) and refreshed automatically every time the
ABS releases new Labour Force data.

The methodology replicates **Evans, Moore & Rees (2018), *The Cyclical
Behaviour of Labour Force Participation*** ([RBA Bulletin, September
2018](https://www.rba.gov.au/publications/bulletin/2018/sep/the-cyclical-behaviour-of-labour-force-participation.html)),
translating their EViews state-space program into a Stan model fitted
by HMC.

## Latest chart

The most recent run is published as a chart on the project's
[GitHub Pages site](https://davidastephan.github.io/PartRateDecomp/).
Raw outputs (`states.csv`, `parameters.csv`, PNG + PDF) live in
[`outputs/`](outputs/).

## Model

State vector `x[t] = (ystar, unrstar, prtstar, cycle, cyclelag)`.

**Signal equations**

```
lrgdppc[t] = ystar[t] + cycle[t]
unr[t]    = unrstar[t] + kappa1 * cycle[t] + kappa2 * cyclelag[t]
prt[t]    = prtstar[t] + theta1 * cycle[t] + theta2 * cyclelag[t]
```

**State equations**

```
ystar[t]    = delta + ystar[t-1]                    + e1[t]
unrstar[t]  = unrstar[t-1]                          + e2[t]
prtstar[t]  = prtstar[t-1]                          + e3[t]
cycle[t]    = phi1 * cycle[t-1] + phi2 * cyclelag[t-1] + e4[t]
cyclelag[t] = cycle[t-1]
```

Innovations `(e1, e2, e3, e4)` are jointly Normal with `cov(e1, e4) =
sigma1 * sigma4 * rho`. Variances and `rho` follow a regime switch at
1983Q4 (the floating of the AUD), with separate parameters for the
pre- and post-break periods.

Compared with the original EViews fit (ML), the Stan implementation
adds weakly informative priors and recovers the full posterior, so
estimates ship with credible intervals rather than asymptotic standard
errors. It also makes three changes that materially affect the
decomposition on the extended sample:

1. **Hard stationarity on the cycle.** The AR(2) cycle is
   reparameterised via partial autocorrelations
   (Barndorff-Nielsen-Schou / Jones 1987), so the cycle is guaranteed
   stationary. On the extended sample the unconstrained EViews ML
   estimate sits at `phi1 + phi2 ~ 1`, which makes the cycle a
   near-random-walk and causes it to absorb the persistent post-2000
   drift in Australian output and participation. Forcing stationarity
   pushes that drift back into the random-walk trends where it
   belongs.
2. **Sign constraints break the cycle's sign-flip identification
   ambiguity.** `kappa1 <= 0` (Okun's law) and `theta1 >= 0`
   (procyclical participation) pin down the cycle's sign so all chains
   converge to the same mode.
3. **Student-t cycle innovations** via the standard
   scale-mixture-of-normals (`lambda[t] ~ inv_gamma(nu/2, nu/2)`)
   absorb one-off shocks (COVID-2020 in particular) without distorting
   the surrounding trend or cycle.

The smoothed states (`states.csv`) are well-identified by the data;
some parameter-level posteriors (cycle persistence and the post-1984
innovation variances) remain weakly identified -- a classic feature of
multivariate UC models -- so `parameters.csv` reports wider credible
intervals on those.

## Data

| Series           | Source                                                | Series ID    |
|------------------|-------------------------------------------------------|--------------|
| Unemployment rate (SA, 15+, persons)        | ABS 6202.0 Table 001    | A84423050A   |
| Participation rate (SA, 15+, persons)       | ABS 6202.0 Table 001    | A84423051C   |
| Civilian population aged 15+ (original)     | ABS 6202.0 Table 001    | A84423091W   |
| Real GDP (chain volume, SA, $m)             | ABS 5206.0 Table 1      | A2304402X    |

Monthly Labour Force series are averaged within calendar quarters
(complete quarters only). Real GDP per capita is GDP / civilian
population. For the 1966Q3-1977Q4 period (where the ABS Labour Force
table starts in 1978M2), historical spliced values from
[`data/lf_hist.xlsx`](data/lf_hist.xlsx) are used; a multiplicative
rebasing factor (fixed at the 1978Q2 splice point) keeps the level of
log real GDP per capita continuous.

## Repository layout

```
.
|-- stan/lf_uc.stan       # Stan model with Kalman filter + smoother
|-- R/abs_fetch.R         # Direct download of ABS xlsx tables
|-- R/data_prep.R         # Splice + quarterly aggregation
|-- R/fit.R               # cmdstanr driver
|-- R/outputs.R           # CSV + chart writers
|-- R/run.R               # End-to-end pipeline entrypoint
|-- data/lf_hist.xlsx     # Bundled historical data (1964Q1-2019Q2)
|-- outputs/              # Generated each run (committed)
|-- docs/                 # GitHub Pages site (chart, downloads)
|-- tests/                # Pipeline tests
|-- reference/eviews/     # Original EViews program (for provenance)
`-- .github/workflows/    # Monthly refresh action
```

## Running locally

Requirements: R >= 4.2, [CmdStan](https://mc-stan.org/users/interfaces/cmdstan)
installed via `cmdstanr::install_cmdstan()`, and the R packages listed in
[`DESCRIPTION`](DESCRIPTION).

```bash
Rscript R/run.R
```

This downloads the latest ABS workbooks into `data/raw/`, fits four
HMC chains (~5-10 minutes on a modern laptop), and writes refreshed
outputs to `outputs/` and `docs/`.

Tuning knobs via env vars: `PRD_CHAINS`, `PRD_WARMUP`, `PRD_SAMPLING`,
`PRD_CACHE_DIR`.

## Automatic updates

[`.github/workflows/update.yml`](.github/workflows/update.yml) runs the
pipeline at 02:30 UTC on the 14th-22nd of each month -- the window
during which the ABS publishes the monthly Labour Force release
(typically the third Thursday). When the run produces a new dataset
the workflow commits the refreshed `outputs/` and `docs/` and
redeploys GitHub Pages.

You can also trigger it manually from the **Actions** tab via
"Run workflow".

## Tests

```bash
Rscript -e 'testthat::test_dir("tests")'
```

The tests cover (i) the ABS xlsx parser against a small bundled
fixture and (ii) the quarterly aggregation / splice logic.

## License

MIT. The original methodology is due to Evans, Moore & Rees (2018);
this repository contains an independent open-source implementation.
