// Unobserved-components state-space model of the Australian labour market.
//
// Reference: Evans, Moore & Rees (2018), "The Cyclical Behaviour of Labour Force
// Participation", RBA Bulletin (September).
//
// State vector x[t] = (ystar, unrstar, prtstar, cycle, cyclelag).
//
// Signal equations (no measurement error -- a tiny diagonal regularises the
// forecast covariance for numerical stability):
//   lrgdppc[t] = ystar[t] + cycle[t]
//   unr[t]    = unrstar[t] + kappa1 * cycle[t] + kappa2 * cyclelag[t]
//   prt[t]    = prtstar[t] + theta1 * cycle[t] + theta2 * cyclelag[t]
//
// State equations:
//   ystar[t]    = delta + ystar[t-1] + e1[t]
//   unrstar[t]  = unrstar[t-1] + e2[t]
//   prtstar[t]  = prtstar[t-1] + e3[t]
//   cycle[t]    = phi1 * cycle[t-1] + phi2 * cyclelag[t-1] + e4[t]
//   cyclelag[t] = cycle[t-1]
//
// Innovations (e1, e2, e3, e4) are jointly Normal with diagonal covariance
// except cov(e1, e4) = sigma1 * sigma4 * rho. Variances and the rho correlation
// follow a regime switch at t = break_idx (defaults to 1983Q4): observations
// with t <= break_idx use the "pre" parameters, the rest use "post".

// (Kalman filter is inlined below for speed; no functions block.)

data {
  int<lower=1> T;                     // number of quarterly observations
  vector[T] y;                        // log real GDP per capita * 100
  vector[T] u;                        // unemployment rate (%)
  vector[T] p;                        // participation rate (%)
  int<lower=1, upper=T> break_idx;    // last "pre-break" index (1983Q4)
  vector[5] m0;                       // prior mean of initial state
  cov_matrix[5] P0;                   // prior cov  of initial state

  // Optional informative prior centres for hyperparameters (use 0 for diffuse)
  real prior_delta_mean;
  real<lower=0> prior_delta_sd;
}

transformed data {
  matrix[3, 3] H = diag_matrix(rep_vector(1e-6, 3));
}

parameters {
  real delta;                         // drift in trend output
  real<lower=-2, upper=2> phi1;       // AR(2) coeffs on cycle
  real<lower=-1, upper=1> phi2;
  real kappa1;                        // Okun-type loadings
  real kappa2;
  real theta1;                        // participation loadings on cycle
  real theta2;

  // [1] = pre-break (t <= break_idx), [2] = post-break
  array[2] vector<lower=0>[4] sigma;
  array[2] real<lower=-0.99, upper=0.99> rho;
}

transformed parameters {
  // Stationarity of the AR(2) for the cycle: roots outside the unit circle
  // require phi1 + phi2 < 1, phi2 - phi1 < 1, |phi2| < 1.
  // We don't enforce this hard (the prior will pull) -- but you could.

  // Build state-space matrices for the two regimes.
  matrix[5, 5] T_mat;
  matrix[3, 5] Z;
  array[2] matrix[5, 5] Q;
  vector[5] c_vec;

  T_mat = rep_matrix(0, 5, 5);
  T_mat[1, 1] = 1;
  T_mat[2, 2] = 1;
  T_mat[3, 3] = 1;
  T_mat[4, 4] = phi1;
  T_mat[4, 5] = phi2;
  T_mat[5, 4] = 1;

  Z = rep_matrix(0, 3, 5);
  Z[1, 1] = 1;  Z[1, 4] = 1;
  Z[2, 2] = 1;  Z[2, 4] = kappa1;  Z[2, 5] = kappa2;
  Z[3, 3] = 1;  Z[3, 4] = theta1;  Z[3, 5] = theta2;

  c_vec = rep_vector(0, 5);
  c_vec[1] = delta;

  for (r in 1:2) {
    Q[r] = rep_matrix(0, 5, 5);
    Q[r][1, 1] = square(sigma[r][1]);
    Q[r][2, 2] = square(sigma[r][2]);
    Q[r][3, 3] = square(sigma[r][3]);
    Q[r][4, 4] = square(sigma[r][4]);
    real cov14 = sigma[r][1] * sigma[r][4] * rho[r];
    Q[r][1, 4] = cov14;
    Q[r][4, 1] = cov14;
  }
}

model {
  // ----------------- Priors -----------------
  // Drift: average quarterly log-pc-GDP growth ~ 0.5% (i.e. ~2% annual)
  delta ~ normal(prior_delta_mean, prior_delta_sd);

  // AR(2) cycle: weakly favours persistence, second lag small/negative
  phi1 ~ normal(1.3, 0.5);
  phi2 ~ normal(-0.4, 0.5);

  // Okun / participation loadings: weakly informative around zero
  kappa1 ~ normal(0, 1);
  kappa2 ~ normal(0, 1);
  theta1 ~ normal(0, 1);
  theta2 ~ normal(0, 1);

  // Innovation scales: half-Student-t -- robust, weakly informative
  for (r in 1:2) {
    sigma[r] ~ student_t(4, 0, 1);
    rho[r] ~ normal(0, 0.5);
  }

  // ----------------- Kalman filter likelihood -----------------
  // Inlined for speed: Cholesky of forecast covariance avoids the explicit
  // matrix solves used by the textbook recursion.
  vector[5] m = m0;
  matrix[5, 5] P = P0;
  matrix[5, 3] Zt = Z';

  for (t in 1:T) {
    int r = (t <= break_idx) ? 1 : 2;

    // Predict
    vector[5] m_pred = T_mat * m + c_vec;
    matrix[5, 5] P_pred = quad_form_sym(P, T_mat') + Q[r];
    P_pred = 0.5 * (P_pred + P_pred');

    // Forecast and innovation
    vector[3] y_t = [y[t], u[t], p[t]]';
    vector[3] v = y_t - Z * m_pred;
    matrix[3, 3] F = quad_form_sym(P_pred, Zt) + H;
    F = 0.5 * (F + F');
    matrix[3, 3] L = cholesky_decompose(F);

    // Log marginal likelihood contribution
    target += -0.5 * (3 * log(2 * pi())
                      + 2 * sum(log(diagonal(L)))
                      + dot_self(mdivide_left_tri_low(L, v)));

    // Update -- K = P_pred * Z' * F^{-1}; using the Cholesky factor of F.
    matrix[3, 5] ZP = Z * P_pred;
    matrix[5, 3] K = mdivide_left_spd(F, ZP)';
    m = m_pred + K * v;
    P = P_pred - K * F * K';
    P = 0.5 * (P + P');
  }
}

generated quantities {
  // Exposed outputs only: smoothed state means and standard errors.
  vector[T] ystar_sm;
  vector[T] unrstar_sm;
  vector[T] prtstar_sm;
  vector[T] cycle_sm;
  vector[T] cyclelag_sm;
  vector[T] prt_cycle_sm;
  vector[T] ystar_se;
  vector[T] unrstar_se;
  vector[T] prtstar_se;
  vector[T] cycle_se;
  vector[T] prt_cycle_se;

  {
    // All filter/smoother working memory stays local so it does not bloat the
    // posterior draws output.
    array[T] vector[5] m_filt;
    array[T] matrix[5, 5] P_filt;
    array[T] vector[5] m_pred_arr;
    array[T] matrix[5, 5] P_pred_arr;
    array[T] vector[5] m_smooth_arr;
    array[T] matrix[5, 5] P_smooth_arr;

    vector[5] m = m0;
    matrix[5, 5] P = P0;
    matrix[5, 3] Zt = Z';

    // Forward filter
    for (t in 1:T) {
      int r = (t <= break_idx) ? 1 : 2;
      vector[5] m_pred = T_mat * m + c_vec;
      matrix[5, 5] P_pred = quad_form_sym(P, T_mat') + Q[r];
      P_pred = 0.5 * (P_pred + P_pred');
      m_pred_arr[t] = m_pred;
      P_pred_arr[t] = P_pred;

      vector[3] y_t = [y[t], u[t], p[t]]';
      vector[3] v = y_t - Z * m_pred;
      matrix[3, 3] F = quad_form_sym(P_pred, Zt) + H;
      F = 0.5 * (F + F');
      matrix[5, 3] PZt = P_pred * Zt;
      matrix[5, 3] K = PZt / F;
      m = m_pred + K * v;
      P = P_pred - K * F * K';
      P = 0.5 * (P + P');
      m_filt[t] = m;
      P_filt[t] = P;
    }

    // Backward (RTS) smoother
    m_smooth_arr[T] = m_filt[T];
    P_smooth_arr[T] = P_filt[T];
    for (k in 1:(T - 1)) {
      int t = T - k;
      matrix[5, 5] J = P_filt[t] * T_mat' / P_pred_arr[t + 1];
      m_smooth_arr[t] = m_filt[t]
                       + J * (m_smooth_arr[t + 1] - m_pred_arr[t + 1]);
      P_smooth_arr[t] = P_filt[t]
                       + J * (P_smooth_arr[t + 1] - P_pred_arr[t + 1]) * J';
      P_smooth_arr[t] = 0.5 * (P_smooth_arr[t] + P_smooth_arr[t]');
    }

    for (t in 1:T) {
      ystar_sm[t]    = m_smooth_arr[t][1];
      unrstar_sm[t]  = m_smooth_arr[t][2];
      prtstar_sm[t]  = m_smooth_arr[t][3];
      cycle_sm[t]    = m_smooth_arr[t][4];
      cyclelag_sm[t] = m_smooth_arr[t][5];
      prt_cycle_sm[t] = theta1 * cycle_sm[t] + theta2 * cyclelag_sm[t];
      ystar_se[t]    = sqrt(fmax(0.0, P_smooth_arr[t][1, 1]));
      unrstar_se[t]  = sqrt(fmax(0.0, P_smooth_arr[t][2, 2]));
      prtstar_se[t]  = sqrt(fmax(0.0, P_smooth_arr[t][3, 3]));
      cycle_se[t]    = sqrt(fmax(0.0, P_smooth_arr[t][4, 4]));
      real v44 = P_smooth_arr[t][4, 4];
      real v55 = P_smooth_arr[t][5, 5];
      real v45 = P_smooth_arr[t][4, 5];
      prt_cycle_se[t] = sqrt(fmax(0.0,
                                  square(theta1) * v44
                                  + square(theta2) * v55
                                  + 2 * theta1 * theta2 * v45));
    }
  }
}
