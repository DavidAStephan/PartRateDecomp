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
// Differences from the EViews replica:
//
// 1. The AR(2) cycle is reparameterised via partial autocorrelations
//    (Barndorff-Nielsen-Schou / Jones 1987) so the cycle is guaranteed
//    stationary. The original EViews ML estimate sits at phi1+phi2~1 on
//    the extended sample, which makes the cycle a near-random-walk and
//    causes it to absorb the persistent post-2000 drift in Australian
//    output and participation. Hard stationarity forces that drift back
//    into the random-walk trends, where it belongs.
//
// 2. Cycle innovations follow a Student-t via the standard
//    scale-mixture-of-normals representation: e4[t] | lambda[t] ~ N(0,
//    sigma4^2 * lambda[t]) with lambda[t] ~ inv_gamma(nu/2, nu/2). This
//    makes the model robust to one-off shocks (COVID-19 in particular).
//    The correlation cov(e1, e4) scales as sigma1 sigma4 rho sqrt(lambda).
//
// Variances and the rho correlation follow a regime switch at t =
// break_idx (1983Q4 default); observations with t <= break_idx use the
// "pre" parameters, the rest use "post".

data {
  int<lower=1> T;                     // number of quarterly observations
  vector[T] y;                        // log real GDP per capita * 100
  vector[T] u;                        // unemployment rate (%)
  vector[T] p;                        // participation rate (%)
  int<lower=1, upper=T> break_idx;    // last "pre-break" index (1983Q4)
  vector[5] m0;                       // prior mean of initial state
  cov_matrix[5] P0;                   // prior cov  of initial state

  real prior_delta_mean;
  real<lower=0> prior_delta_sd;
}

transformed data {
  matrix[3, 3] H = diag_matrix(rep_vector(1e-6, 3));
}

parameters {
  real delta;
  // Stationary AR(2) via partial autocorrelations on (-1, 1).
  real<lower=-1, upper=1> pacf1;
  real<lower=-1, upper=1> pacf2;

  // Sign constraints break the (cycle <-> -cycle) identification ambiguity.
  // Okun's law: a positive output cycle lowers unemployment, so kappa1 <= 0.
  // Participation is mildly procyclical in Australia, so theta1 >= 0.
  // kappa2 and theta2 are free in sign (lag terms can go either way).
  real<upper=0> kappa1;
  real kappa2;
  real<lower=0> theta1;
  real theta2;

  // [1] = pre-break (t <= break_idx), [2] = post-break
  array[2] vector<lower=0>[4] sigma;
  array[2] real<lower=-0.99, upper=0.99> rho;

  // Cycle Student-t -> scale-mixture-of-normals
  real<lower=2> nu;                   // degrees of freedom
  vector<lower=0>[T] lambda;          // per-period scale factor
}

transformed parameters {
  // Levinson-Durbin: AR(2) coefficients from partial autocorrelations.
  // For any pacf1, pacf2 in (-1, 1) the resulting AR(2) is stationary.
  real phi1 = pacf1 * (1 - pacf2);
  real phi2 = pacf2;

  matrix[5, 5] T_mat;
  matrix[3, 5] Z;
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
}

model {
  // ----------------- Priors -----------------
  delta ~ normal(prior_delta_mean, prior_delta_sd);

  // Stationarity-respecting cycle priors. Centred on moderate persistence
  // so the cycle mean-reverts between recessions (matching the published
  // Evans/Moore/Rees Graph 3 behaviour). The unconstrained likelihood has
  // a second mode at pacf1 ~ 1 that produces too-slow recoveries on the
  // extended sample; the prior here puts most mass at the lower-persistence
  // mode that matches the paper's between-recession dynamics.
  pacf1 ~ normal(0.55, 0.15);
  pacf2 ~ normal(-0.2, 0.2);

  // Half-normal priors on the sign-constrained loadings.
  kappa1 ~ normal(0, 1);   // truncated to (-inf, 0]
  kappa2 ~ normal(0, 1);
  theta1 ~ normal(0, 1);   // truncated to [0, +inf)
  theta2 ~ normal(0, 1);

  for (r in 1:2) {
    sigma[r] ~ student_t(4, 0, 1);
    rho[r] ~ normal(0, 0.5);
  }

  // Robust cycle: scale mixture giving Student-t with df = nu. Tight prior
  // (gamma(16, 0.8): mean 20, sd 5) keeps the model close to Gaussian for
  // typical quarters so the cycle innovation scale is pinned down by data;
  // genuine outliers still pull lambda[t] large without collapsing nu.
  nu ~ gamma(16, 0.8);
  lambda ~ inv_gamma(nu / 2, nu / 2); // E[lambda]=1 once nu>2

  // ----------------- Kalman filter likelihood -----------------
  vector[5] m = m0;
  matrix[5, 5] P = P0;
  matrix[5, 3] Zt = Z';

  for (t in 1:T) {
    int r = (t <= break_idx) ? 1 : 2;

    // Time-varying Q: cycle innovation variance scaled by lambda[t]; the
    // (1,4) cross-term scales by sqrt(lambda[t]) so the joint (e1, e4) is
    // a valid (scale-mixed) bivariate normal.
    matrix[5, 5] Q_t = rep_matrix(0, 5, 5);
    real sl = sqrt(lambda[t]);
    Q_t[1, 1] = square(sigma[r][1]);
    Q_t[2, 2] = square(sigma[r][2]);
    Q_t[3, 3] = square(sigma[r][3]);
    Q_t[4, 4] = square(sigma[r][4]) * lambda[t];
    real cov14 = sigma[r][1] * sigma[r][4] * rho[r] * sl;
    Q_t[1, 4] = cov14;
    Q_t[4, 1] = cov14;

    vector[5] m_pred = T_mat * m + c_vec;
    matrix[5, 5] P_pred = quad_form_sym(P, T_mat') + Q_t;
    P_pred = 0.5 * (P_pred + P_pred');

    vector[3] y_t = [y[t], u[t], p[t]]';
    vector[3] v = y_t - Z * m_pred;
    matrix[3, 3] F = quad_form_sym(P_pred, Zt) + H;
    F = 0.5 * (F + F');
    matrix[3, 3] L = cholesky_decompose(F);

    target += -0.5 * (3 * log(2 * pi())
                      + 2 * sum(log(diagonal(L)))
                      + dot_self(mdivide_left_tri_low(L, v)));

    matrix[3, 5] ZP = Z * P_pred;
    matrix[5, 3] K = mdivide_left_spd(F, ZP)';
    m = m_pred + K * v;
    P = P_pred - K * F * K';
    P = 0.5 * (P + P');
  }
}

generated quantities {
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
    array[T] vector[5] m_filt;
    array[T] matrix[5, 5] P_filt;
    array[T] vector[5] m_pred_arr;
    array[T] matrix[5, 5] P_pred_arr;
    array[T] vector[5] m_smooth_arr;
    array[T] matrix[5, 5] P_smooth_arr;

    vector[5] m = m0;
    matrix[5, 5] P = P0;
    matrix[5, 3] Zt = Z';

    for (t in 1:T) {
      int r = (t <= break_idx) ? 1 : 2;

      matrix[5, 5] Q_t = rep_matrix(0, 5, 5);
      real sl = sqrt(lambda[t]);
      Q_t[1, 1] = square(sigma[r][1]);
      Q_t[2, 2] = square(sigma[r][2]);
      Q_t[3, 3] = square(sigma[r][3]);
      Q_t[4, 4] = square(sigma[r][4]) * lambda[t];
      real cov14 = sigma[r][1] * sigma[r][4] * rho[r] * sl;
      Q_t[1, 4] = cov14;
      Q_t[4, 1] = cov14;

      vector[5] m_pred = T_mat * m + c_vec;
      matrix[5, 5] P_pred = quad_form_sym(P, T_mat') + Q_t;
      P_pred = 0.5 * (P_pred + P_pred');
      m_pred_arr[t] = m_pred;
      P_pred_arr[t] = P_pred;

      vector[3] y_t = [y[t], u[t], p[t]]';
      vector[3] v = y_t - Z * m_pred;
      matrix[3, 3] F = quad_form_sym(P_pred, Zt) + H;
      F = 0.5 * (F + F');
      matrix[3, 5] ZP = Z * P_pred;
      matrix[5, 3] K = mdivide_left_spd(F, ZP)';
      m = m_pred + K * v;
      P = P_pred - K * F * K';
      P = 0.5 * (P + P');
      m_filt[t] = m;
      P_filt[t] = P;
    }

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
