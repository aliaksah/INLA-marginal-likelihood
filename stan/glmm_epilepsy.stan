// Stan model for Longitudinal Poisson Mixed Model (Epileptic Seizures)
// Matching Chib & Jeliazkov (2001) / INLA Section 6 formulation exactly
data {
  int<lower=1> N;                 // total observations (290)
  int<lower=1> J;                 // number of patients (58)
  array[N] int<lower=0> y;        // seizure count
  vector[N] ofs;                  // log offset
  vector[N] trt;                  // treatment indicator
  vector[N] bsl;                  // baseline visit indicator
  array[N] int<lower=1, upper=J> id; // patient id
}

parameters {
  vector[4] beta;                 // fixed effects: intercept, trt, trt*bsl, bsl
  matrix[2, J] u;                 // random effects for each patient
  cov_matrix[2] Sigma;            // covariance matrix of random effects
}

transformed parameters {
  vector[N] eta;
  for (n in 1:N) {
    eta[n] = ofs[n] + beta[1] + beta[2] * trt[n] + beta[3] * (trt[n] * bsl[n]) +
             beta[4] * bsl[n] + u[1, id[n]] + u[2, id[n]] * bsl[n];
  }
}

model {
  // Priors on fixed effects: beta ~ N(0, 100) (prec = 0.01)
  beta ~ normal(0, 10);
  
  // Prior on precision matrix: Q = Sigma^-1 ~ Wishart(4, I_2)
  // Corresponding to Sigma ~ InvWishart(4, I_2)
  Sigma ~ inv_wishart(4, diag_matrix(rep_vector(1.0, 2)));
  
  // Random effects distribution
  for (j in 1:J) {
    u[, j] ~ multi_normal(rep_vector(0.0, 2), Sigma);
  }
  
  // Poisson likelihood
  y ~ poisson_log(eta);
}
