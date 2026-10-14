# ==============================================================================
# 05_stress_tests.R: Stress Testing INLA Against Numerically Exact Benchmarks
# ==============================================================================
source("R/common.R")
suppressPackageStartupMessages({
  library(statmod) # For Gauss-Hermite quadrature nodes and weights
})

cat_header("05: Stress Testing INLA (Non-Gaussian & Small-Cluster Regimes)")

# Ensure statmod is available
if (!requireNamespace("statmod", quietly = TRUE)) {
  install.packages("statmod", repos = "https://cloud.r-project.org")
  library(statmod)
}

# ------------------------------------------------------------------------------
# SCENARIO S1: Low-Count Poisson with Observation-Level Random Effects
# ------------------------------------------------------------------------------
cat("\n=== SCENARIO S1: Low-Count Poisson with Observation-Level Random Effects ===\n")
cat("Description: n = 50 observations, low rate (mostly 0s and 1s), observation-level iid latent field.\n")

set.seed(12345)
n_s1 <- 50
x_s1 <- seq(-1, 1, length.out = n_s1)
b0_true_s1 <- -1.5
b1_true_s1 <- 0.8
sigma_u_s1 <- 0.7
u_true_s1  <- rnorm(n_s1, 0, sigma_u_s1)
eta_s1     <- b0_true_s1 + b1_true_s1 * x_s1 + u_true_s1
y_s1       <- rpois(n_s1, exp(eta_s1))

cat(sprintf("Observed counts: min = %d, max = %d, median = %d, percent zeros = %.1f%%\n",
            min(y_s1), max(y_s1), median(y_s1), 100 * mean(y_s1 == 0)))
d_s1 <- data.frame(y = y_s1, x = x_s1, idx = seq_len(n_s1))

# Gauss-Hermite quadrature rule for 1D inner integral (50 points)
gh50 <- statmod::gauss.quad(50, kind = "hermite")
# Transform standard Hermite (weight exp(-x^2)) to standard Normal (weight exp(-z^2/2) / sqrt(2*pi))
# z = sqrt(2) * x, weights scaled by 1/sqrt(pi)
gh_nodes_std <- sqrt(2) * gh50$nodes
gh_weights_std <- gh50$weights / sqrt(pi)

# Inner likelihood for one observation y_i given (beta0, beta1, sigma_u)
quad_pois_1obs <- function(yi, xi, b0, b1, su) {
  mu_eval <- exp(b0 + b1 * xi + su * gh_nodes_std)
  # dpois on each node
  lik_nodes <- dpois(yi, mu_eval)
  sum(lik_nodes * gh_weights_std)
}

# Log conditional likelihood of y given (b0, b1, su)
ll_s1_cond <- function(b0, b1, su) {
  sum(vapply(seq_len(n_s1), function(i) {
    log(max(quad_pois_1obs(y_s1[i], x_s1[i], b0, b1, su), 1e-300))
  }, 0.0))
}

# Exact 3D integration over (b0, b1, log_tau)
cat("Computing numerically exact marginal likelihood via adaptive 3D quadrature... ")
# Find mode of (b0, b1, log_tau)
lpost_s1 <- function(par) {
  b0 <- par[1]; b1 <- par[2]; ltau <- par[3]
  tau <- exp(ltau); su <- 1 / sqrt(tau)
  ll <- ll_s1_cond(b0, b1, su)
  lprior_b <- dnorm(b0, 0, 10, log = TRUE) + dnorm(b1, 0, 10, log = TRUE)
  lprior_tau <- dgamma(tau, 1, 1, log = TRUE) + ltau # Jacobian
  ll + lprior_b + lprior_tau
}

opt_s1 <- optim(c(-1.5, 0.8, log(2)), function(p) -lpost_s1(p), method = "BFGS", hessian = TRUE)
mode_s1 <- opt_s1$par
H_s1 <- opt_s1$hessian
cov_s1 <- solve(H_s1)
sd_s1 <- sqrt(diag(cov_s1))

# Fine 3D grid integration centered at mode (+/- 5 standard deviations, 35 points per dimension)
n_grid_s1 <- 35
g_b0   <- seq(mode_s1[1] - 5 * sd_s1[1], mode_s1[1] + 5 * sd_s1[1], length.out = n_grid_s1)
g_b1   <- seq(mode_s1[2] - 5 * sd_s1[2], mode_s1[2] + 5 * sd_s1[2], length.out = n_grid_s1)
g_ltau <- seq(mode_s1[3] - 5 * sd_s1[3], mode_s1[3] + 5 * sd_s1[3], length.out = n_grid_s1)

h_b0 <- g_b0[2] - g_b0[1]
h_b1 <- g_b1[2] - g_b1[1]
h_lt <- g_ltau[2] - g_ltau[1]

grid_s1_pts <- expand.grid(b0 = g_b0, b1 = g_b1, ltau = g_ltau)
vals_s1 <- unlist(mclapply(seq_len(nrow(grid_s1_pts)), function(row_idx) {
  lpost_s1(c(grid_s1_pts$b0[row_idx], grid_s1_pts$b1[row_idx], grid_s1_pts$ltau[row_idx]))
}, mc.cores = N_CORES))

m_max_s1 <- max(vals_s1)
exact_mlik_s1 <- m_max_s1 + log(sum(exp(vals_s1 - m_max_s1)) * h_b0 * h_b1 * h_lt)
cat(sprintf("Exact = %.4f\n", exact_mlik_s1))

# Fit INLA on S1 under various configurations
f_s1 <- y ~ x + f(idx, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 1))))
configs_s1 <- list(
  list(name = "INLA Default (Simplified Laplace + CCD)", strat = "simplified.laplace", int = "ccd"),
  list(name = "INLA Gaussian Approximation",            strat = "gaussian",           int = "ccd"),
  list(name = "INLA Full Laplace + CCD",                strat = "laplace",            int = "ccd"),
  list(name = "INLA Full Laplace + Grid",               strat = "laplace",            int = "grid")
)

res_s1 <- list()
for (cfg in configs_s1) {
  m_inla <- inla(f_s1, data = d_s1, family = "poisson",
                 control.compute = list(mlik = TRUE),
                 control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01, prec.intercept = 0.01),
                 control.inla = list(strategy = cfg$strat, int.strategy = cfg$int))
  res_s1[[length(res_s1) + 1]] <- data.frame(
    Scenario = "S1: Low-Count Poisson GLMM",
    Method = cfg$name,
    Strategy = cfg$strat,
    Integration = cfg$int,
    LogMargLik = round(m_inla$mlik[1, 1], 4),
    LogMargLik_Gau = round(m_inla$mlik[2, 1], 4),
    Exact_Ref = round(exact_mlik_s1, 4),
    Discrepancy = round(m_inla$mlik[1, 1] - exact_mlik_s1, 4)
  )
}

# ------------------------------------------------------------------------------
# SCENARIO S2: Small-Cluster Binary Logistic GLMM
# ------------------------------------------------------------------------------
cat("\n=== SCENARIO S2: Small-Cluster Binary Logistic GLMM ===\n")
cat("Description: J = 20 clusters of size K = 3 (N = 60), binary outcome, cluster random intercept.\n")

set.seed(54321)
J_s2 <- 20
K_s2 <- 3
N_s2 <- J_s2 * K_s2
cluster_id <- rep(seq_len(J_s2), each = K_s2)
x_s2 <- runif(N_s2, -1.5, 1.5)
b0_true_s2 <- -0.3
b1_true_s2 <- 1.0
sigma_u_s2 <- 1.2
u_cluster  <- rnorm(J_s2, 0, sigma_u_s2)
eta_s2     <- b0_true_s2 + b1_true_s2 * x_s2 + u_cluster[cluster_id]
p_s2       <- plogis(eta_s2)
y_s2       <- rbinom(N_s2, 1, p_s2)

cat(sprintf("Observed binary outcomes: %d ones, %d zeros (overall prevalence = %.1f%%)\n",
            sum(y_s2), sum(1 - y_s2), 100 * mean(y_s2)))
d_s2 <- data.frame(y = y_s2, x = x_s2, cluster = cluster_id)

# Inner likelihood for cluster j given (b0, b1, su)
quad_binom_1cluster <- function(j, b0, b1, su) {
  idx <- which(cluster_id == j)
  yj <- y_s2[idx]; xj <- x_s2[idx]
  # For each quadrature node, calculate product over k = 1..K
  lik_nodes <- vapply(seq_along(gh_nodes_std), function(q) {
    z <- gh_nodes_std[q]
    etaj <- b0 + b1 * xj + su * z
    prod(plogis(etaj)^yj * (1 - plogis(etaj))^(1 - yj))
  }, 0.0)
  sum(lik_nodes * gh_weights_std)
}

ll_s2_cond <- function(b0, b1, su) {
  sum(vapply(seq_len(J_s2), function(j) {
    log(max(quad_binom_1cluster(j, b0, b1, su), 1e-300))
  }, 0.0))
}

cat("Computing numerically exact marginal likelihood for S2 via 3D quadrature... ")
lpost_s2 <- function(par) {
  b0 <- par[1]; b1 <- par[2]; ltau <- par[3]
  tau <- exp(ltau); su <- 1 / sqrt(tau)
  ll <- ll_s2_cond(b0, b1, su)
  lprior_b <- dnorm(b0, 0, 10, log = TRUE) + dnorm(b1, 0, 10, log = TRUE)
  lprior_tau <- dgamma(tau, 1, 1, log = TRUE) + ltau
  ll + lprior_b + lprior_tau
}

opt_s2 <- optim(c(-0.3, 1.0, log(1/1.2^2)), function(p) -lpost_s2(p), method = "BFGS", hessian = TRUE)
mode_s2 <- opt_s2$par
cov_s2 <- solve(opt_s2$hessian)
sd_s2 <- sqrt(diag(cov_s2))

n_grid_s2 <- 35
g_b0_2   <- seq(mode_s2[1] - 5 * sd_s2[1], mode_s2[1] + 5 * sd_s2[1], length.out = n_grid_s2)
g_b1_2   <- seq(mode_s2[2] - 5 * sd_s2[2], mode_s2[2] + 5 * sd_s2[2], length.out = n_grid_s2)
g_ltau_2 <- seq(mode_s2[3] - 5 * sd_s2[3], mode_s2[3] + 5 * sd_s2[3], length.out = n_grid_s2)

h_b0_2 <- g_b0_2[2] - g_b0_2[1]
h_b1_2 <- g_b1_2[2] - g_b1_2[1]
h_lt_2 <- g_ltau_2[2] - g_ltau_2[1]

grid_s2_pts <- expand.grid(b0 = g_b0_2, b1 = g_b1_2, ltau = g_ltau_2)
vals_s2 <- unlist(mclapply(seq_len(nrow(grid_s2_pts)), function(row_idx) {
  lpost_s2(c(grid_s2_pts$b0[row_idx], grid_s2_pts$b1[row_idx], grid_s2_pts$ltau[row_idx]))
}, mc.cores = N_CORES))

m_max_s2 <- max(vals_s2)
exact_mlik_s2 <- m_max_s2 + log(sum(exp(vals_s2 - m_max_s2)) * h_b0_2 * h_b1_2 * h_lt_2)
cat(sprintf("Exact = %.4f\n", exact_mlik_s2))

# Fit INLA on S2
f_s2 <- y ~ x + f(cluster, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 1))))
configs_s2 <- list(
  list(name = "INLA Default (Simplified Laplace + CCD)", strat = "simplified.laplace", int = "ccd"),
  list(name = "INLA Gaussian Approximation",            strat = "gaussian",           int = "ccd"),
  list(name = "INLA Full Laplace + CCD",                strat = "laplace",            int = "ccd"),
  list(name = "INLA Full Laplace + Grid",               strat = "laplace",            int = "grid")
)

res_s2 <- list()
for (cfg in configs_s2) {
  m_inla <- inla(f_s2, data = d_s2, family = "binomial", Ntrials = rep(1, N_s2),
                 control.compute = list(mlik = TRUE),
                 control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01, prec.intercept = 0.01),
                 control.inla = list(strategy = cfg$strat, int.strategy = cfg$int))
  res_s2[[length(res_s2) + 1]] <- data.frame(
    Scenario = "S2: Small-Cluster Binary Logistic GLMM",
    Method = cfg$name,
    Strategy = cfg$strat,
    Integration = cfg$int,
    LogMargLik = round(m_inla$mlik[1, 1], 4),
    LogMargLik_Gau = round(m_inla$mlik[2, 1], 4),
    Exact_Ref = round(exact_mlik_s2, 4),
    Discrepancy = round(m_inla$mlik[1, 1] - exact_mlik_s2, 4)
  )
}

df_stress <- rbind(do.call(rbind, res_s1), do.call(rbind, res_s2))
write.csv(df_stress, "results/stress_tests.csv", row.names = FALSE)
print(df_stress[, c("Scenario", "Method", "LogMargLik", "Exact_Ref", "Discrepancy")])

cat("\n05_stress_tests.R completed successfully. Results saved to results/stress_tests.csv\n")
