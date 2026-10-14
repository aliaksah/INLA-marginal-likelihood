# ==============================================================================
# 01_tables.R: Reproduce Tables 1, 2, 3, 4 and Decision-Scale Summary (Table 5)
# ==============================================================================
source("R/common.R")

cat_header("01: Tables 1-4 Reproduction and Decision-Scale Analysis")

# ------------------------------------------------------------------------------
# 1. TABLE 1: Simple linear latent model (Section 2)
# ------------------------------------------------------------------------------
cat("\n--- TABLE 1: Simple Linear Latent Model (Section 2) ---\n")
spe <- 2 # D = 2 is the observed single datum y = 2
tab1_s0 <- c(1000, 10, 0.1)
tab1_exact <- numeric(length(tab1_s0))
tab1_inla  <- numeric(length(tab1_s0))
tab1_hm    <- matrix(NA, nrow = length(tab1_s0), ncol = N_REPS)

true_marg_lik <- function(x, s0, s1 = 1) log(dnorm(x, 0, sqrt(s0^2 + s1^2)))
hm_marg_lik <- function(x, s0, s1 = 1, n = if (QUICK) 1e6 else 1e7) {
  post_prec <- 1/s0^2 + 1/s1^2
  t <- rnorm(n, (x/s1^2)/post_prec, sqrt(1/post_prec))
  lik <- dnorm(x, t, s1)
  log(1 / mean(1 / lik))
}

for (k in seq_along(tab1_s0)) {
  s0 <- tab1_s0[k]
  dat <- data.frame(y = spe, x = 1)
  out <- inla(y ~ x - 1, data = dat, control.compute = list(mlik = TRUE),
              control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
              control.fixed = list(mean = 0, prec = (1/s0)^2))
  tab1_exact[k] <- true_marg_lik(spe, s0, 1)
  tab1_inla[k]  <- out$mlik[1, 1] # Row 1 and Row 2 are identical (0 hyperparameters)
  
  cat(sprintf("  Harmonic Mean for sigma0 = %-6g (%d reps)... ", s0, N_REPS))
  # Run reps with deterministic seeds
  hm_vals <- mclapply(seq_len(N_REPS), function(r) {
    seed <- as.integer(1000 * k + r)
    set.seed(seed)
    hm_marg_lik(spe, s0, 1)
  }, mc.cores = min(N_REPS, N_CORES))
  for (r in seq_len(N_REPS)) {
    tab1_hm[k, r] <- hm_vals[[r]]
    note_seed(sprintf("tab1_hm_s0_%g_rep%d", s0, r), 1000 * k + r)
  }
  cat("done.\n")
}

tab1_hm_mean <- rowMeans(tab1_hm)
tab1_hm_sd   <- apply(tab1_hm, 1, sd)
tab1_hm_err  <- tab1_hm_mean - tab1_exact # Error against exact analytical value

df_tab1 <- data.frame(
  sigma0 = tab1_s0,
  sigma1 = rep(1, length(tab1_s0)),
  D = rep(spe, length(tab1_s0)), # D column explained as observed datum y = 2
  Exact = round(tab1_exact, 4),
  INLA = round(tab1_inla, 4),
  HM_mean = round(tab1_hm_mean, 4),
  HM_sd   = round(tab1_hm_sd, 4),
  HM_error = round(tab1_hm_err, 4)
)
for (r in seq_len(N_REPS)) {
  df_tab1[[paste0("HM_rep", r)]] <- round(tab1_hm[, r], 4)
}
write.csv(df_tab1, "results/table1_linear.csv", row.names = FALSE)
print(df_tab1[, c("sigma0", "D", "Exact", "INLA", "HM_mean", "HM_sd", "HM_error")])

# ------------------------------------------------------------------------------
# 2. TABLE 2: Gaussian Bayesian Regression on US Crime Data (Section 3)
# ------------------------------------------------------------------------------
cat("\n--- TABLE 2: US Crime Data (Section 3) ---\n")
crime <- load_crime()
tab2_models <- c("M1", "M1", "M1", "M2", "M2", "M2")
tab2_sbeta  <- c(1000, 10, 0.1, 1000, 10, 0.1)
tab2_exact  <- numeric(6)
tab2_inla1  <- numeric(6) # Row 1 (integration)
tab2_inla2  <- numeric(6) # Row 2 (Gaussian)
tab2_chib   <- matrix(NA, nrow = 6, ncol = N_REPS)

mcmc_iters <- if (QUICK) 1e4 else 1e5
idx <- 1
for (mod in 1:2) {
  f_curr <- if (mod == 1) crime$f1 else crime$f2
  mod_name <- if (mod == 1) "M1" else "M2"
  for (s_beta in c(1000, 10, 0.1)) {
    precp <- (1/s_beta)^2
    
    # Exact analytical benchmark
    tab2_exact[idx] <- exact_gauss(f_curr, crime$data, s_beta)
    
    # INLA
    m <- inla_gauss(f_curr, crime$data, precp)
    tab2_inla1[idx] <- m$mlik[1, 1]
    tab2_inla2[idx] <- m$mlik[2, 1]
    
    cat(sprintf("  Chib's MCMCregress for %s, sbeta = %-6g (%d reps)... ", mod_name, s_beta, N_REPS))
    # Parallel Chib runs with distinct recorded seeds
    seeds_curr <- 2000 * idx + seq_len(N_REPS)
    chib_vals <- mclapply(seq_len(N_REPS), function(r) {
      chib_gauss(f_curr, crime$data, precp, seeds_curr[r], mcmc = mcmc_iters)
    }, mc.cores = min(N_REPS, N_CORES))
    for (r in seq_len(N_REPS)) {
      tab2_chib[idx, r] <- chib_vals[[r]]
      note_seed(sprintf("tab2_chib_%s_sb%g_rep%d", mod_name, s_beta, r), seeds_curr[r])
    }
    cat("done.\n")
    idx <- idx + 1
  }
}

tab2_chib_mean <- rowMeans(tab2_chib)
tab2_chib_sd   <- apply(tab2_chib, 1, sd)
tab2_chib_se   <- tab2_chib_sd / sqrt(N_REPS)
diff_inla_exact <- tab2_inla1 - tab2_exact
diff_chib_exact <- tab2_chib_mean - tab2_exact
diff_inla_chib  <- tab2_inla1 - tab2_chib_mean

df_tab2 <- data.frame(
  Model = tab2_models,
  sigma_beta = tab2_sbeta,
  Exact = round(tab2_exact, 4),
  INLA_int = round(tab2_inla1, 4),
  INLA_gau = round(tab2_inla2, 4),
  Chib_mean = round(tab2_chib_mean, 4),
  Chib_sd = round(tab2_chib_sd, 4),
  Chib_mcse = round(tab2_chib_se, 4),
  Diff_INLA_Exact = round(diff_inla_exact, 4),
  Diff_Chib_Exact = round(diff_chib_exact, 4),
  Diff_INLA_Chib = round(diff_inla_chib, 4)
)
for (r in seq_len(N_REPS)) {
  df_tab2[[paste0("Chib_rep", r)]] <- round(tab2_chib[, r], 4)
}
write.csv(df_tab2, "results/table2_crime.csv", row.names = FALSE)
print(df_tab2[, c("Model", "sigma_beta", "Exact", "INLA_int", "Chib_mean", "Chib_sd", "Diff_INLA_Exact", "Diff_INLA_Chib")])

# ------------------------------------------------------------------------------
# 3. TABLE 3: Probit Bayesian Regression (Section 4)
# ------------------------------------------------------------------------------
cat("\n--- TABLE 3: Simulated Probit Data (Section 4) ---\n")
probit <- load_probit()
tab3_models <- c("M1", "M1", "M1", "M2", "M2", "M2")
tab3_sbeta  <- c(1000, 10, 0.1, 1000, 10, 0.1)
tab3_inla   <- numeric(6)
tab3_is_ev  <- numeric(6)
tab3_is_se  <- numeric(6)
tab3_chib   <- matrix(NA, nrow = 6, ncol = N_REPS)

n_is_samples <- if (QUICK) 2e4 else 2e5
idx <- 1
for (mod in 1:2) {
  f_curr <- if (mod == 1) probit$f1 else probit$f2
  mod_name <- if (mod == 1) "M1" else "M2"
  for (s_beta in c(1000, 10, 0.1)) {
    precp <- (1/s_beta)^2
    
    # INLA (classic mode, row 1)
    m <- inla_binom(f_curr, probit$data, "probit", precp)
    tab3_inla[idx] <- m$mlik[1, 1]
    
    # Importance sampling benchmark
    is_seed <- 3000 * idx + 99
    is_res <- is_logev(f_curr, probit$data, link = "probit", sb = s_beta,
                       N = n_is_samples, seed = is_seed)
    tab3_is_ev[idx] <- is_res$logev
    tab3_is_se[idx] <- is_res$se
    note_seed(sprintf("tab3_is_%s_sb%g", mod_name, s_beta), is_seed)
    
    cat(sprintf("  Chib's MCMCprobit for %s, sbeta = %-6g (%d reps)... ", mod_name, s_beta, N_REPS))
    seeds_curr <- 4000 * idx + seq_len(N_REPS)
    chib_vals <- mclapply(seq_len(N_REPS), function(r) {
      chib_probit(f_curr, probit$data, precp, seeds_curr[r], mcmc = mcmc_iters)
    }, mc.cores = min(N_REPS, N_CORES))
    for (r in seq_len(N_REPS)) {
      tab3_chib[idx, r] <- chib_vals[[r]]
      note_seed(sprintf("tab3_chib_%s_sb%g_rep%d", mod_name, s_beta, r), seeds_curr[r])
    }
    cat("done.\n")
    idx <- idx + 1
  }
}

tab3_chib_mean <- rowMeans(tab3_chib)
tab3_chib_sd   <- apply(tab3_chib, 1, sd)
tab3_chib_se   <- tab3_chib_sd / sqrt(N_REPS)
diff_inla_is   <- tab3_inla - tab3_is_ev
diff_chib_is   <- tab3_chib_mean - tab3_is_ev
diff_inla_chib <- tab3_inla - tab3_chib_mean

df_tab3 <- data.frame(
  Model = tab3_models,
  sigma_beta = tab3_sbeta,
  INLA = round(tab3_inla, 4),
  IS_ref = round(tab3_is_ev, 4),
  IS_mcse = round(tab3_is_se, 4),
  Chib_mean = round(tab3_chib_mean, 4),
  Chib_sd = round(tab3_chib_sd, 4),
  Chib_mcse = round(tab3_chib_se, 4),
  Diff_INLA_IS = round(diff_inla_is, 4),
  Diff_Chib_IS = round(diff_chib_is, 4),
  Diff_INLA_Chib = round(diff_inla_chib, 4)
)
for (r in seq_len(N_REPS)) {
  df_tab3[[paste0("Chib_rep", r)]] <- round(tab3_chib[, r], 4)
}
write.csv(df_tab3, "results/table3_probit.csv", row.names = FALSE)
print(df_tab3[, c("Model", "sigma_beta", "INLA", "IS_ref", "Chib_mean", "Chib_sd", "Diff_INLA_IS", "Diff_INLA_Chib")])

# ------------------------------------------------------------------------------
# 4. TABLE 4: Logistic Regression with Logit Link (Pima Data, Section 5)
# ------------------------------------------------------------------------------
cat("\n--- TABLE 4: Pima Indians Logit Regression (Section 5) ---\n")
pima <- load_pima()

# Run INLA for M1 and M2 at precisions 0.01 (sbeta = 10) and 1.0 (sbeta = 1)
m1_100 <- inla_binom(pima$f1, pima$data, "logit", 0.01)
m2_100 <- inla_binom(pima$f2, pima$data, "logit", 0.01)
m1_1   <- inla_binom(pima$f1, pima$data, "logit", 1.0)
m2_1   <- inla_binom(pima$f2, pima$data, "logit", 1.0)

# Run independent Importance Sampling reference to establish ground truth
cat("  Running Importance Sampling benchmark for Pima models... ")
is_pima_m1_100 <- is_logev(pima$f1, pima$data, "logit", 10, N = n_is_samples, seed = 5001)
is_pima_m2_100 <- is_logev(pima$f2, pima$data, "logit", 10, N = n_is_samples, seed = 5002)
is_pima_m1_1   <- is_logev(pima$f1, pima$data, "logit", 1,  N = n_is_samples, seed = 5003)
is_pima_m2_1   <- is_logev(pima$f2, pima$data, "logit", 1,  N = n_is_samples, seed = 5004)
cat("done.\n")

# Benchmark literature methods from Friel & Wyse (2012)
df_tab4 <- data.frame(
  Method = c("INLA (classic)",
             "Importance sampling (reference)",
             "Laplace approximation",
             "Chib and Jeliazkov",
             "Laplace approximation MAP",
             "Harmonic mean estimator",
             "Power posteriors",
             "Annealed importance sampling",
             "Nested sampling"),
  M1_s100 = c(round(m1_100$mlik[1,1], 2), round(is_pima_m1_100$logev, 2), -257.26, -257.23, -257.28, -279.47, -257.98, -257.87, -258.82),
  M2_s100 = c(round(m2_100$mlik[1,1], 2), round(is_pima_m2_100$logev, 2), -259.89, -259.84, -259.90, -284.78, -260.59, -260.43, -261.38),
  M1_s1   = c(round(m1_1$mlik[1,1], 2),   round(is_pima_m1_1$logev, 2),   -247.33, -247.31, -247.33, -259.84, -247.57, -247.30, -246.82),
  M2_s1   = c(round(m2_1$mlik[1,1], 2),   round(is_pima_m2_1$logev, 2),   -247.59, -247.58, -247.62, -260.55, -247.84, -247.59, -246.97)
)

# Compute Bayes factor log BF_12 = log m(y|M1) - log m(y|M2)
df_tab4$logBF_s100 <- round(df_tab4$M1_s100 - df_tab4$M2_s100, 2)
df_tab4$logBF_s1   <- round(df_tab4$M1_s1   - df_tab4$M2_s1,   2)

write.csv(df_tab4, "results/table4_pima.csv", row.names = FALSE)
print(df_tab4)

# ------------------------------------------------------------------------------
# 5. TABLE 5: Decision-Scale Evaluation (Bayes Factors and Model Selection)
# ------------------------------------------------------------------------------
cat("\n--- TABLE 5: Decision-Scale Evaluation (log BF and Posterior Probability) ---\n")

# Function to build decision-scale comparison row
make_dec_row <- function(dataset, prior_str, mlik_ref_m1, mlik_ref_m2, mlik_inla_m1, mlik_inla_m2) {
  lbf_ref  <- mlik_ref_m1 - mlik_ref_m2
  lbf_inla <- mlik_inla_m1 - mlik_inla_m2
  delta_lbf <- lbf_inla - lbf_ref
  
  p_m1_ref  <- post_m1(lbf_ref)
  p_m1_inla <- post_m1(lbf_inla)
  delta_prob <- p_m1_inla - p_m1_ref
  
  cat_ref  <- kr_class(lbf_ref)
  cat_inla <- kr_class(lbf_inla)
  same_decision <- (cat_ref == cat_inla) && (sign(lbf_ref) == sign(lbf_inla))
  
  data.frame(
    Dataset = dataset,
    Prior = prior_str,
    logBF_Ref = round(lbf_ref, 3),
    logBF_INLA = round(lbf_inla, 3),
    Delta_logBF = round(delta_lbf, 4),
    P_M1_Ref = round(p_m1_ref, 4),
    P_M1_INLA = round(p_m1_inla, 4),
    Delta_P_M1 = round(delta_prob, 4),
    Evidence_Category = cat_ref,
    Decision_Altered = !same_decision
  )
}

dec_rows <- list(
  # Crime (Gaussian)
  make_dec_row("US Crime", "sigma_beta = 1000", tab2_exact[1], tab2_exact[4], tab2_inla1[1], tab2_inla1[4]),
  make_dec_row("US Crime", "sigma_beta = 10",   tab2_exact[2], tab2_exact[5], tab2_inla1[2], tab2_inla1[5]),
  make_dec_row("US Crime", "sigma_beta = 0.1",  tab2_exact[3], tab2_exact[6], tab2_inla1[3], tab2_inla1[6]),
  # Probit (Binary)
  make_dec_row("Sim Probit", "sigma_beta = 1000", tab3_is_ev[1], tab3_is_ev[4], tab3_inla[1], tab3_inla[4]),
  make_dec_row("Sim Probit", "sigma_beta = 10",   tab3_is_ev[2], tab3_is_ev[5], tab3_inla[2], tab3_inla[5]),
  make_dec_row("Sim Probit", "sigma_beta = 0.1",  tab3_is_ev[3], tab3_is_ev[6], tab3_inla[3], tab3_inla[6]),
  # Pima (Logit)
  make_dec_row("Pima Logit", "sigma_beta^2 = 100", is_pima_m1_100$logev, is_pima_m2_100$logev, m1_100$mlik[1,1], m2_100$mlik[1,1]),
  make_dec_row("Pima Logit", "sigma_beta^2 = 1",   is_pima_m1_1$logev,   is_pima_m2_1$logev,   m1_1$mlik[1,1],   m2_1$mlik[1,1])
)

df_dec <- do.call(rbind, dec_rows)
write.csv(df_dec, "results/table5_decision_scale.csv", row.names = FALSE)
print(df_dec)

# Save all seeds recorded during execution
write_seeds("results/seeds.csv")
cat("\n01_tables.R completed successfully. Results saved in results/\n")
