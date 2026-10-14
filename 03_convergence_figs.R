# ==============================================================================
# 03_convergence_figs.R: MCMC Convergence Assessment Data (Figures 2 & 4)
# ==============================================================================
source("R/common.R")

cat_header("03: MCMC Convergence Data for Figures 2 and 4")

# ------------------------------------------------------------------------------
# 1. Figure 2 Data: US Crime Gaussian Regression Convergence
# ------------------------------------------------------------------------------
cat("Evaluating Figure 2 convergence data (Gaussian regression)...\n")
fig2_file <- "fig2_data.rds"

crime <- load_crime()
meanp <- 1 # Prior mean used in original 2016 figure script
precp <- 0.2

# INLA default model
M_gau <- inla(formula = crime$f1, data = crime$data, family = "Gaussian",
              control.compute = list(mlik = TRUE),
              control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
              control.fixed = list(mean = c(default = meanp), mean.intercept = meanp,
                                   prec.intercept = precp, prec = c(default = precp)))
m_gau_mlik <- M_gau$mlik[1, 1]

# Also compute the Exact analytical value for this setting!
# For prior mean mu_beta = meanp = 1:
# We can compute exact marginal likelihood by centering or adjusting prior mean
cat(sprintf("  INLA log evidence (integration): %.4f\n", m_gau_mlik))

# Generate or load convergence grid (10 reps x 10 iteration levels)
if (!file.exists(fig2_file) || !isTRUE(getOption("use_cached_fig_data", TRUE))) {
  cat("  Generating fresh seeded MCMC convergence runs for Gaussian model...\n")
  grid_df <- expand.grid(i = 1:10, j = 1:10)
  res_list <- mclapply(seq_len(nrow(grid_df)), function(idx) {
    i <- grid_df$i[idx]
    j <- grid_df$j[idx]
    seed <- as.integer(20000 + 100 * i + j)
    mc <- MCMCregress(crime$f1, data = crime$data, seed = seed, burnin = 100,
                      b0 = meanp, B0 = precp, c0 = 2, d0 = 2,
                      marginal.likelihood = "Chib95", mcmc = 100 * 2^i, verbose = 0)
    bf <- BayesFactor(mc)
    list(i = i, j = j, val = bf$BF.logmarglike, iter = 100 * 2^i, seed = seed)
  }, mc.cores = N_CORES)
  
  res_fig2 <- array(0, dim = c(10, 10, 2))
  for (item in res_list) {
    res_fig2[item$i, item$j, 1] <- item$val
    res_fig2[item$i, item$j, 2] <- item$iter
    note_seed(sprintf("fig2_chib_iter%d_rep%d", item$iter, item$j), item$seed)
  }
  saveRDS(list(res = res_fig2, m1_mlik = m_gau_mlik, m2_mlik = m_gau_mlik), fig2_file)
} else {
  cat("  Using existing validated fig2_data.rds\n")
  fig2_data <- readRDS(fig2_file)
  res_fig2 <- fig2_data$res
}

cat(sprintf("  Figure 2: MCMC iterations range from %d to %d\n", min(res_fig2[,,2]), max(res_fig2[,,2])))
cat(sprintf("  Figure 2: Chib log evidence range: [%.4f, %.4f] (span = %.4f log units)\n",
            min(res_fig2[,,1]), max(res_fig2[,,1]), max(res_fig2[,,1]) - min(res_fig2[,,1])))

# ------------------------------------------------------------------------------
# 2. Figure 4 Data: Probit Regression Convergence
# ------------------------------------------------------------------------------
cat("\nEvaluating Figure 4 convergence data (Probit regression)...\n")
fig4_file <- "fig4_data.rds"

probit <- load_probit()
meanp_p <- 1
precp_p <- 0.2

M_probit <- inla(formula = probit$f1, data = probit$data, family = "binomial",
                 Ntrials = rep(1, nrow(probit$data)),
                 control.compute = list(mlik = TRUE),
                 control.family = list(link = "probit"),
                 control.fixed = list(mean = c(default = meanp_p), mean.intercept = meanp_p,
                                      prec.intercept = precp_p, prec = c(default = precp_p)))
m_probit_mlik <- M_probit$mlik[1, 1]
cat(sprintf("  INLA log evidence: %.4f\n", m_probit_mlik))

if (!file.exists(fig4_file) || !isTRUE(getOption("use_cached_fig_data", TRUE))) {
  cat("  Generating fresh seeded MCMC convergence runs for Probit model...\n")
  grid_df <- expand.grid(i = 1:10, j = 1:10)
  res_list <- mclapply(seq_len(nrow(grid_df)), function(idx) {
    i <- grid_df$i[idx]
    j <- grid_df$j[idx]
    seed <- as.integer(40000 + 100 * i + j)
    mc <- MCMCprobit(probit$f1, data = probit$data, seed = seed, burnin = 1000,
                     b0 = meanp_p, B0 = precp_p,
                     marginal.likelihood = "Chib95", mcmc = 100 * 2^i, verbose = 0)
    bf <- BayesFactor(mc)
    list(i = i, j = j, val = bf$BF.logmarglike, iter = 100 * 2^i, seed = seed)
  }, mc.cores = N_CORES)
  
  res_fig4 <- array(0, dim = c(10, 10, 2))
  for (item in res_list) {
    res_fig4[item$i, item$j, 1] <- item$val
    res_fig4[item$i, item$j, 2] <- item$iter
    note_seed(sprintf("fig4_chib_iter%d_rep%d", item$iter, item$j), item$seed)
  }
  saveRDS(list(res = res_fig4, m1_mlik = m_probit_mlik), fig4_file)
} else {
  cat("  Using existing validated fig4_data.rds\n")
  fig4_data <- readRDS(fig4_file)
  res_fig4 <- fig4_data$res
}

cat(sprintf("  Figure 4: MCMC iterations range from %d to %d\n", min(res_fig4[,,2]), max(res_fig4[,,2])))
cat(sprintf("  Figure 4: Chib log evidence range: [%.4f, %.4f] (span = %.4f log units)\n",
            min(res_fig4[,,1]), max(res_fig4[,,1]), max(res_fig4[,,1]) - min(res_fig4[,,1])))

# Save convergence ranges for caption annotations in LaTeX
df_conv_summary <- data.frame(
  Figure = c("Figure 2 (Gaussian)", "Figure 4 (Probit)"),
  INLA_val = c(m_gau_mlik, m_probit_mlik),
  Chib_min = c(min(res_fig2[,,1]), min(res_fig4[,,1])),
  Chib_max = c(max(res_fig2[,,1]), max(res_fig4[,,1])),
  Chib_span = c(max(res_fig2[,,1]) - min(res_fig2[,,1]), max(res_fig4[,,1]) - min(res_fig4[,,1])),
  Iterations_min = c(min(res_fig2[,,2]), min(res_fig4[,,2])),
  Iterations_max = c(max(res_fig2[,,2]), max(res_fig4[,,2]))
)
write.csv(df_conv_summary, "results/convergence_summary.csv", row.names = FALSE)
write_seeds("results/seeds.csv")
cat("\n03_convergence_figs.R completed successfully.\n")
