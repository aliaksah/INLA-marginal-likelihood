# ==============================================================================
# 04_section6_glmm.R: Longitudinal Poisson GLMM Systematic Comparison
# ==============================================================================
source("R/common.R")
suppressPackageStartupMessages({
  library(rstan)
  library(bridgesampling)
})

cat_header("04: Longitudinal Poisson GLMM Systematic INLA & Reference Benchmark")

# Load seizure data
d_glmm <- glmm_inla_data()
cat("Dataset: 58 patients, 290 observations total.\n")

# ------------------------------------------------------------------------------
# 1. Systematic INLA Evaluation across Strategies and Integration Schemes
# ------------------------------------------------------------------------------
cat("\nEvaluating INLA across strategies and integration schemes...\n")

configs <- list(
  list(name = "Default (Simplified Laplace + CCD)", strat = "simplified.laplace", int = "ccd", extra = list()),
  list(name = "Gaussian Approximation + CCD",      strat = "gaussian",           int = "ccd", extra = list()),
  list(name = "Full Laplace + CCD",                strat = "laplace",            int = "ccd", extra = list()),
  list(name = "Simplified Laplace + Grid",         strat = "simplified.laplace", int = "grid", extra = list()),
  list(name = "Full Laplace + Grid",               strat = "laplace",            int = "grid", extra = list()),
  list(name = "Empirical Bayes (EB mode)",         strat = "simplified.laplace", int = "eb",   extra = list()),
  list(name = "Fine Grid (dz = 1.0, Full Laplace)",strat = "laplace",            int = "grid", extra = list(dz = 1.0))
)

inla_res <- list()
for (cfg in configs) {
  cat(sprintf("  Running INLA: %-35s ... ", cfg$name))
  ctrl_inla <- c(list(strategy = cfg$strat, int.strategy = cfg$int), cfg$extra)
  m <- inla(GLMM_FORMULA, data = d_glmm, family = "poisson",
            control.compute = list(mlik = TRUE),
            control.fixed = ctrl_fixed(0.01),
            control.family = list(link = "log"),
            control.inla = ctrl_inla)
  cat(sprintf("done. (Int = %.2f, Gau = %.2f)\n", m$mlik[1, 1], m$mlik[2, 1]))
  inla_res[[length(inla_res) + 1]] <- data.frame(
    Configuration = cfg$name,
    Strategy = cfg$strat,
    Integration = cfg$int,
    LogMargLik_Integration = round(m$mlik[1, 1], 2),
    LogMargLik_Gaussian = round(m$mlik[2, 1], 2),
    stringsAsFactors = FALSE
  )
}
df_inla_glmm <- do.call(rbind, inla_res)

# ------------------------------------------------------------------------------
# 2. Reference Benchmarks: Chib & Jeliazkov (2001), Chib et al. (1998), and Stan
# ------------------------------------------------------------------------------
cat("\nEvaluating Reference Benchmarks for Section 6 GLMM...\n")
bs_mean <- NA; bs_sd <- NA; bs_se <- NA

has_stan <- FALSE
tryCatch({
  rstan_options(auto_write = TRUE)
  stan_model_file <- "stan/glmm_epilepsy.stan"
  cat("  Attempting to compile Stan model... ")
  sm <- stan_model(file = stan_model_file)
  has_stan <- TRUE
  cat("done.\n")
}, error = function(e) {
  cat("\n  [Note: C++ compilation via Stan unavailable on this system (requires Xcode license agreement).]\n")
  cat("  Using published literature benchmarks: Chib & Jeliazkov (2001) and Chib et al. (1998).\n")
})

if (has_stan) {
  stan_dat <- list(
    N   = nrow(d_glmm),
    J   = length(unique(d_glmm$id2)),
    y   = as.integer(d_glmm$y),
    ofs = d_glmm$ofs,
    trt = d_glmm$trt,
    bsl = d_glmm$bsl,
    id  = as.integer(d_glmm$id2)
  )
  n_bridge_reps <- if (QUICK) 2 else 5
  bridge_results <- numeric(n_bridge_reps)
  for (r in seq_len(n_bridge_reps)) {
    seed_r <- 60000 + r
    fit <- sampling(sm, data = stan_dat, iter = if (QUICK) 3000 else 10000,
                    warmup = if (QUICK) 1000 else 2000,
                    chains = 2, cores = min(2, N_CORES), seed = seed_r,
                    control = list(adapt_delta = 0.95), refresh = 0)
    bs <- bridge_sampler(fit, silent = TRUE)
    bridge_results[r] <- bs$logml
    note_seed(sprintf("glmm_bridge_rep%d", r), seed_r)
  }
  bs_mean <- mean(bridge_results)
  bs_sd   <- sd(bridge_results)
  bs_se   <- bs_sd / sqrt(n_bridge_reps)
}

# Literature values from published papers
chib2001_val <- -915.49
chib2001_se  <- 0.12 # Reported Monte Carlo standard error in Chib & Jeliazkov (2001)
chib1998_val <- -915.23

ref_rows <- list(
  data.frame(
    Method_or_Setting = "Chib & Jeliazkov (2001) [JASA]",
    Type = "MCMC (Chib-Jeliazkov)",
    Log_Evidence = chib2001_val,
    Log_Evidence_Gaussian = NA,
    Uncertainty = sprintf("SE = %.2f (10 MCMC reps)", chib2001_se),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Method_or_Setting = "Chib, Greenberg & Winkelmann (1998) [J. Econometrics]",
    Type = "MCMC (Kernel density)",
    Log_Evidence = chib1998_val,
    Log_Evidence_Gaussian = NA,
    Uncertainty = "Reported point estimate",
    stringsAsFactors = FALSE
  )
)
if (!is.na(bs_mean)) {
  ref_rows[[length(ref_rows) + 1]] <- data.frame(
    Method_or_Setting = "Bridge sampling (Stan)",
    Type = "MCMC Bridge",
    Log_Evidence = round(bs_mean, 2),
    Log_Evidence_Gaussian = NA,
    Uncertainty = sprintf("SD = %.2f (SE = %.2f)", bs_sd, bs_se),
    stringsAsFactors = FALSE
  )
}

df_inla_part <- data.frame(
  Method_or_Setting = df_inla_glmm$Configuration,
  Type = "INLA",
  Log_Evidence = df_inla_glmm$LogMargLik_Integration,
  Log_Evidence_Gaussian = df_inla_glmm$LogMargLik_Gaussian,
  Uncertainty = "Deterministic",
  stringsAsFactors = FALSE
)

df_all <- rbind(df_inla_part, do.call(rbind, ref_rows))
write.csv(df_all, "results/table5_glmm_comparison.csv", row.names = FALSE)
write.csv(df_all, "results/table6_glmm_comparison.csv", row.names = FALSE) # alias
write_seeds("results/seeds.csv")
print(df_all)
cat("\n04_section6_glmm.R completed successfully. Results saved to results/table5_glmm_comparison.csv (and alias table6_glmm_comparison.csv)\n")
