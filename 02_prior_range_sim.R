# ==============================================================================
# 02_prior_range_sim.R: Prior Range Numerical Evaluation (Figures 1 & 3)
# ==============================================================================
source("R/common.R")

cat_header("02: Prior-Range Numerical Error Evaluation (Figures 1 and 3)")

# ------------------------------------------------------------------------------
# 1. Figure 1: US Crime Data (Gaussian Regression, Prior Precision in [0, 100])
# ------------------------------------------------------------------------------
cat("Evaluating Figure 1 (US Crime Gaussian Regression across 100 prior precisions)...\n")
p5 <- read.csv("prob5.csv")
# Columns: V1 = Chib M1, V2 = Chib M2, V3 = INLA M1, V4 = INLA M2, V5 = Precision

# Compute Exact analytical marginal likelihood across all 100 points for M1 and M2
crime <- load_crime()
cat("  Computing Exact analytical benchmark for all 100 prior precision values... ")
exact_m1 <- numeric(nrow(p5))
exact_m2 <- numeric(nrow(p5))

exact_vals <- mclapply(seq_len(nrow(p5)), function(i) {
  prec <- p5$V5[i]
  sb <- 1 / sqrt(prec)
  c(exact_gauss(crime$f1, crime$data, sb),
    exact_gauss(crime$f2, crime$data, sb))
}, mc.cores = N_CORES)

for (i in seq_len(nrow(p5))) {
  exact_m1[i] <- exact_vals[[i]][1]
  exact_m2[i] <- exact_vals[[i]][2]
}
cat("done.\n")

# Errors against Exact benchmark
err_inla_m1 <- abs(p5$V3 - exact_m1)
err_inla_m2 <- abs(p5$V4 - exact_m2)
err_chib_m1 <- abs(p5$V1 - exact_m1)
err_chib_m2 <- abs(p5$V2 - exact_m2)

# Discrepancy between INLA and Chib
diff_inla_chib_m1 <- abs(p5$V3 - p5$V1)
diff_inla_chib_m2 <- abs(p5$V4 - p5$V2)

# Log Bayes Factor comparison (M1 vs M2)
lbf_exact <- exact_m1 - exact_m2
lbf_inla  <- p5$V3 - p5$V4
lbf_chib  <- p5$V1 - p5$V2
diff_lbf_inla_exact <- abs(lbf_inla - lbf_exact)
diff_lbf_chib_exact <- abs(lbf_chib - lbf_exact)

summarize_dist <- function(x) {
  c(Median = median(x),
    Mean   = mean(x),
    P90    = quantile(x, 0.90, names = FALSE),
    P95    = quantile(x, 0.95, names = FALSE),
    Max    = max(x))
}

cat("\nSummary of Absolute Discrepancy: INLA vs. Exact (Gaussian, M1):\n")
print(round(summarize_dist(err_inla_m1), 4))
cat("Summary of Absolute Discrepancy: INLA vs. Exact (Gaussian, M2):\n")
print(round(summarize_dist(err_inla_m2), 4))
cat("Summary of Absolute Difference in log BF (M1 vs M2, INLA vs Exact):\n")
print(round(summarize_dist(diff_lbf_inla_exact), 4))

# ------------------------------------------------------------------------------
# 2. Figure 3: Simulated Probit Data (Prior Precision in [0, 10])
# ------------------------------------------------------------------------------
cat("\nEvaluating Figure 3 (Simulated Probit Regression across 100 prior precisions)...\n")
p2 <- read.csv("prob2.csv")
# Columns: V1 = Chib M1, V2 = Chib M2, V3 = INLA M1, V4 = INLA M2, V5 = Precision

diff_probit_m1 <- abs(p2$V3 - p2$V1)
diff_probit_m2 <- abs(p2$V4 - p2$V2)

lbf_probit_inla <- p2$V3 - p2$V4
lbf_probit_chib <- p2$V1 - p2$V2
diff_lbf_probit <- abs(lbf_probit_inla - lbf_probit_chib)

cat("Summary of Absolute Discrepancy: INLA vs. Chib (Probit, M1):\n")
print(round(summarize_dist(diff_probit_m1), 4))
cat("Summary of Absolute Discrepancy: INLA vs. Chib (Probit, M2):\n")
print(round(summarize_dist(diff_probit_m2), 4))
cat("Summary of Absolute Difference in log BF (M1 vs M2, INLA vs Chib):\n")
print(round(summarize_dist(diff_lbf_probit), 4))

# Combine into a structured data frame and save
df_summary <- data.frame(
  Setting = c("Gaussian M1 (|INLA - Exact|)",
              "Gaussian M2 (|INLA - Exact|)",
              "Gaussian M1 (|INLA - Chib|)",
              "Gaussian M2 (|INLA - Chib|)",
              "Gaussian log BF (|INLA - Exact|)",
              "Probit M1 (|INLA - Chib|)",
              "Probit M2 (|INLA - Chib|)",
              "Probit log BF (|INLA - Chib|)"),
  rbind(
    summarize_dist(err_inla_m1),
    summarize_dist(err_inla_m2),
    summarize_dist(diff_inla_chib_m1),
    summarize_dist(diff_inla_chib_m2),
    summarize_dist(diff_lbf_inla_exact),
    summarize_dist(diff_probit_m1),
    summarize_dist(diff_probit_m2),
    summarize_dist(diff_lbf_probit)
  )
)
write.csv(df_summary, "results/prior_range_summary.csv", row.names = FALSE)
cat("\nPrior range summary saved to results/prior_range_summary.csv\n")
