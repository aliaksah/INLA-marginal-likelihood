# ==============================================================================
# Master Reproduction Script for:
# "Estimating the marginal likelihood with Integrated nested Laplace approximation (INLA)"
# Authors: Aliaksandr Hubin & Geir Storvik
# Journal of Reproducible Statistics
# ==============================================================================
# This master script sequentially executes all modular benchmarks:
# - 00_environment.R:     Logs environment details, compares classic vs compact modes
# - 01_tables.R:          Reproduces Tables 1-4 with exact/IS benchmarks & Table 5
# - 02_prior_range_sim.R: Evaluates numerical error across 100 prior precisions
# - 03_convergence_figs.R:Computes MCMC convergence data for Figures 2 & 4
# - 04_section6_glmm.R:   Evaluates longitudinal Poisson GLMM (INLA vs Stan Bridge)
# - 05_stress_tests.R:    Runs S1 & S2 stress tests against Gauss-Hermite references
# - generate_figures.R:   Renders all publication-ready vector PDF figures
# ==============================================================================

rm(list = ls(all = TRUE))

if (dir.exists("reproducibility package")) setwd("reproducibility package")

cat("==============================================================================\n")
cat("Executing Master Reproduction Suite: INLA Marginal Likelihood (Hubin & Storvik)\n")
cat("R version:    ", R.version.string, "\n")
cat("Platform:     ", R.version$platform, "\n")
cat("QUICK mode:   ", if (identical(Sys.getenv("QUICK"), "1")) "ENABLED (Fast CI run)" else "DISABLED (Full publication run)", "\n")
cat("==============================================================================\n\n")

# Step 0: Software Environment and Mode Assessment
source("00_environment.R")

# Step 1: Core Comparison Tables (Tables 1-4) and Decision Scale (Table 5)
source("01_tables.R")

# Step 2: Prior Range Simulation Error Quantification (Figures 1 and 3)
source("02_prior_range_sim.R")

# Step 3: MCMC Convergence Assessment Data (Figures 2 and 4)
source("03_convergence_figs.R")

# Step 4: Longitudinal Poisson GLMM Benchmark (Section 6 & Table 6)
source("04_section6_glmm.R")

# Step 5: Stress Testing Non-Gaussian and Small-Cluster Limits (Table 7)
source("05_stress_tests.R")

# Step 6: Vector PDF Figure Generation (apppic/)
source("generate_figures.R")

cat("\n==============================================================================\n")
cat("SUCCESS: All models, comparisons, stress tests, and figures executed!\n")
cat("Results are archived in results/, logs in logs/, and figures in apppic/.\n")
cat("==============================================================================\n")
