# Reproducibility Package: Estimating Marginal Likelihoods with INLA

**Paper:** *Estimating the marginal likelihood with Integrated nested Laplace approximation (INLA)*  
**Authors:** Aliaksandr Hubin and Geir Storvik  
**Journal:** *Journal of Reproducible Statistics* (FAIR Press Journals)

---

## 1. Overview

This self-contained reproducibility package provides the complete code, datasets, and execution pipeline to reproduce all marginal likelihood estimation benchmarks, comparative evaluations, and figures presented in the manuscript.

The study evaluates **Integrated Nested Laplace Approximations (INLA)** against:
- Analytical ground truth (Section 2)
- Harmonic mean estimator (Sections 2, 5)
- Chib's (1995) Gibbs sampling method (Sections 3, 4)
- Chib and Jeliazkov's (2001) Metropolis–Hastings / GLMM method (Sections 5, 6)
- Standard Laplace and Laplace MAP approximations (Section 5)
- Power posteriors, Annealed Importance Sampling, and Nested Sampling (Section 5)

---

## 2. System and Software Requirements

- **Language:** R ($\ge 4.4.0$)
- **Required Packages:**
  - `INLA` ($\ge 24.06.27$)
  - `MCMCpack` ($\ge 1.7.0$)
  - `MASS`
  - `geepack`
  - `coda`

To install the required packages in R:
```R
# Install CRAN packages
install.packages(c("MASS", "geepack", "MCMCpack", "coda"))

# Install R-INLA
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

---

## 3. Quick Start: One-Click Full Reproduction

To reproduce all results from scratch in a single run:

```bash
Rscript run_all.R
```

This master script:
1. Executes all INLA models and exact analytical benchmarks across Sections 2–6.
2. Computes the Monte Carlo baseline methods (all 5 replications of the Harmonic Mean estimator and Chib's MCMC algorithms).
3. Compares all reproduced values against the manuscript's reported values in formatted tables.
4. Generates high-resolution vector PDF figures for Figures 1–4 in `../apppic/`.

---

## 4. Package Structure and Inventory

This package is designed to be fully self-contained and run completely offline without external internet dependencies:

| File | Purpose |
|:---|:---|
| `run_all.R` | Master reproduction script executing all models, baselines, and tables (Sections 2–6). |
| `generate_figures.R` | Generates publication-ready vector PDF figures (Figures 1–4). |
| `simcen-x1.txt`, `simcen-y1.txt` | US crime covariates and response data (Section 3). |
| `sim3-X.txt`, `sim3-Y.txt` | Simulated binary response and covariates for probit regression (Section 4). |
| `epileptic.txt` | Longitudinal epileptic seizure counts from Diggle et al. (1994) (Section 6). |
| `prob1.csv`–`prob5.csv` | Full simulation arrays across 100 perturbed prior precision replications. |
| `fig2_data.rds`, `fig4_data.rds` | Convergence tracking data across MCMC iteration counts. |

*Historical Code Reference:* The original 2016 exploratory and raw draft scripts are still reproducible and archived in the historical GitHub repository: [https://github.com/aliaksah/EMJMCMC2016](https://github.com/aliaksah/EMJMCMC2016).

---

## 5. Expected Numerical Precision Across Architectures

When reproducing the results across different operating systems, CPU architectures, and R-INLA versions, readers and reviewers should expect the following numerical precision:

1. **Exact Analytical Models (Table 1):**
   - Matches to **4 decimal places** (`-7.8267`, `-3.2463`, `-2.9041`) across all platforms.
2. **Standardized GLMs (Table 4, Pima Indians):**
   - Matches to **2 decimal places** (`-257.25`, `-259.89`, `-247.32`, `-247.59`) across all platforms.
3. **Multi-dimensional Numerical Laplace Integrations (Tables 2 & 3):**
   - Agrees within **$\Delta \le 0.003$ log units** (e.g., Table 3 matches within $\Delta = 0.0002$, Table 2 matches within $\Delta = 0.0033$).
   - Minor variations in the 3rd or 4th decimal place are expected and normal when moving between legacy **Intel x86_64** hardware (used in 2016) and modern **ARM64** hardware (e.g., Apple Silicon M-series), as well as across upstream releases of the core GMRFLib solver (quadrature weights, sparse Cholesky tolerances, and compiler vector math).
4. **Longitudinal GLMM (Section 6):**
   - Gaussian log-marginal likelihood reproduces at $-915.67$ (paper reported $-915.61$, Chib reported $-915.49 \pm 0.12$).

---

## 6. Prior Specifications and Errata

- **Section 5 (Table 4, Pima Indians logit model):** In the original 2016 arXiv preprint text, the prior mean $\mu_\beta$ was inadvertently listed as $1$. As verified against the code implementation and the benchmark study of Friel & Wyse (2012), the true prior mean applied was $\mu_\beta = 0$ with precisions $\sigma_\beta^{-2} \in \{0.01, 1\}$. The manuscript text has been corrected.
- **Cross-Platform INLA Stability:** On modern ARM architectures (such as Apple Silicon macOS), `inla.setOption("inla.mode", "classic")` and `inla.setOption("num.threads", "1:1")` are set in `run_all.R` to ensure cross-platform stability and numerical determinism.

---

## 7. License

- Code: MIT License
- Documentation & Data: Creative Commons Attribution 4.0 International (CC-BY 4.0)
