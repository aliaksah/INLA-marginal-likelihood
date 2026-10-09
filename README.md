# Reproducibility Package: Estimating Marginal Likelihoods with INLA

**Paper:** *Estimating the marginal likelihood with Integrated nested Laplace approximation (INLA)*  
**Authors:** Aliaksandr Hubin and Geir Storvik  


---

## 1. Overview

This package contains the code, data, and scripts required to reproduce all marginal likelihood estimation results, simulation benchmarks, and comparative evaluations presented in the manuscript.

The study compares the approximations of the marginal likelihood obtained with **Integrated Nested Laplace Approximations (INLA)** against:
- Analytical ground truth (Section 2)
- Harmonic mean estimator (Sections 2, 5)
- Chib's (1995) Gibbs-based method (Sections 3, 4)
- Chib and Jeliazkov's (2001) Metropolis–Hastings / GLMM method (Sections 5, 6)
- Standard Laplace and Laplace MAP approximations (Section 5)
- Power posteriors (Section 5)
- Annealed importance sampling (Section 5)
- Nested sampling (Section 5)

---

## 2. System and Software Requirements

- **Language:** R (tested on R $\ge$ 4.4.0)
- **Key R Packages:**
  - `INLA` (tested on version 24.06.27 and 26.08.07)
  - `MCMCpack` (for Chib's marginal likelihood via MCMCregress / MCMCprobit)
  - `MASS` (for Pima Indians data and distribution utilities)
  - `geepack` (for the epileptic seizure count data in Section 6)
  - `coda` (MCMC diagnostics and output analysis)
  - `RCurl` (optional, for remote downloads if local files are absent)

To install the required CRAN packages in R:
```R
install.packages(c("MASS", "geepack", "MCMCpack", "coda", "RCurl"))
```
To install R-INLA:
```R
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

---

## 3. Quick Start: One-Click Full Reproduction

To run the full suite reproducing all INLA results and baseline methods (Harmonic Mean, Chib's MCMC, and literature benchmarks) for Tables 1, 2, 3, 4, Section 6, and Figures 1–4 in one command:

```bash
Rscript run_all.R
```

This master script executes the models, computes the baseline replications, generates high-resolution vector PDF figures, and prints aligned comparison tables matching the manuscript.

---

## 4. Script and Artifact Inventory

| Script | Paper Section | Output / Artifact | Description |
|:---|:---|:---|:---|
| `run_all.R` | All Sections | Tables 1, 2, 3, 4, Sec 6 | Master reproduction script for all INLA calculations |
| `example1_mlikers.R` | Section 2 | Table 1 | Simple Gaussian latent model: INLA vs. Exact vs. Harmonic Mean |
| `example2final_mlikers.R` | Section 3 | Table 2, Figures 1 & 2 | Gaussian Bayesian regression (US crime data): INLA vs. Chib (1995) |
| `example3final_mlikers.R` | Section 4 | Table 3, Figures 3 & 4 | Logistic regression with probit link: INLA vs. Chib (1995) |
| `example_logit_mlikers.R` | Section 5 | Table 4 | Logistic regression with logit link (Pima data): INLA vs. Friel & Wyse benchmarks |
| `mixed_model.R` | Section 6 | Section 6 | Poisson GLMM with longitudinal random effects (Epileptic seizures): INLA vs. Chib & Jeliazkov (2001) |
| `compare mlikers.R` | Sections 3–4 | `prob1.csv`–`prob5.csv` | Full simulation studies evaluating INLA and Chib under perturbed prior precisions |
| `generate_figures.R` | Sections 3–4 | Figures 1, 2, 3, 4 (`apppic/*.pdf`) | Generates vector PDF figures for the manuscript |

---

## 5. Datasets

All datasets required by the scripts are included locally within this directory:
- `simcen-x1.txt`, `simcen-y1.txt`: Standardized US crime data covariates and response (Section 3).
- `sim3-X.txt`, `sim3-Y.txt`: Simulated binary data for probit logistic regression (Section 4).
- `epileptic.txt`: Longitudinal epileptic seizure counts from Diggle et al. (1994) / Chib and Jeliazkov (2001) (Section 6).
- `Pima.tr`, `Pima.te`: Built into the R `MASS` package (Section 5).
- `prob1.csv`–`prob5.csv`: Simulation output arrays across 100 perturbed prior precision replications.

---

## 6. Notes on Prior Specifications and Errata

- **Section 5 (Table 4, Pima Indians logit model):** In the original arXiv preprint text, the prior mean $\mu_\beta$ was inadvertently written as 1. As noted in editorial checks and verified against the implementation code (`example_logit_mlikers.R`) and the benchmark study of Friel & Wyse (2012), the true prior mean applied was $\mu_\beta = 0$ with precisions $\sigma_\beta^{-2} \in \{0.01, 1\}$. The manuscript text has been updated accordingly.
- **Cross-Platform INLA Stability:** For probit regression with large sample sizes on modern ARM architectures (e.g., Apple Silicon macOS), `inla.setOption("inla.mode", "classic")` and `inla.setOption("num.threads", "1:1")` are recommended to guarantee identical numerical convergence.

---

## 7. License

- Code: MIT License
- Manuscript & Documentation: Creative Commons Attribution 4.0 International (CC-BY 4.0)
