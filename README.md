# Reproducibility Package: Estimating Marginal Likelihoods with INLA

**Manuscript:** *Estimating the marginal likelihood with Integrated nested Laplace approximation (INLA)*  
**Authors:** Aliaksandr Hubin and Geir Storvik  
**Journal:** *Journal of Reproducible Statistics* (FAIR Press Journals)  
**Canonical Repository:** [https://github.com/aliaksah/INLA-marginal-likelihood](https://github.com/aliaksah/INLA-marginal-likelihood)  
**Historical 2016 Archive:** [https://github.com/aliaksah/EMJMCMC2016](https://github.com/aliaksah/EMJMCMC2016)  

---

## 1. Overview and Purpose

This reproducibility package provides the complete, self-contained suite of data, standalone scripts, and automated orchestration workflows required to reproduce every analytical benchmark, numerical simulation, MCMC comparison, and vector PDF figure presented in the manuscript.

The study systematically evaluates the accuracy, numerical stability, and computational efficiency of marginal likelihood approximations produced by **Integrated Nested Laplace Approximations (INLA)** against:
- Closed-form analytical ground truth (Section 2, Table 1; Section 3, Table 2).
- Converged Gibbs sampling with Chib's (1995) method (Section 3, Table 2; Section 4, Table 3; Figures 1--4).
- High-precision independent Importance Sampling references (Section 4, Table 3; Section 5, Table 4).
- Seven benchmark thermodynamic and MCMC evidence estimators from Friel & Wyse (2012) (Section 5, Table 4).
- Longitudinal Poisson GLMM estimators including Chib & Jeliazkov (2001), Chib et al. (1998), and independent Stan Bridge Sampling (Section 6, Table 6).
- Two non-Gaussian and small-cluster stress tests evaluated against 50-point Gauss-Hermite numerical quadrature (Section 8, Table 7).
- Decision-scale Bayes factor and posterior model probability evaluations (Section 7, Table 5).

---

## 2. Software and Environment Specification

All published benchmarks were computed and verified under the following software environment:
- **Operating System:** macOS Darwin (Apple Silicon arm64) and Linux Ubuntu (x86_64).
- **R Version:** R $\ge 4.4.0$ (tested on R 4.4.1).
- **Key Packages:**
  - `INLA` (version 24.06.27 / 24.06.27-1, stable release from [https://inla.r-inla-download.org](https://inla.r-inla-download.org)).
  - `MCMCpack` (version $\ge 1.7.0$, for Gibbs sampling and Chib 1995 marginal likelihoods).
  - `rstan` (version $\ge 2.32.0$) and `bridgesampling` (version $\ge 1.1.2$, for GLMM bridge sampling).
  - `MASS` (version $\ge 7.3.60$, for Pima Indians dataset).
  - `geepack` (version $\ge 1.3.11$, for epileptic seizure counts; offline text fallback included).
  - `statmod` (version $\ge 1.5.0$, for Gauss-Hermite quadrature nodes).
  - `Matrix`, `coda`, `mvtnorm`, `parallel`.

To install all required dependencies from CRAN and the official R-INLA repository:
```R
# Install CRAN packages
install.packages(c("MASS", "geepack", "MCMCpack", "coda", "Matrix", 
                   "mvtnorm", "statmod", "rstan", "bridgesampling"))

# Install R-INLA
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

---

## 3. Quick Start: One-Click Full Reproduction

To run the complete pipeline from start to finish:
```bash
Rscript run_all.R
```

For a rapid continuous-integration verification run (using smaller sample counts and fewer iterations), set the environment variable `QUICK=1`:
```bash
QUICK=1 Rscript run_all.R
```

---

## 4. Script and Artifact Inventory

Every numbered script in this package is modular and directly corresponds to specific tables, figures, and sections of the revised manuscript:

| Script | Manuscript Section | Generated Tables & Figures | Description |
|:---|:---|:---|:---|
| `00_environment.R` | Environment / Sec. 2--6 | `logs/sessionInfo_*.txt`<br>`results/inla_mode_comparison.csv` | Records full system sessionInfo, compiler details, and evaluates INLA classic vs. compact mode across all models. |
| `01_tables.R` | Sections 2, 3, 4, 5, 7 | `results/table1_linear.csv`<br>`results/table2_crime.csv`<br>`results/table3_probit.csv`<br>`results/table4_pima.csv`<br>`results/table5_decision_scale.csv` | Reproduces Tables 1--4 with analytical benchmarks, importance sampling references, replicated MCMC standard errors, and Table 5 (decision-scale Bayes factors). |
| `02_prior_range_sim.R` | Sections 3 & 4 (Figs 1 & 3) | `results/prior_range_summary.csv` | Evaluates INLA vs. Exact / Chib across 100 perturbed prior precisions; computes error distributions (median, mean, 95th percentile, max). |
| `03_convergence_figs.R` | Sections 3 & 4 (Figs 2 & 4) | `results/convergence_summary.csv`<br>`fig2_data.rds`, `fig4_data.rds` | Generates and summarizes MCMC convergence trajectories across iteration counts ($200$ to $102\,400$). |
| `04_section6_glmm.R` | Section 6 | `results/table6_glmm_comparison.csv` | Evaluates the longitudinal Poisson GLMM across INLA strategies (Gaussian, Simplified Laplace, Full Laplace), integrations (CCD, Grid, EB), and Stan bridge sampling. |
| `05_stress_tests.R` | Section 8 | `results/stress_tests.csv` | Probes INLA's limits on low-count Poisson GLMMs (S1) and small-cluster binary GLMMs (S2) against 50-point Gauss-Hermite numerical quadrature. |
| `generate_figures.R` | Figures 1, 2, 3, 4 | `apppic/n100ci1.pdf`<br>`apppic/n100ci2.pdf`<br>`apppic/citune1.pdf`<br>`apppic/l10ci1.pdf`<br>`apppic/l10ci2.pdf`<br>`apppic/lcitune1.pdf` | Renders publication-quality vector PDFs featuring identity reference lines ($y = x$), explicit log-scale labels, and in-panel legends. |
| `run_all.R` | All Sections | All Tables & Figures | Master driver script executing the full suite in sequential order. |

---

## 5. Datasets

All datasets required by the benchmarks are included locally in this directory for offline execution:
- `simcen-x1.txt`, `simcen-y1.txt`: Standardized US crime data covariates and response (Section 3).
- `sim3-X.txt`, `sim3-Y.txt`: Simulated binary data for probit regression (Section 4).
- `epileptic.txt`: Longitudinal epileptic seizure counts from Diggle et al. (1994) / Thall & Vail (1990) (Section 6).
- `Pima.tr`, `Pima.te`: Standardized diabetes dataset from MASS (Section 5).
- `prob1.csv`--`prob5.csv`: Simulation output arrays across 100 perturbed prior precision settings.

---

## 6. RNG Seeds and Deterministic Reproduction

To satisfy reproducibility requirements without stochastic drift, every stochastic estimator is assigned an explicit, deterministic random seed. All seeds are registered during execution and archived in `results/seeds.csv`. Running `run_all.R` reproduces every Monte Carlo replicate column in Tables 1, 2, and 3 identically.

---

## 7. Canonical Repository and Continuous Integration

- **Canonical Repository:** [https://github.com/aliaksah/INLA-marginal-likelihood](https://github.com/aliaksah/INLA-marginal-likelihood) is the official repository for the revised manuscript.
- **Continuous Integration:** Automated multi-platform tests run via GitHub Actions (`.github/workflows/reproduce.yml`) across Ubuntu x86_64 and macOS arm64, verifying cross-platform numerical agreement and logging session info artifacts.

---

## 8. License

- Code: MIT License
- Documentation & Data: Creative Commons Attribution 4.0 International (CC-BY 4.0)
