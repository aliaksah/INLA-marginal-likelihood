# ==============================================================================
# Shared helpers for the reproducibility package
# "Estimating the marginal likelihood with INLA" (Hubin & Storvik)
# ==============================================================================
# Sourced by every numbered script. Provides: package checks, INLA settings,
# data loaders, exact / reference evidence calculators, seed bookkeeping and
# decision-scale (Bayes factor) summaries.
# ==============================================================================

# ---- run mode ---------------------------------------------------------------
# QUICK=1 (environment variable) shrinks every Monte Carlo job so that the
# whole pipeline finishes in minutes (used by the continuous-integration
# workflow). The published numbers are produced with QUICK unset.
QUICK  <- identical(Sys.getenv("QUICK"), "1")
N_REPS <- if (QUICK) 2 else 5        # replications of every stochastic estimator
N_CORES <- max(1, min(8, parallel::detectCores() - 1))
if (QUICK) N_CORES <- min(N_CORES, 2)

# ---- directories ------------------------------------------------------------
if (dir.exists("reproducibility package")) setwd("reproducibility package")
dir.create("results", showWarnings = FALSE)
dir.create("logs",    showWarnings = FALSE)
APPPIC <- if (dir.exists("../apppic")) "../apppic" else "apppic"
dir.create(APPPIC, showWarnings = FALSE, recursive = TRUE)

# ---- packages ---------------------------------------------------------------
ensure_packages <- function(cran = c("MASS", "geepack", "MCMCpack", "coda", "Matrix", "mvtnorm")) {
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    message("Installing R-INLA from the official repository ...")
    install.packages("INLA", dep = TRUE,
                     repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
  }
}
ensure_packages()
suppressPackageStartupMessages({
  library(INLA); library(MASS); library(MCMCpack); library(coda); library(parallel)
})

# INLA settings used for ALL published numbers (see README):
#  * inla.mode = "classic": the integration-based marginal likelihood of the
#    R-INLA 24.06.27 default ("compact") is not used because it (i) differs
#    by about +0.38 in the Gaussian example (see 00_environment.R) and (ii)
#    aborts for the probit model on arm64 macOS.
#  * one thread, so that results do not depend on the number of cores.
init_inla <- function() {
  inla.setOption(inla.mode = "classic")
  inla.setOption(num.threads = "1:1")
  invisible(NULL)
}
init_inla()

# ---- seed bookkeeping -------------------------------------------------------
.SEEDS <- new.env()
.SEEDS$tab <- data.frame(task = character(), seed = integer(), stringsAsFactors = FALSE)
use_seed <- function(task, seed) {              # set.seed + record
  .SEEDS$tab[nrow(.SEEDS$tab) + 1, ] <- list(task, as.integer(seed))
  set.seed(seed)
  invisible(seed)
}
note_seed <- function(task, seed) {             # record only (MCMCpack seed= argument)
  .SEEDS$tab[nrow(.SEEDS$tab) + 1, ] <- list(task, as.integer(seed))
  invisible(seed)
}
write_seeds <- function(file) {
  old <- if (file.exists(file)) read.csv(file, stringsAsFactors = FALSE) else NULL
  new <- unique(rbind(old, .SEEDS$tab))
  write.csv(new, file, row.names = FALSE)
}

# ---- data -------------------------------------------------------------------
load_crime <- function() {
  simx <- read.table("simcen-x1.txt", sep = ","); simy <- read.table("simcen-y1.txt")
  d <- cbind(simy, simx); names(d)[1] <- "Y"
  list(data = d,
       f1 = as.formula(paste("Y ~ 1 +", paste0(colnames(d)[-c(1,3,5,7,9,11,13,15)], collapse = "+"))),
       f2 = as.formula(paste("Y ~ 1 +", paste0(colnames(d)[-c(1,3,7,9,13)],        collapse = "+"))))
}
load_probit <- function() {
  simx <- read.table("sim3-X.txt", sep = ","); simy <- read.table("sim3-Y.txt", sep = ",")
  d <- cbind(simy, simx); names(d)[1] <- "Y1"
  d$V2 <- (d$V10 + d$V14) * d$V9
  d$V5 <- (d$V11 + d$V15) * d$V12
  list(data = d,
       f1 = as.formula(paste("Y1 ~ 1 +", paste0(colnames(d)[-c(1,3,5,7,9,11,13,15,17,19)], collapse = "+"))),
       f2 = as.formula(paste("Y1 ~ 1 +", paste0(colnames(d)[-c(1,3,7,9,13,15,17,19)],       collapse = "+"))))
}
load_pima <- function() {
  d <- rbind(MASS::Pima.tr, MASS::Pima.te)
  d$type <- as.integer(d$type == "Yes")
  for (i in 1:7) d[, i] <- (d[, i] - mean(d[, i])) / sqrt(var(d[, i]))
  list(data = d,
       f1 = type ~ 1 + npreg + glu + bmi + ped,
       f2 = type ~ 1 + npreg + glu + bmi + ped + age)
}
# Epileptic seizure data (Thall & Vail 1990). The patient with baseline count 151
# (id 207, counts 102/65/72/63) is the well-known outlier that is removed.
# The fallback built from epileptic.txt is identical to geepack::seizure (checked
# in 00_environment.R), so the result does not depend on geepack being installed.
load_seizure <- function() {
  if (requireNamespace("geepack", quietly = TRUE)) {
    e <- new.env(); utils::data("seizure", package = "geepack", envir = e); s <- e$seizure
    s <- s[, c("base", "y1", "y2", "y3", "y4", "trt")]
  } else {
    raw <- read.table("epileptic.txt", sep = ";", header = TRUE)
    w <- reshape(raw[, c("id", "visit", "seizure")], timevar = "visit", idvar = "id", direction = "wide")
    s <- data.frame(base = w$seizure.0, y1 = w$seizure.1, y2 = w$seizure.2,
                    y3 = w$seizure.3, y4 = w$seizure.4,
                    trt = raw$trt[match(w$id, raw$id)])
  }
  stopifnot(nrow(s) == 59, sum(s$base == 151) == 1)
  s <- s[s$base != 151, ]
  s$id <- seq_len(nrow(s))
  rownames(s) <- NULL
  s
}
# long-format data used by R-INLA (same construction as in the original scripts)
glmm_inla_data <- function(s = load_seizure()) {
  n <- nrow(s)
  y   <- c(s$base, s$y1, s$y2, s$y3, s$y4)
  trt <- rep(s$trt, 5)
  bsl <- c(rep(0, n), rep(1, 4 * n))
  ofs <- c(rep(8, n), rep(2, 4 * n))
  id  <- rep(s$id, 5)
  data.frame(y = y, trt = trt, bsl = bsl, ofs = ofs, id2 = id, id3 = id + max(s$id), bsl2 = bsl)
}
GLMM_FORMULA <- y ~ offset(ofs) + I(trt) + I(trt * bsl) + I(bsl) +
  f(id2, model = "iid2d", n = 2 * 290) + f(id3, bsl, copy = "id2")

# ---- INLA wrappers ----------------------------------------------------------
ctrl_fixed <- function(prec) list(mean = c(default = 0), mean.intercept = 0,
                                  prec.intercept = prec, prec = c(default = prec))
inla_gauss <- function(f, d, prec, ...) inla(f, data = d, family = "gaussian",
  control.compute = list(mlik = TRUE),
  control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
  control.fixed = ctrl_fixed(prec), ...)
inla_binom <- function(f, d, link, prec, ...) inla(f, data = d, family = "binomial",
  Ntrials = rep(1, nrow(d)), control.compute = list(mlik = TRUE),
  control.family = list(link = link), control.fixed = ctrl_fixed(prec), ...)
# Row 1 of $mlik = "log marginal-likelihood (integration)", row 2 = "(Gaussian)"
mlik_int <- function(m) unname(m$mlik[1, 1])
mlik_gau <- function(m) unname(m$mlik[2, 1])

# ---- exact evidence: Gaussian regression (Section 3) ------------------------
# beta ~ N(0, sb^2 I) (incl. intercept) is integrated analytically; tau = 1/sigma^2 ~
# Gamma(1, 1) by trapezoidal quadrature on log tau (spectrally accurate, checked
# against a 3x finer grid in 01_tables.R).
exact_gauss <- function(formula, data, sb, npts = 20001, lim = c(-20, 15)) {
  X <- model.matrix(formula, data); y <- data[[all.vars(formula)[1]]]
  n <- length(y); p <- ncol(X)
  XtX <- crossprod(X); Xty <- crossprod(X, y); yty <- sum(y^2); s2 <- sb^2
  lg <- function(u) {
    tau <- exp(u); A <- diag(p) + tau * s2 * XtX
    ldA <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
    q <- tau * (yty - tau * s2 * as.numeric(crossprod(Xty, solve(A, Xty))))
    -0.5 * n * log(2 * pi) + 0.5 * n * u - 0.5 * ldA - 0.5 * q + dgamma(tau, 1, 1, log = TRUE) + u
  }
  u <- seq(lim[1], lim[2], length.out = npts); v <- vapply(u, lg, 0)
  m <- max(v); m + log(sum(exp(v - m)) * (u[2] - u[1]))
}

# ---- Chib (1995) via MCMCpack ----------------------------------------------
# MCMCregress/MCMCprobit with marginal.likelihood = "Chib95"; the prior is
# beta ~ N(0, precp^-1 I) (incl. intercept) and, for the Gaussian model,
# sigma^2 ~ IG(c0/2, d0/2) = IG(1, 1)  <=>  tau ~ Gamma(1, 1).
chib_gauss <- function(f, d, precp, seed, mcmc = 1e5, b0 = 0) {
  mc <- MCMCregress(f, data = d, seed = seed, burnin = 100, b0 = b0, B0 = precp, c0 = 2, d0 = 2,
                    marginal.likelihood = "Chib95", mcmc = mcmc, verbose = 0)
  BayesFactor(mc)$BF.logmarglike
}
chib_probit <- function(f, d, precp, seed, mcmc = 1e5, b0 = 0) {
  mc <- MCMCprobit(f, data = d, seed = seed, burnin = 1000, b0 = b0, B0 = precp,
                   marginal.likelihood = "Chib95", mcmc = mcmc, verbose = 0)
  BayesFactor(mc)$BF.logmarglike
}

# ---- Importance-sampling reference for binary-response GLMs -----------------
# Proposal: multivariate t_6 centred at the posterior mode with the (inflated)
# inverse Hessian. Returns log evidence, Monte Carlo SE and the relative ESS.
# Prior: beta ~ N(0, sb^2 I) including the intercept. Constants are included.
is_logev <- function(f, d, link, sb, N = 2e5, df = 6, seed, batch = 2000) {
  X <- model.matrix(f, d); y <- d[[all.vars(f)[1]]]; p <- ncol(X)
  ll <- function(eta) if (link == "probit")
    sum(y * pnorm(eta, log.p = TRUE) + (1 - y) * pnorm(-eta, log.p = TRUE)) else
    sum(y * plogis(eta, log.p = TRUE) + (1 - y) * plogis(-eta, log.p = TRUE))
  lpr <- function(b) sum(dnorm(b, 0, sb, log = TRUE))
  lp1 <- function(b) ll(drop(X %*% b)) + lpr(b)
  op <- optim(rep(0, p), function(b) -lp1(b), method = "BFGS", hessian = TRUE,
              control = list(reltol = 1e-14, maxit = 2000))
  mode <- op$par; S <- solve(op$hessian) * 1.3; L <- t(chol(S))
  set.seed(seed)
  lw <- numeric(0)
  for (b in seq_len(ceiling(N / batch))) {
    Z <- matrix(rnorm(batch * p), p, batch)
    u <- rchisq(batch, df)
    B <- mode + L %*% sweep(Z, 2, sqrt(u / df), "/")          # p x batch
    eta <- X %*% B                                            # n x batch
    lik <- if (link == "probit")
      colSums(y * pnorm(eta, log.p = TRUE) + (1 - y) * pnorm(-eta, log.p = TRUE)) else
      colSums(y * plogis(eta, log.p = TRUE) + (1 - y) * plogis(-eta, log.p = TRUE))
    lprior <- colSums(matrix(dnorm(B, 0, sb, log = TRUE), p))
    lq <- mvtnorm::dmvt(t(B), delta = mode, sigma = S, df = df, log = TRUE)
    lw <- c(lw, lik + lprior - lq)
  }
  m <- max(lw); w <- exp(lw - m); Nn <- length(w)
  list(logev = m + log(mean(w)),
       se    = sd(w) / (mean(w) * sqrt(Nn)),
       ess   = sum(w)^2 / sum(w^2) / Nn)
}

# ---- decision-scale summaries -----------------------------------------------
# Kass & Raftery (1995) categories for the strength of evidence in ln BF:
# <1 not worth more than a bare mention, 1-3 positive, 3-5 strong, >5 very strong.
kr_class <- function(lnbf) as.character(cut(abs(lnbf), c(-Inf, 1, 3, 5, Inf),
  labels = c("not worth a mention", "positive", "strong", "very strong"), right = FALSE))
# posterior probability of M1 for two equally likely models
post_m1 <- function(lnbf) plogis(lnbf)

rep_summary <- function(x) c(mean = mean(x), sd = sd(x), se = sd(x) / sqrt(length(x)))

cat_header <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
