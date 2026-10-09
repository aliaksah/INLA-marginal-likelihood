# ==============================================================================
# Generate High-Resolution Vector PDF Figures for JRS Manuscript
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(MCMCpack)
  library(INLA)
})

inla.setOption("inla.mode", "classic")
inla.setOption("num.threads", "1:1")

output_dir <- if (dir.exists("../apppic")) "../apppic" else "apppic"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Figure 1: US Crime Gaussian regression (prob5.csv, precision in [0, 100])
# ------------------------------------------------------------------------------
cat("Generating Figure 1 (n100ci1.pdf and n100ci2.pdf)...\n")
p5 <- read.csv("prob5.csv")

# Left: Model M1
pdf(file.path(output_dir, "n100ci1.pdf"), width = 6.2, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(x = p5$V1, y = p5$V3, xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", pch = 1)
mtext("INLA", side = 2, line = 2.8, cex = 1.7)
mtext("Chib's method", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
dev.off()

# Right: Model M2
pdf(file.path(output_dir, "n100ci2.pdf"), width = 6.2, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(x = p5$V2, y = p5$V4, xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", pch = 1)
mtext("INLA", side = 2, line = 2.8, cex = 1.7)
mtext("Chib's method", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
dev.off()

# ------------------------------------------------------------------------------
# Figure 2: Convergence for Gaussian model with MCMC iterations
# ------------------------------------------------------------------------------
cat("Generating Figure 2 (citune1.pdf)...\n")
fig2_file <- "fig2_data.rds"

if (file.exists(fig2_file)) {
  fig2_data <- readRDS(fig2_file)
  res_fig2 <- fig2_data$res
  m1_mlik <- fig2_data$m1_mlik
  m2_mlik <- fig2_data$m2_mlik
} else {
  simx <- read.table("simcen-x1.txt", sep = ",")
  simy <- read.table("simcen-y1.txt")
  data.example <- cbind(simy, simx)
  names(data.example)[1] = "Y"
  formula1 <- as.formula(paste(colnames(data.example)[1], "~ 1 +", paste0(colnames(data.example)[-c(1,3,5,7,9,11,13,15)], collapse = "+")))
  
  meanp <- 1
  precp <- 0.2
  res_fig2 <- array(0, dim = c(10, 10, 2))
  
  set.seed(42)
  for (i in 1:10) {
    for (j in 1:10) {
      seed <- i + i * j + j * runif(n = 1, min = 0, max = 1000)
      m1C <- MCMCregress(formula1, seed = seed, burnin = 100, b0 = meanp, data = data.example,
                         B0 = precp, c0 = 2, d0 = 2, marginal.likelihood = "Chib95", mcmc = 100 * 2^i, verbose = 0)
      BF <- BayesFactor(m1C)
      res_fig2[i, j, 1] <- BF$BF.logmarglike
      res_fig2[i, j, 2] <- 100 * 2^i
    }
  }
  
  M1 <- inla(formula = formula1, data = data.example, family = "Gaussian", control.predictor = list(compute = TRUE),
             control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
             control.fixed = list(mean = c(default = meanp), mean.intercept = meanp, prec.intercept = precp, prec = c(default = precp)),
             control.inla = list(strategy = "laplace", npoints = 100, diff.logdens = 2.50, int.strategy = "grid", dz = 1.913, interpolator = "gaussian")
  )
  M2 <- inla(formula = formula1, data = data.example, family = "Gaussian", control.predictor = list(compute = TRUE),
             control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
             control.fixed = list(mean = c(default = meanp), mean.intercept = meanp, prec.intercept = precp, prec = c(default = precp))
  )
  m1_mlik <- M1$mlik[1]
  m2_mlik <- M2$mlik[1]
  saveRDS(list(res = res_fig2, m1_mlik = m1_mlik, m2_mlik = m2_mlik), fig2_file)
}

pdf(file.path(output_dir, "citune1.pdf"), width = 12.5, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(y = res_fig2[, , 1], x = res_fig2[, , 2], type = "p", col = 2,
     ylim = c(min(min(res_fig2[, , 1]), m2_mlik), max(max(res_fig2[, , 1]), m2_mlik)),
     xlab = "", ylab = "", yaxt = "n", xaxt = "n")
mtext("MLIK", side = 2, line = 2.8, cex = 1.7)
mtext("Number of MCMC iterations", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
abline(a = m2_mlik, b = 0, col = 4) # default INLA (dark blue)
dev.off()

# ------------------------------------------------------------------------------
# Figure 3: Probit regression (prob2.csv, precision in [0, 10])
# ------------------------------------------------------------------------------
cat("Generating Figure 3 (l10ci1.pdf and l10ci2.pdf)...\n")
p2 <- read.csv("prob2.csv")

# Left: Model M1
pdf(file.path(output_dir, "l10ci1.pdf"), width = 6.2, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(x = p2$V1, y = p2$V3, xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", pch = 1)
mtext("INLA", side = 2, line = 2.8, cex = 1.7)
mtext("Chib's method", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
dev.off()

# Right: Model M2
pdf(file.path(output_dir, "l10ci2.pdf"), width = 6.2, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(x = p2$V2, y = p2$V4, xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", pch = 1)
mtext("INLA", side = 2, line = 2.8, cex = 1.7)
mtext("Chib's method", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
dev.off()

# ------------------------------------------------------------------------------
# Figure 4: Probit regression convergence with MCMC iterations
# ------------------------------------------------------------------------------
cat("Generating Figure 4 (lcitune1.pdf)...\n")
fig4_file <- "fig4_data.rds"

if (file.exists(fig4_file)) {
  fig4_data <- readRDS(fig4_file)
  res_fig4 <- fig4_data$res
  m1_probit_mlik <- fig4_data$m1_mlik
} else {
  simx <- read.table("sim3-X.txt", sep = ",")
  simy <- read.table("sim3-Y.txt", sep = ",")
  data.example <- cbind(simy, simx)
  names(data.example)[1] = "Y1"
  data.example$V2 <- (data.example$V10 + data.example$V14) * data.example$V9
  data.example$V5 <- (data.example$V11 + data.example$V15) * data.example$V12
  formula1 <- as.formula(paste(colnames(data.example)[1], "~ 1 +", paste0(colnames(data.example)[-c(1,3,5,7,9,11,13,15,17,19)], collapse = "+")))
  
  meanp <- 1
  precp <- 0.2
  n <- 2000
  
  grid_df <- expand.grid(i = 1:10, j = 1:10)
  cores <- min(detectCores(), 10)
  
  res_list <- mclapply(1:nrow(grid_df), function(idx) {
    i <- grid_df$i[idx]
    j <- grid_df$j[idx]
    seed <- i + i * j + j * (idx * 3.1415)
    m1C <- MCMCprobit(formula1, seed = seed, burnin = 1000, b0 = meanp, data = data.example,
                      B0 = precp, marginal.likelihood = "Chib95", mcmc = 100 * 2^i, verbose = 0)
    BF <- BayesFactor(m1C)
    list(i = i, j = j, val = BF$BF.logmarglike, iter = 100 * 2^i)
  }, mc.cores = cores)
  
  res_fig4 <- array(0, dim = c(10, 10, 2))
  for (item in res_list) {
    res_fig4[item$i, item$j, 1] <- item$val
    res_fig4[item$i, item$j, 2] <- item$iter
  }
  
  M1 <- inla(formula = formula1, data = data.example, family = "binomial", Ntrials = rep(1, n),
             control.family = list(link = "probit"), control.predictor = list(compute = TRUE),
             control.fixed = list(mean = c(default = meanp), mean.intercept = meanp, prec.intercept = precp, prec = c(default = precp)))
  m1_probit_mlik <- M1$mlik[1]
  saveRDS(list(res = res_fig4, m1_mlik = m1_probit_mlik), fig4_file)
}

pdf(file.path(output_dir, "lcitune1.pdf"), width = 12.5, height = 5.0)
par(mar = c(5, 5, 2, 2))
plot(y = res_fig4[, , 1], x = res_fig4[, , 2], type = "p", col = 2,
     ylim = c(min(min(res_fig4[, , 1]), m1_probit_mlik), max(max(res_fig4[, , 1]), m1_probit_mlik)),
     xlab = "", ylab = "", yaxt = "n", xaxt = "n")
mtext("MLIK", side = 2, line = 2.8, cex = 1.7)
mtext("Number of MCMC iterations", side = 1, line = 3.5, cex = 1.7)
axis(2, cex.axis = 1.5)
axis(1, cex.axis = 1.5)
abline(a = m1_probit_mlik, b = 0, col = 4) # INLA default (dark blue)
dev.off()

cat("All 4 figures successfully generated in PDF format!\n")
