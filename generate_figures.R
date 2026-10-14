# ==============================================================================
# generate_figures.R: Generate Publication-Quality Vector PDF Figures for JRS
# ==============================================================================
# Incorporates all reviewer recommendations:
# - Identity line (y = x) added to Figures 1 and 3
# - Explicit scale labels: "Estimated log marginal likelihood"
# - In-panel legends for Figures 2 and 4 distinguishing markers and line styles
# ==============================================================================

source("R/common.R")

cat_header("Generating Publication-Quality Vector PDF Figures (apppic/)")

output_dir <- APPPIC
cat(sprintf("Output directory for figure PDFs: %s\n", output_dir))

# ------------------------------------------------------------------------------
# Figure 1: US Crime Gaussian regression (prob5.csv, precision in [0, 100])
# ------------------------------------------------------------------------------
cat("Generating Figure 1 (n100ci1.pdf and n100ci2.pdf)...\n")
p5 <- read.csv("prob5.csv")

# Left panel: Model M1
pdf(file.path(output_dir, "n100ci1.pdf"), width = 6.2, height = 5.2)
par(mar = c(5, 5, 2, 2))
lim_m1 <- range(c(p5$V1, p5$V3))
plot(x = p5$V1, y = p5$V3, xlim = lim_m1, ylim = lim_m1,
     xlab = "", ylab = "", main = "", pch = 1, col = "black", cex = 1.1)
abline(a = 0, b = 1, col = "gray40", lty = 2, lwd = 1.5) # Identity line y = x
mtext("INLA (log marginal likelihood)", side = 2, line = 2.8, cex = 1.3)
mtext("Chib's method (log marginal likelihood)", side = 1, line = 3.2, cex = 1.3)
legend("topleft", legend = c("Estimates", "Identity line (y = x)"),
       pch = c(1, NA), lty = c(NA, 2), col = c("black", "gray40"),
       lwd = c(NA, 1.5), bty = "n", cex = 1.1)
dev.off()

# Right panel: Model M2
pdf(file.path(output_dir, "n100ci2.pdf"), width = 6.2, height = 5.2)
par(mar = c(5, 5, 2, 2))
lim_m2 <- range(c(p5$V2, p5$V4))
plot(x = p5$V2, y = p5$V4, xlim = lim_m2, ylim = lim_m2,
     xlab = "", ylab = "", main = "", pch = 1, col = "black", cex = 1.1)
abline(a = 0, b = 1, col = "gray40", lty = 2, lwd = 1.5) # Identity line y = x
mtext("INLA (log marginal likelihood)", side = 2, line = 2.8, cex = 1.3)
mtext("Chib's method (log marginal likelihood)", side = 1, line = 3.2, cex = 1.3)
legend("topleft", legend = c("Estimates", "Identity line (y = x)"),
       pch = c(1, NA), lty = c(NA, 2), col = c("black", "gray40"),
       lwd = c(NA, 1.5), bty = "n", cex = 1.1)
dev.off()

# ------------------------------------------------------------------------------
# Figure 2: Convergence for Gaussian model with MCMC iterations
# ------------------------------------------------------------------------------
cat("Generating Figure 2 (citune1.pdf)...\n")
fig2_file <- "fig2_data.rds"
if (!file.exists(fig2_file)) {
  source("03_convergence_figs.R")
}
fig2_data <- readRDS(fig2_file)
res_fig2 <- fig2_data$res
m2_mlik <- fig2_data$m2_mlik

pdf(file.path(output_dir, "citune1.pdf"), width = 11.5, height = 5.2)
par(mar = c(5, 5, 2, 2))
y_range_f2 <- range(c(res_fig2[, , 1], m2_mlik))
# Add small padding
y_pad_f2 <- 0.05 * diff(y_range_f2)
ylim_f2 <- c(y_range_f2[1] - y_pad_f2, y_range_f2[2] + y_pad_f2)

plot(y = res_fig2[, , 1], x = res_fig2[, , 2], type = "p", col = 2, pch = 1, cex = 1.2,
     ylim = ylim_f2, xlab = "", ylab = "", log = "x")
abline(h = m2_mlik, col = 4, lwd = 2) # default INLA horizontal line
mtext("Estimated log marginal likelihood", side = 2, line = 2.8, cex = 1.3)
mtext("MCMC iteration count (log scale)", side = 1, line = 3.2, cex = 1.3)
legend("topright", legend = c("Chib (1995) replicate", "INLA default estimate"),
       pch = c(1, NA), lty = c(NA, 1), col = c(2, 4),
       pt.cex = 1.2, lwd = c(NA, 2), bty = "n", cex = 1.2)
dev.off()

# ------------------------------------------------------------------------------
# Figure 3: Probit regression (prob2.csv, precision in [0, 10])
# ------------------------------------------------------------------------------
cat("Generating Figure 3 (l10ci1.pdf and l10ci2.pdf)...\n")
p2 <- read.csv("prob2.csv")

# Left panel: Model M1
pdf(file.path(output_dir, "l10ci1.pdf"), width = 6.2, height = 5.2)
par(mar = c(5, 5, 2, 2))
lim_p_m1 <- range(c(p2$V1, p2$V3))
plot(x = p2$V1, y = p2$V3, xlim = lim_p_m1, ylim = lim_p_m1,
     xlab = "", ylab = "", main = "", pch = 1, col = "black", cex = 1.1)
abline(a = 0, b = 1, col = "gray40", lty = 2, lwd = 1.5) # Identity line
mtext("INLA (log marginal likelihood)", side = 2, line = 2.8, cex = 1.3)
mtext("Chib's method (log marginal likelihood)", side = 1, line = 3.2, cex = 1.3)
legend("topleft", legend = c("Estimates", "Identity line (y = x)"),
       pch = c(1, NA), lty = c(NA, 2), col = c("black", "gray40"),
       lwd = c(NA, 1.5), bty = "n", cex = 1.1)
dev.off()

# Right panel: Model M2
pdf(file.path(output_dir, "l10ci2.pdf"), width = 6.2, height = 5.2)
par(mar = c(5, 5, 2, 2))
lim_p_m2 <- range(c(p2$V2, p2$V4))
plot(x = p2$V2, y = p2$V4, xlim = lim_p_m2, ylim = lim_p_m2,
     xlab = "", ylab = "", main = "", pch = 1, col = "black", cex = 1.1)
abline(a = 0, b = 1, col = "gray40", lty = 2, lwd = 1.5) # Identity line
mtext("INLA (log marginal likelihood)", side = 2, line = 2.8, cex = 1.3)
mtext("Chib's method (log marginal likelihood)", side = 1, line = 3.2, cex = 1.3)
legend("topleft", legend = c("Estimates", "Identity line (y = x)"),
       pch = c(1, NA), lty = c(NA, 2), col = c("black", "gray40"),
       lwd = c(NA, 1.5), bty = "n", cex = 1.1)
dev.off()

# ------------------------------------------------------------------------------
# Figure 4: Probit regression convergence with MCMC iterations
# ------------------------------------------------------------------------------
cat("Generating Figure 4 (lcitune1.pdf)...\n")
fig4_file <- "fig4_data.rds"
if (!file.exists(fig4_file)) {
  source("03_convergence_figs.R")
}
fig4_data <- readRDS(fig4_file)
res_fig4 <- fig4_data$res
m1_probit_mlik <- fig4_data$m1_mlik

pdf(file.path(output_dir, "lcitune1.pdf"), width = 11.5, height = 5.2)
par(mar = c(5, 5, 2, 2))
y_range_f4 <- range(c(res_fig4[, , 1], m1_probit_mlik))
y_pad_f4 <- 0.05 * diff(y_range_f4)
ylim_f4 <- c(y_range_f4[1] - y_pad_f4, y_range_f4[2] + y_pad_f4)

plot(y = res_fig4[, , 1], x = res_fig4[, , 2], type = "p", col = 2, pch = 1, cex = 1.2,
     ylim = ylim_f4, xlab = "", ylab = "", log = "x")
abline(h = m1_probit_mlik, col = 4, lwd = 2) # INLA line
mtext("Estimated log marginal likelihood", side = 2, line = 2.8, cex = 1.3)
mtext("MCMC iteration count (log scale)", side = 1, line = 3.2, cex = 1.3)
legend("topright", legend = c("Chib (1995) replicate", "INLA default estimate"),
       pch = c(1, NA), lty = c(NA, 1), col = c(2, 4),
       pt.cex = 1.2, lwd = c(NA, 2), bty = "n", cex = 1.2)
dev.off()

cat("All 4 figures successfully regenerated and saved to apppic/!\n")
