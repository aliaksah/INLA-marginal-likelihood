# ==============================================================================
# 00_environment.R: Record Software Environment and Inspect INLA Modes
# ==============================================================================
source("R/common.R")

cat_header("00: Software Environment and INLA Strategy / Mode Assessment")

# 1. Capture sessionInfo
os_info <- Sys.info()
log_fname <- sprintf("logs/sessionInfo_%s_%s.txt", tolower(os_info["sysname"]), tolower(os_info["machine"]))
sink(log_fname)
cat("=== R SESSION INFO ===\n")
print(sessionInfo())
cat("\n=== INLA VERSION & BUILD INFO ===\n")
cat("INLA Version:", as.character(packageVersion("INLA")), "\n")
tryCatch({
  cat("INLA OS:", INLA:::inla.os(), "\n")
  cat("INLA 64bit:", INLA:::inla.os.64bit(), "\n")
}, error = function(e) NULL)
cat("\n=== CRAN PACKAGES ===\n")
pkgs <- c("MCMCpack", "MASS", "geepack", "coda", "Matrix", "mvtnorm")
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("%-12s : %s\n", p, as.character(packageVersion(p))))
  }
}
sink()
cat(sprintf("Session info logged to %s\n", log_fname))

# 2. Check seizure data loader consistency
cat("Checking seizure dataset consistency between geepack and epileptic.txt...\n")
if (requireNamespace("geepack", quietly = TRUE)) {
  s_gee <- {
    e <- new.env(); utils::data("seizure", package = "geepack", envir = e); s <- e$seizure
    s <- s[s$base != 151, c("base", "y1", "y2", "y3", "y4", "trt")]
    s$id <- seq_len(nrow(s))
    rownames(s) <- NULL
    s
  }
  raw <- read.table("epileptic.txt", sep = ";", header = TRUE)
  w <- reshape(raw[, c("id", "visit", "seizure")], timevar = "visit", idvar = "id", direction = "wide")
  s_txt <- data.frame(base = w$seizure.0, y1 = w$seizure.1, y2 = w$seizure.2,
                      y3 = w$seizure.3, y4 = w$seizure.4,
                      trt = raw$trt[match(w$id, raw$id)])
  s_txt <- s_txt[s_txt$base != 151, ]
  s_txt$id <- seq_len(nrow(s_txt))
  rownames(s_txt) <- NULL
  
  ord_gee <- order(s_gee$trt, s_gee$base, s_gee$y1, s_gee$y2, s_gee$y3, s_gee$y4)
  ord_txt <- order(s_txt$trt, s_txt$base, s_txt$y1, s_txt$y2, s_txt$y3, s_txt$y4)
  match_ok <- isTRUE(all.equal(s_gee[ord_gee, c("base","y1","y2","y3","y4","trt")],
                               s_txt[ord_txt, c("base","y1","y2","y3","y4","trt")],
                               check.attributes = FALSE))
  cat(sprintf("  geepack vs fallback text identical: %s\n", match_ok))
  if (!match_ok) stop("Discrepancy detected in seizure data loader!")
}

# 3. Test INLA classic vs compact mode across models
cat("\nComparing INLA inla.mode = 'classic' vs 'compact' across models...\n")
comp_res <- list()

# (a) Linear Gaussian (Table 1, sigma0 = 1000)
dat1 <- data.frame(y = 2, x = 1)
inla.setOption("inla.mode", "classic")
m1_cls <- inla(y ~ x - 1, data = dat1, control.compute = list(mlik = TRUE),
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.fixed = list(mean = 0, prec = (1/1000)^2))
inla.setOption("inla.mode", "compact")
m1_cmp <- inla(y ~ x - 1, data = dat1, control.compute = list(mlik = TRUE),
               control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
               control.fixed = list(mean = 0, prec = (1/1000)^2))
comp_res[[length(comp_res) + 1]] <- data.frame(
  Example = "Ex 1: Linear Normal (s0=1000)",
  Classic_Int = round(m1_cls$mlik[1,1], 4), Classic_Gau = round(m1_cls$mlik[2,1], 4),
  Compact_Int = round(m1_cmp$mlik[1,1], 4), Compact_Gau = round(m1_cmp$mlik[2,1], 4)
)

# (b) US Crime (Table 2, M1, sbeta = 1000)
crime <- load_crime()
inla.setOption("inla.mode", "classic")
m2_cls <- inla_gauss(crime$f1, crime$data, (1/1000)^2)
inla.setOption("inla.mode", "compact")
m2_cmp <- inla_gauss(crime$f1, crime$data, (1/1000)^2)
comp_res[[length(comp_res) + 1]] <- data.frame(
  Example = "Ex 2: Crime M1 (sbeta=1000)",
  Classic_Int = round(m2_cls$mlik[1,1], 4), Classic_Gau = round(m2_cls$mlik[2,1], 4),
  Compact_Int = round(m2_cmp$mlik[1,1], 4), Compact_Gau = round(m2_cmp$mlik[2,1], 4)
)

# (c) Probit (Table 3, M1, sbeta = 10)
probit <- load_probit()
inla.setOption("inla.mode", "classic")
m3_cls <- inla_binom(probit$f1, probit$data, "probit", (1/10)^2)
# Compact mode known to fail with index assertion on some platforms
inla.setOption("inla.mode", "compact")
m3_cmp_int <- NA; m3_cmp_gau <- NA
tryCatch({
  m3_cmp <- inla_binom(probit$f1, probit$data, "probit", (1/10)^2)
  m3_cmp_int <- round(m3_cmp$mlik[1,1], 4)
  m3_cmp_gau <- round(m3_cmp$mlik[2,1], 4)
}, error = function(e) {
  cat("  [Note: Probit compact mode threw error as expected on arm64]\n")
})
comp_res[[length(comp_res) + 1]] <- data.frame(
  Example = "Ex 3: Probit M1 (sbeta=10)",
  Classic_Int = round(m3_cls$mlik[1,1], 4), Classic_Gau = round(m3_cls$mlik[2,1], 4),
  Compact_Int = m3_cmp_int, Compact_Gau = m3_cmp_gau
)

# (d) Pima Logit (Table 4, M1, sbeta = 10)
pima <- load_pima()
inla.setOption("inla.mode", "classic")
m4_cls <- inla_binom(pima$f1, pima$data, "logit", 0.01)
inla.setOption("inla.mode", "compact")
m4_cmp <- inla_binom(pima$f1, pima$data, "logit", 0.01)
comp_res[[length(comp_res) + 1]] <- data.frame(
  Example = "Ex 4: Pima Logit M1 (sbeta=10)",
  Classic_Int = round(m4_cls$mlik[1,1], 4), Classic_Gau = round(m4_cls$mlik[2,1], 4),
  Compact_Int = round(m4_cmp$mlik[1,1], 4), Compact_Gau = round(m4_cmp$mlik[2,1], 4)
)

# (e) Seizure GLMM (Section 6)
d_glmm <- glmm_inla_data()
inla.setOption("inla.mode", "classic")
m5_cls <- inla(GLMM_FORMULA, data = d_glmm, family = "poisson",
               control.compute = list(mlik = TRUE),
               control.fixed = ctrl_fixed(0.01),
               control.family = list(link = "log"))
inla.setOption("inla.mode", "compact")
m5_cmp <- inla(GLMM_FORMULA, data = d_glmm, family = "poisson",
               control.compute = list(mlik = TRUE),
               control.fixed = ctrl_fixed(0.01),
               control.family = list(link = "log"))
comp_res[[length(comp_res) + 1]] <- data.frame(
  Example = "Ex 5: Seizure Poisson GLMM",
  Classic_Int = round(m5_cls$mlik[1,1], 4), Classic_Gau = round(m5_cls$mlik[2,1], 4),
  Compact_Int = round(m5_cmp$mlik[1,1], 4), Compact_Gau = round(m5_cmp$mlik[2,1], 4)
)

df_comp <- do.call(rbind, comp_res)
write.csv(df_comp, "results/inla_mode_comparison.csv", row.names = FALSE)
print(df_comp)

# Reset INLA to classic for all subsequent scripts
init_inla()
cat("\n00_environment.R completed successfully.\n")
