# ---------------------------------------------------------------------------
# site_dgp_check.R -- provenance for the site vignettes' DGPs (dev/site_plan.md).
# Defines DGP-A (time-invariant R, 3 periods) and DGP-B (time-varying R) exactly
# as the vignettes hand-code them, and checks every scenario claim: closed-form
# vs 4e6-draw oracle for the composition term C, per-scenario test rejection
# rates over 20 seeds, and that biases are visible at the vignette n. Run with
# Rscript from this directory; ~3 min. Not part of the package build.
# Result (2026-09-09): all claims hold; S3 uses kappa = 2, sort_p = 0.5.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(rddid))
options(warn = 1)
tic <- function() assign(".t0", proc.time()[["elapsed"]], envir = globalenv())
toc <- function(lab) cat(sprintf("   [%s: %.1fs]\n", lab, proc.time()[["elapsed"]] - .t0))

# ---------------- DGP-A: time-invariant running variable, 3 periods -----------
m_a <- function(r, theta) r + (theta / 2) * r^2 * (r >= 0)
dgp_a <- function(n = 2000, alpha = c(1, 1, 1), theta = c(2, 2, 2), tau = 1,
                  scheme = c("pc", "cs", "pv"), seed = 1) {
  scheme <- match.arg(scheme); set.seed(seed)
  R0 <- runif(n, -1, 1); u0 <- rnorm(n, 0, 0.5)
  do.call(rbind, lapply(1:3, function(t) {
    if (scheme == "cs") { R <- runif(n, -1, 1); u <- rnorm(n, 0, 0.5); id <- (t - 1) * n + seq_len(n) }
    else { R <- if (scheme == "pv") R0 + rnorm(n, 0, 0.1) else R0; u <- u0; id <- seq_len(n) }
    V <- as.integer(R >= 0)
    data.frame(id = id, t = t, R = R,
               Y = m_a(R, theta[t]) + alpha[t] * V + tau * V * (t == 3) + u + rnorm(n, 0, 0.5))
  }))
}

# ---------------- DGP-B: time-varying running variable ------------------------
m_b <- function(r) r + r^2 * (r >= 0)
dgp_b <- function(n = 4000, d = 0, alpha = c(1, 1), kappa = 0, sort_p = 0, sort_w = 0.5,
                  tau = 1, periods = 2, gamma = 0, seed = 1) {
  set.seed(seed)
  eta <- rnorm(n)
  R <- sapply(seq_len(periods), function(t) eta + d * (t - 1) + rnorm(n, 0, 0.5))
  if (sort_p > 0) {                       # sorting in the RD period on the period-1 side
    V1 <- R[, 1] >= 0; r2 <- R[, periods]
    mover <- r2 > -sort_w & r2 < 0 & V1 & runif(n) < sort_p
    R[mover, periods] <- -R[mover, periods]
  }
  V <- R >= 0
  do.call(rbind, lapply(seq_len(periods), function(t) {
    other <- if (periods == 2) V[, 3 - t] else V[, periods]   # 2p: other side; 3p: RD-period side (unused when alpha homog.)
    a <- alpha[other + 1] + gamma * t
    data.frame(id = seq_len(n), t = t, R = R[, t],
               Y = m_b(R[, t]) + a * V[, t] + tau * V[, t] * (t == periods) + kappa * eta + rnorm(n, 0, 0.5))
  }))
}

fit <- function(dat, weights = "constant", bwselect = "cct", scheme = "auto", ...) {
  r <- rddid(dat, y = "Y", x = "R", time = "t", id = "id", t_rd = max(dat$t),
             weights = weights, bwselect = bwselect, scheme = scheme, ...)
  e <- r$estimates
  c(est = e["Conventional", "est"], se = e["Conventional", "se"],
    est_bc = e["Robust", "est"], se_bc = e["Robust", "se"], scheme = r$scheme)
}

cat("=========== DGP-A ===========\n")
tic(); a1 <- dgp_a(); f <- fit(a1); toc("rddid cct pc")
cat("A const, constant w:", round(as.numeric(f[1:4]), 3), f[5], " (truth 1)\n")
f <- fit(a1, weights = "linear"); cat("A const, linear w  :", round(as.numeric(f[1:4]), 3), "\n")
a2 <- dgp_a(alpha = 0.5 + 0.5 * (1:3))
f <- fit(a2); cat("A lin, constant w  :", round(as.numeric(f[1:4]), 3), " (bias truth 0.75)\n")
f <- fit(a2, weights = "linear"); cat("A lin, linear w    :", round(as.numeric(f[1:4]), 3), " (truth 1)\n")
# per-period table
for (t in 1:3) { d <- a1[a1$t == t, ]; bw <- rd_bw_cct(d$Y, d$R); p <- rd_period(d$Y, d$R, h = bw["h"], b = bw["b"], id = d$id)
  cat(sprintf("  t=%d D=%.3f (%.3f) Dbc=%.3f (%.3f) h=%.3f n=%d\n", t, p$D, sqrt(p$V_D), p$D_bc, sqrt(p$V_D_bc), p$h, p$n)) }
# 20 seeds
res <- t(sapply(1:20, function(s) c(as.numeric(fit(dgp_a(seed = s))[1:2]), as.numeric(fit(dgp_a(alpha = 0.5 + 0.5*(1:3), seed = s), weights = "linear")[1:2]))))
cat("A 20 seeds: const mean est", round(mean(res[,1]),3), "sd", round(sd(res[,1]),3), "mean se", round(mean(res[,2]),3),
    "| linear mean est", round(mean(res[,3]),3), "sd", round(sd(res[,3]),3), "mean se", round(mean(res[,4]),3), "\n")

cat("\n--- V2 bandwidths, theta=(1,6,3) ---\n")
a3 <- dgp_a(theta = c(1, 6, 3))
for (rule in c("cct", "joint", "iter")) { tic()
  r <- rddid(a3, y = "Y", x = "R", time = "t", id = "id", t_rd = 3, bwselect = rule)
  hs <- vapply(r$fits, function(f) unname(f$h), 1); hs <- hs[order(as.numeric(names(hs)))]  # r$fits lists the RD period first
  e <- r$estimates
  cat(sprintf("%-5s h=%s  est=%.3f (%.3f) bc=%.3f (%.3f)", rule, paste(round(hs, 3), collapse = "/"), e[1,"est"], e[1,"se"], e[2,"est"], e[2,"se"])); toc(rule) }
cat("\n--- V2 schemes ---\n")
for (sc in c("cs", "pc", "pv")) { r <- rddid(dgp_a(scheme = sc), y = "Y", x = "R", time = "t", id = if (sc == "cs") NULL else "id", t_rd = 3, bwselect = "cct")
  e <- r$estimates; cat(sprintf("data %s -> detected %s | est %.3f  se_cs %.3f se_pc %.3f se_pv %.3f\n", sc, r$scheme, e[1,"est"], e[1,"se_cs"], e[1,"se_pc"], e[1,"se_pv"])) }

cat("\n=========== DGP-B oracle (4e6 draws) ===========\n")
oracle <- function(d, alpha, kappa = 0, sort_p = 0, N = 4e6, eps = 0.01) {
  set.seed(99); dat <- dgp_b(n = N, d = d, alpha = alpha, kappa = kappa, sort_p = sort_p, tau = 1)
  R1 <- dat$R[dat$t == 1]; R2 <- dat$R[dat$t == 2]; Y1 <- dat$Y[dat$t == 1]; Y2 <- dat$Y[dat$t == 2]
  V1 <- R1 >= 0; V2 <- R2 >= 0
  lim <- function(z, r, side) mean(z[if (side == "+") r >= 0 & r < eps else r < 0 & r > -eps])
  pi2p1 <- lim(V1, R2, "+"); pi2m1 <- lim(V1, R2, "-"); pi1p1 <- lim(V2, R1, "+"); pi1m1 <- lim(V2, R1, "-")
  D1 <- lim(Y1, R1, "+") - lim(Y1, R1, "-"); D2 <- lim(Y2, R2, "+") - lim(Y2, R2, "-")
  C <- (alpha[2] - alpha[1]) * (pi2p1 - pi1p1)
  cat(sprintf("  pi_2,(+)(1)=%.3f pi_2,(-)(1)=%.3f | pi_1,(+)(1)=%.3f pi_1,(-)(1)=%.3f | D2-D1=%.3f (tau=1) => bias %.3f | C=%.3f\n",
              pi2p1, pi2m1, pi1p1, pi1m1, D2 - D1, D2 - D1 - 1, C))
  if (d > 0) { s <- sqrt(0.25/1.25 + 0.25)
    cat(sprintf("  closed form: pi_1,(+)(1)=%.3f  pi_2,(+)(1)=%.3f\n", pnorm(d / s), 1 - pnorm((d / 1.25) / s))) }
}
cat("S0:"); oracle(0, c(1, 1)); cat("S1:"); oracle(0.5, c(1, 1)); cat("S2:"); oracle(0.5, c(0.5, 1.5)); cat("S3:"); oracle(0, c(1, 1), kappa = 1, sort_p = 0.5)

cat("\n=========== DGP-B tests at seed 1 + 20-seed rejection rates ===========\n")
tests <- function(dat, t_rd = 2) {
  comps <- setdiff(unique(dat$t), t_rd)
  tc <- rd_typecont(dat, x = "R", time = "t", id = "id", bwselect = "cct", S = 199L, bc = FALSE)
  cs <- rd_compstable(dat, x = "R", time = "t", id = "id", t_rd = t_rd, comparisons = comps, bwselect = "cct", S = 199L, bc = FALSE)
  hg <- rd_homog(dat, y = "Y", x = "R", time = "t", id = "id", comparisons = comps, t_rd = t_rd, bwselect = "cct", type_by = "rd_side", bc = FALSE)
  tr <- if (length(comps) >= 2) rd_trendcell(dat, y = "Y", x = "R", time = "t", id = "id", comparisons = comps, t_rd = t_rd, bwselect = "cct", type_by = "rd_side", bc = FALSE)$p_value else NA
  f <- fit(dat)
  c(typecont = tc$ll_wald$p, compstable = cs$pairs[[1]]$ll_wald$p, homog = hg$p_value, trendcell = tr, bias = as.numeric(f["est"]) - 1, se = as.numeric(f["se"]))
}
scen <- list(S0 = list(), S1 = list(d = 0.5), S2 = list(d = 0.5, alpha = c(0.5, 1.5)), S3 = list(kappa = 1, sort_p = 0.5),
             S0_3p = list(periods = 3), S4 = list(periods = 3, gamma = 0.3))
for (nm in names(scen)) { tic(); r1 <- do.call(tests, list(dat = do.call(dgp_b, c(scen[[nm]], seed = 1)), t_rd = if (!is.null(scen[[nm]]$periods)) 3 else 2))
  cat(nm, "seed1:", paste(names(r1), round(r1, 3), sep = "=", collapse = "  ")); toc("tests")
  rr <- t(sapply(2:21, function(s) do.call(tests, list(dat = do.call(dgp_b, c(scen[[nm]], seed = s)), t_rd = if (!is.null(scen[[nm]]$periods)) 3 else 2))))
  cat(sprintf("   20 seeds: rej@5%%  typecont %.2f compstable %.2f homog %.2f trendcell %.2f | bias mean %.3f sd %.3f mean se %.3f\n",
              mean(rr[,1] < .05), mean(rr[,2] < .05), mean(rr[,3] < .05), mean(rr[,4] < .05, na.rm = TRUE), mean(rr[,5]), sd(rr[,5]), mean(rr[,6]))) }
if (TRUE) { d4 <- dgp_b(periods = 3, gamma = 0.3); f <- fit(d4, weights = "linear"); cat("S4 linear weights est:", round(as.numeric(f[1:2]), 3), "(truth 1)\n") }
