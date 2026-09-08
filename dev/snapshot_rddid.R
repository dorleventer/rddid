# Numerical regression snapshot of rddid outputs. Run from the package root:
#   Rscript dev/snapshot_rddid.R out.rds > out.txt
# before and after any change to rd_period / aggregate.R / bandwidth.R, then diff the .txt files.
# Under the App.-B-synced package, p = 1 paths are identical to the pre-sync snapshot except the
# PC cross term of bwselect = "iter" under scheme "pc" (see dev/appB_map.md §3).
suppressMessages(pkgload::load_all(".", quiet = TRUE))
out_path <- commandArgs(trailingOnly = TRUE)[1]

make_panel <- function(scheme, n = 1500, seed = 7) {
  set.seed(seed)
  Tn <- 3; t_rd <- 3
  R0 <- runif(n, -1, 1)
  rows <- list()
  for (t in 1:Tn) {
    R <- switch(scheme,
      cs = runif(n, -1, 1),
      pc = R0,
      pv = pmin(1, pmax(-1, R0 + rnorm(n, 0, 0.15))))
    id <- if (scheme == "cs") (t - 1) * n + seq_len(n) else seq_len(n)
    # confounding jump 0.4 in every period, treatment effect 1 in t_rd,
    # period-specific curvature so the biases do not cancel
    m <- 0.3 * t + 0.5 * R + (0.4 + 0.3 * t) * R^2 + 0.4 * (R >= 0) + (t == t_rd) * 1 * (R >= 0)
    u <- if (scheme == "cs") rnorm(n, 0, 0.5) else 0.35 * rnorm(n) + 0.35 * rep(rnorm(n, 0, 1), 1)  # unit component reused below
    rows[[t]] <- data.frame(id = id, year = t, R = R, m = m, u = u)
  }
  d <- do.call(rbind, rows)
  if (scheme != "cs") {              # shared unit effect across periods -> cross-period error covariance
    set.seed(seed + 1)
    ue <- rnorm(n, 0, 0.4); d$u <- 0.35 * rnorm(nrow(d)) + ue[d$id]
  }
  d$Y <- d$m + d$u
  d
}

res <- list()
for (sc in c("cs", "pc", "pv")) {
  d <- make_panel(sc)
  for (bw in c("joint", "iter", "cct")) for (p in c(1L, 2L)) {
    key <- sprintf("%s|%s|p%d", sc, bw, p)
    r <- tryCatch(rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                        bwselect = bw, scheme = "auto", p = p, q = p + 1L),
                  error = function(e) e)
    if (inherits(r, "error")) { res[[key]] <- list(error = conditionMessage(r)); next }
    bwi <- r$bandwidth
    res[[key]] <- list(
      scheme = r$scheme,
      estimates = r$estimates,
      bws = lapply(r$fits, function(f) c(h = f$h, b = f$b)),
      bw_info = bwi[setdiff(names(bwi), c("bws", "pilot"))],
      consts = t(sapply(r$fits, function(f) c(b_const = f$b_const, v_const = f$v_const, n = f$n))))
  }
  # iter seeded from cct too (p = 1)
  key <- sprintf("%s|iter-cct|p1", sc)
  r <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
             bwselect = "iter", start = "cct", scheme = "auto")
  res[[key]] <- list(scheme = r$scheme, estimates = r$estimates,
                     bws = lapply(r$fits, function(f) c(h = f$h, b = f$b)),
                     niter = r$bandwidth$niter)
  # forced-scheme variants of iter (p = 1): pc and cs on the same data
  for (fs in c("cs", "pc")) {
    key <- sprintf("%s|iter-force%s|p1", sc, fs)
    r <- rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
               bwselect = "iter", scheme = fs)
    res[[key]] <- list(estimates = r$estimates,
                       bws = lapply(r$fits, function(f) c(h = f$h, b = f$b)),
                       niter = r$bandwidth$niter)
  }
}
# rd_period constants at a fixed bandwidth, p = 1 and p = 2
d <- make_panel("pc"); d3 <- d[d$year == 3, ]
for (p in 1:2) {
  f <- rd_period(d3$Y, d3$R, h = 0.35, b = 0.5, id = d3$id, p = p, q = p + 1L)
  res[[sprintf("rd_period|p%d", p)]] <- f[c("D", "D_bc", "V_D", "V_D_bc", "b_const", "v_const", "n")]
}
saveRDS(res, out_path)
# human-readable dump
for (k in names(res)) {
  cat("==", k, "\n")
  x <- res[[k]]
  if (!is.null(x$error)) { cat("  ERROR:", x$error, "\n"); next }
  if (!is.null(x$estimates)) print(round(as.matrix(x$estimates), 6))
  if (!is.null(x$bws)) print(round(do.call(rbind, x$bws), 6))
  if (!is.null(x$bw_info)) str(x$bw_info)
  if (!is.null(x$consts)) print(signif(x$consts, 6))
  if (!is.null(x$niter)) cat("  niter", x$niter, "\n")
  if (is.null(x$estimates)) print(unlist(x))
}
