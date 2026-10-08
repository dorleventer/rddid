# Golden-master snapshot of every exported function (UX sweep gate G1, T1644).
#
#   Rscript dev/snapshot_all.R dev/snapshots/<name>.rds        # record
#   Rscript dev/snapshot_compare.R <baseline.rds> <new.rds>    # compare (exit 1 on any change)
#
# Every exported function is called on seeded data under every option it accepts, and the FULL
# return object is stored. The compare script requires every leaf of the baseline to be
# identical() in the new object; new leaves are allowed (listed), changed or missing ones fail.
# Field `call` is dropped (it records the call text, not a number).
#
# The calls below use only argument names that survive the sweep (renamed arguments keep the old
# name as an alias), so this script must run unchanged on every commit of the sweep.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
out_path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_path)) stop("usage: Rscript dev/snapshot_all.R <out.rds>")

# ---- data ----------------------------------------------------------------------------------
# One panel per sampling scheme. pv has real switchers (movement sd 0.15 against a bandwidth
# of ~0.3), so every composition test has populated cells.
make_panel <- function(scheme, n = 1500, n_periods = 3, seed = 7) {
  set.seed(seed)
  t_rd <- n_periods
  R0 <- runif(n, -1, 1)
  ue <- rnorm(n, 0, 0.4)                       # unit effect (panel schemes)
  rows <- lapply(seq_len(n_periods), function(t) {
    R <- switch(scheme,
      cs = runif(n, -1, 1),
      pc = R0,
      pv = pmin(1, pmax(-1, R0 + rnorm(n, 0, 0.15))))
    id <- if (scheme == "cs") (t - 1) * n + seq_len(n) else seq_len(n)
    # confounding jump 0.4 + 0.1 t in every period, treatment effect 1 in t_rd,
    # period-specific curvature so the biases do not cancel
    m <- 0.3 * t + 0.5 * R + (0.4 + 0.3 * t) * R^2 + (0.4 + 0.1 * t) * (R >= 0) +
      (t == t_rd) * 1 * (R >= 0)
    u <- if (scheme == "cs") rnorm(n, 0, 0.5) else 0.35 * rnorm(n) + ue
    data.frame(id = id, year = t, R = R, Y = m + u)
  })
  do.call(rbind, rows)
}
dat <- list(cs = make_panel("cs"), pc = make_panel("pc"), pv = make_panel("pv"),
            pv4 = make_panel("pv", n_periods = 4, seed = 11))

res <- list()
keep <- function(key, expr) {
  r <- tryCatch(expr, error = function(e) structure(list(error = conditionMessage(e)),
                                                      class = "snapshot_error"))
  if (!inherits(r, "snapshot_error")) r$call <- NULL
  res[[key]] <<- r
  cat(sprintf("%-60s %s\n", key, if (inherits(r, "snapshot_error")) paste("ERROR:", r$error) else "ok"))
}

# ---- rddid() -------------------------------------------------------------------------------
for (sc in c("cs", "pc", "pv")) {
  d <- dat[[sc]]
  for (bw in c("joint", "iter", "cct")) {
    keep(sprintf("rddid|%s|%s|constant", sc, bw),
         rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = bw))
    keep(sprintf("rddid|%s|%s|linear", sc, bw),
         rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = bw,
               weights = "linear"))
  }
  keep(sprintf("rddid|%s|joint|p2", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "joint",
             p = 2L, q = 3L))
  keep(sprintf("rddid|%s|cct|p2", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
             p = 2L, q = 3L))
  keep(sprintf("rddid|%s|iter|start-cct", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "iter",
             start = "cct"))
  keep(sprintf("rddid|%s|joint|noreg", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "joint",
             regularize = FALSE))
  keep(sprintf("rddid|%s|fixed-h", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = 0.3, b = 0.45))
  keep(sprintf("rddid|%s|joint|atu", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "joint",
             estimand = "atu"))
  keep(sprintf("rddid|%s|joint|numeric-w", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "joint",
             comparisons = c(2, 1), weights = c(0.3, 0.7)))
  for (fs in c("cs", "pc", "pv"))
    keep(sprintf("rddid|%s|joint|scheme-%s", sc, fs),
         rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "joint",
               scheme = fs))
  for (k in c("epanechnikov", "uniform"))
    keep(sprintf("rddid|%s|cct|%s", sc, k),
         rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
               kernel = k))
  keep(sprintf("rddid|%s|cct|level90", sc),
       rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = "cct",
             level = 0.90))
}
keep("rddid|pc|no-id",
     rddid(dat$pc, y = "Y", x = "R", time = "year", t_rd = 3, bwselect = "cct"))
keep("rddid|pv4|joint|linear-3comp",
     rddid(dat$pv4, y = "Y", x = "R", time = "year", id = "id", t_rd = 4, bwselect = "joint",
           weights = "linear"))

# ---- rd_period(), rd_bw_cct() --------------------------------------------------------------
d3 <- dat$pv[dat$pv$year == 3, ]
for (k in c("triangular", "epanechnikov", "uniform")) {
  keep(sprintf("rd_period|%s|p1", k),
       rd_period(d3$Y, d3$R, h = 0.3, b = 0.45, id = d3$id, c = 0, kernel = k))
  keep(sprintf("rd_period|%s|p2", k),
       rd_period(d3$Y, d3$R, h = 0.4, b = 0.5, id = d3$id, c = 0, p = 2L, q = 3L, kernel = k))
  keep(sprintf("rd_bw_cct|%s|p1", k), as.list(rd_bw_cct(d3$Y, d3$R, c = 0, kernel = k)))
  keep(sprintf("rd_bw_cct|%s|p2", k), as.list(rd_bw_cct(d3$Y, d3$R, c = 0, p = 2L, kernel = k)))
}
keep("rd_period|b-eq-h", rd_period(d3$Y, d3$R, h = 0.3))
keep("rd_period|cutoff-0.1", rd_period(d3$Y, d3$R, h = 0.3, b = 0.4, c = 0.1))

# ---- validation tests ------------------------------------------------------------------------
for (sc in c("pc", "pv")) {
  d <- dat[[sc]]
  for (bw in c("cct", "rot")) for (bc in c(TRUE, FALSE)) {
    tag <- sprintf("%s|%s|bc%d", sc, bw, bc)
    keep(sprintf("rd_typecont|%s", tag),
         rd_typecont(d, x = "R", time = "year", id = "id", bwselect = bw, bc = bc))
    keep(sprintf("rd_compstable|%s", tag),
         rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3, bwselect = bw, bc = bc))
    keep(sprintf("rd_homog|%s", tag),
         rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = bw,
                  bc = bc))
    keep(sprintf("rd_trendcell|%s", tag),
         rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, bwselect = bw,
                      bc = bc))
  }
  keep(sprintf("rd_typecont|%s|atu", sc),
       rd_typecont(d, x = "R", time = "year", id = "id", estimand = "atu"))
  keep(sprintf("rd_compstable|%s|atu", sc),
       rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu"))
  keep(sprintf("rd_homog|%s|atu", sc),
       rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu"))
  keep(sprintf("rd_trendcell|%s|atu", sc),
       rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu"))
  keep(sprintf("rd_typecont|%s|scheme-cs", sc),
       rd_typecont(d, x = "R", time = "year", id = "id", scheme = "cs"))
  keep(sprintf("rd_compstable|%s|scheme-cs", sc),
       rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3, scheme = "cs"))
  keep(sprintf("rd_homog|%s|pattern", sc),
       rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, type_by = "pattern"))
  keep(sprintf("rd_trendcell|%s|pattern", sc),
       rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                    type_by = "pattern"))
  keep(sprintf("rd_homog|%s|fixed-h", sc),
       rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, h = 0.3))
  keep(sprintf("rd_compstable|%s|fixed-h", sc),
       rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3, h = 0.3))
  keep(sprintf("rd_homog|%s|min_n-50", sc),
       rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, min_n = 50L))
  keep(sprintf("rd_homog|%s|p2-dots", sc),
       rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, p = 2L, q = 3L))
  keep(sprintf("rd_compstable|%s|comparisons-2", sc),
       rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3, comparisons = 2))
}
keep("rd_trendcell|pv4|linear",
     rd_trendcell(dat$pv4, y = "Y", x = "R", time = "year", id = "id", t_rd = 4,
                  trend = "linear"))
keep("rd_trendcell|pv4|constant",
     rd_trendcell(dat$pv4, y = "Y", x = "R", time = "year", id = "id", t_rd = 4))
keep("rd_typecont|pv4", rd_typecont(dat$pv4, x = "R", time = "year", id = "id"))

attr(res, "git_head") <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                                  error = function(e) NA_character_)
attr(res, "recorded") <- format(Sys.time())
saveRDS(res, out_path)
cat(sprintf("\n%d cells (%d errors) -> %s  [HEAD %s]\n", length(res),
            sum(vapply(res, inherits, logical(1), "snapshot_error")), out_path,
            attr(res, "git_head")))
