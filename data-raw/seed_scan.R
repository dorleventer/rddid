# Seed scan behind the choice of seed for rddid_sim_pv (data-raw/make_sim.R). Run from the package root:
#   Rscript data-raw/seed_scan.R
# 40 seeds of the moving-R design: estimate, and the p-values of the four tests. In this design
# composition stability genuinely fails (the yearly deviations of R from the base position have a
# spread that grows with the year, so the mix of units near the cutoff differs by year; rejection
# rate ~0.8 across seeds), while type
# continuity, homogeneous confounding and the constant within-type jump hold. The shipped seed is the
# one where the holding assumptions are comfortably not rejected and the failing one clearly is.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("data-raw/make_sim.R")
src <- src[!grepl("use_data|^rddid_sim|^sw <-|^cat[(]", src)]
eval(parse(text = src))
res <- do.call(rbind, lapply(1:40, function(sd) {
  d <- make_sim(moving = TRUE, seed = 20261000 + sd)
  f <- suppressMessages(rddid(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3))
  p <- function(expr) tryCatch(suppressMessages(expr)$p_value, error = function(e) NA_real_)
  data.frame(seed = 20261000 + sd, est = f$estimates$est[1],
             typecont   = p(rd_typecont(d, x = "R", time = "year", id = "id")),
             compstable = p(rd_compstable(d, x = "R", time = "year", id = "id", t_rd = 3)),
             homog      = p(rd_homog(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)),
             trendcell  = p(rd_trendcell(d, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)))
}))
cat(sprintf("rejection rates at 5%%: typecont %.2f  compstable %.2f  homog %.2f  trendcell %.2f\n",
            mean(res$typecont < .05), mean(res$compstable < .05), mean(res$homog < .05, na.rm = TRUE),
            mean(res$trendcell < .05, na.rm = TRUE)))
cat(sprintf("mean est %.3f (truth 1), sd %.3f\n", mean(res$est), sd(res$est)))
ok <- res[res$typecont > .2 & res$homog > .2 & res$trendcell > .2, ]
print(round(ok, 3))
