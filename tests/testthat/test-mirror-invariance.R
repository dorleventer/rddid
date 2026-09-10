# Mirror invariance of the RD-DID estimator and three of the four validation
# tests (dev/atu_estimand_plan.md S1, S4 items 1-3, paper Section 6): mirroring
# the running variable (x -> -x, cutoff 0) negates the rddid() estimate and
# leaves everything else (SEs, bandwidths, scheme detection, the three tests
# other than rd_compstable()) numerically unchanged. Mirroring is applied by
# hand via mirror_x() (helper-mirror.R) -- no `estimand` argument exists in
# this worktree.
#
# On any FAILURE at the stated tolerance, do not loosen the tolerance or skip
# -- report the configuration and the max abs discrepancy.

dat_a <- dgp_a()   # 3 periods, PC scheme (defaults: n = 2000)

# =============================================================================
# A. rddid()
# =============================================================================

rddid_configs <- list(
  list(name = "dgp_a", data = dat_a, t_rd = 3, comparisons = c(1, 2), scheme = "pc"),
  list(name = "S1",    data = S1,    t_rd = 2, comparisons = 1,       scheme = "pv")
)

check_rddid_invariance <- function(cfg, bwselect, weights = "constant") {
  info0 <- sprintf("data=%s bwselect=%s weights=%s", cfg$name, bwselect, weights)

  orig <- rddid(cfg$data, y = "Y", x = "R", time = "t", id = "id",
                t_rd = cfg$t_rd, comparisons = cfg$comparisons,
                scheme = cfg$scheme, bwselect = bwselect, weights = weights)
  mir  <- rddid(mirror_x(cfg$data), y = "Y", x = "R", time = "t", id = "id",
                t_rd = cfg$t_rd, comparisons = cfg$comparisons,
                scheme = cfg$scheme, bwselect = bwselect, weights = weights)

  for (row in c("Conventional", "Robust")) {
    expect_equal(orig$estimates[row, "est"], -mir$estimates[row, "est"],
                 tolerance = 1e-10, info = paste(info0, row, "est negated"))
    for (sefield in c("se", "se_cs", "se_pc", "se_pv")) {
      expect_equal(orig$estimates[row, sefield], mir$estimates[row, sefield],
                   tolerance = 1e-10, info = paste(info0, row, sefield))
    }
    expect_equal(orig$estimates[row, "ci_l"], -mir$estimates[row, "ci_u"],
                 tolerance = 1e-10, info = paste(info0, row, "ci_l == -ci_u(mirror)"))
    expect_equal(orig$estimates[row, "ci_u"], -mir$estimates[row, "ci_l"],
                 tolerance = 1e-10, info = paste(info0, row, "ci_u == -ci_l(mirror)"))
  }

  # per-period bandwidths, indexed BY NAME (fits lists the RD period first)
  pks <- names(orig$fits)
  expect_identical(pks, names(mir$fits), info = paste(info0, "fits period names"))
  for (pk in pks) {
    fo <- orig$fits[[pk]]; fm <- mir$fits[[pk]]
    expect_equal(fo$h, fm$h, tolerance = 1e-8, info = paste(info0, "period", pk, "h"))
    expect_equal(fo$b, fm$b, tolerance = 1e-8, info = paste(info0, "period", pk, "b"))
  }

  expect_identical(orig$scheme_detected, mir$scheme_detected,
                    info = paste(info0, "scheme_detected"))
}

for (cfg in rddid_configs) {
  for (bw in c("cct", "joint", "iter")) {
    test_that(sprintf("rddid() mirror invariance: %s, bwselect=%s", cfg$name, bw), {
      check_rddid_invariance(cfg, bw)
    })
  }
}

test_that("rddid() mirror invariance: dgp_a, weights=linear, bwselect=joint", {
  check_rddid_invariance(rddid_configs[[1]], "joint", weights = "linear")
})

# =============================================================================
# B. rd_typecont()
# =============================================================================
# Per S1 of the plan, on mirrored data the type-continuity jump is IDENTICAL
# (not negated): the outcome 1 - 1{V_s = 1} on -R gives jump
# (1 - beta_hat_-) - (1 - beta_hat_+) = D_hat, same SE. The returned object
# exposes no raw per-(period,type) jump/SE directly; the per-period LL-Wald
# (`per_period[[k]]$ll_wald`, stat/df/p) is the only per-period numeric summary
# it carries (with a single kept binary-type contrast per period here, this is
# the squared standardized jump), so that is what we compare, alongside the
# joint `ll_wald`.

for (bc_val in c(TRUE, FALSE)) {
  test_that(sprintf("rd_typecont() mirror invariance: S1, bc=%s", bc_val), {
    info0 <- sprintf("S1 bwselect=cct bc=%s", bc_val)

    orig <- rd_typecont(S1, x = "R", time = "t", id = "id",
                         bwselect = "cct", bc = bc_val)
    mir  <- rd_typecont(mirror_x(S1), x = "R", time = "t", id = "id",
                         bwselect = "cct", bc = bc_val)

    expect_equal(orig$ll_wald$stat, mir$ll_wald$stat, tolerance = 1e-10,
                 info = paste(info0, "ll_wald$stat"))
    expect_equal(orig$ll_wald$df, mir$ll_wald$df, tolerance = 1e-10,
                 info = paste(info0, "ll_wald$df"))
    expect_equal(orig$ll_wald$p, mir$ll_wald$p, tolerance = 1e-10,
                 info = paste(info0, "ll_wald$p"))

    pks <- names(orig$per_period)
    expect_identical(pks, names(mir$per_period), info = paste(info0, "per_period names"))
    for (pk in pks) {
      lo <- orig$per_period[[pk]]$ll_wald
      lm <- mir$per_period[[pk]]$ll_wald
      expect_equal(lo$stat, lm$stat, tolerance = 1e-10,
                   info = paste(info0, "period", pk, "ll_wald$stat"))
      expect_equal(lo$df, lm$df, tolerance = 1e-10,
                   info = paste(info0, "period", pk, "ll_wald$df"))
      expect_equal(lo$p, lm$p, tolerance = 1e-10,
                   info = paste(info0, "period", pk, "ll_wald$p"))
    }
  })
}

# =============================================================================
# C. rd_homog() and rd_trendcell()
# =============================================================================
# On mirrored data the per-(period, type) jumps are negated with type labels
# swapped ("+" <-> "-"); joint statistics are unchanged. Match rows by
# (period, swapped type).

swap_side <- function(v) ifelse(v == "+", "-", ifelse(v == "-", "+", NA_character_))

# Assert jump (negated) and SE (unchanged) correspondence between an original
# and mirrored per-(period, type/cell) jump table.
check_jump_swap <- function(df_orig, df_mirror, type_col, info0, tol = 1e-10) {
  expect_true(nrow(df_orig) > 0, info = paste(info0, "orig jump table non-empty"))
  for (i in seq_len(nrow(df_orig))) {
    ro <- df_orig[i, ]
    match_row <- df_mirror[df_mirror$period == ro$period &
                            df_mirror[[type_col]] == swap_side(ro[[type_col]]), ]
    expect_equal(nrow(match_row), 1L,
                 info = paste(info0, "period", ro$period, type_col, ro[[type_col]],
                              "-> exactly one swapped-type match in mirrored table"))
    if (nrow(match_row) == 1L) {
      expect_equal(ro$jump, -match_row$jump, tolerance = tol,
                   info = paste(info0, "period", ro$period, type_col, ro[[type_col]],
                                "jump negated"))
      expect_equal(ro$se, match_row$se, tolerance = tol,
                   info = paste(info0, "period", ro$period, type_col, ro[[type_col]],
                                "se unchanged"))
    }
  }
}

test_that("rd_homog() mirror invariance: S0_3, scheme=pc, type_by=rd_side", {
  info0 <- "S0_3 t_rd=3 comparisons=1,2 scheme=pc type_by=rd_side"

  orig <- rd_homog(S0_3, y = "Y", x = "R", time = "t", id = "id",
                    t_rd = 3, comparisons = c(1, 2),
                    type_by = "rd_side", scheme = "pc")
  mir  <- rd_homog(mirror_x(S0_3), y = "Y", x = "R", time = "t", id = "id",
                    t_rd = 3, comparisons = c(1, 2),
                    type_by = "rd_side", scheme = "pc")

  expect_equal(orig$statistic, mir$statistic, tolerance = 1e-10, info = paste(info0, "statistic"))
  expect_equal(orig$df, mir$df, tolerance = 1e-10, info = paste(info0, "df"))
  expect_equal(orig$p_value, mir$p_value, tolerance = 1e-10, info = paste(info0, "p_value"))

  check_jump_swap(orig$period_type_jumps, mir$period_type_jumps, "type", info0)
})

test_that("rd_homog() mirror invariance: S1, scheme=pv, t_rd=2, comparisons=1", {
  info0 <- "S1 t_rd=2 comparisons=1 scheme=pv type_by=rd_side"

  orig <- rd_homog(S1, y = "Y", x = "R", time = "t", id = "id",
                    t_rd = 2, comparisons = 1,
                    type_by = "rd_side", scheme = "pv")
  mir  <- rd_homog(mirror_x(S1), y = "Y", x = "R", time = "t", id = "id",
                    t_rd = 2, comparisons = 1,
                    type_by = "rd_side", scheme = "pv")

  expect_equal(orig$statistic, mir$statistic, tolerance = 1e-10, info = paste(info0, "statistic"))
  expect_equal(orig$df, mir$df, tolerance = 1e-10, info = paste(info0, "df"))
  expect_equal(orig$p_value, mir$p_value, tolerance = 1e-10, info = paste(info0, "p_value"))

  check_jump_swap(orig$period_type_jumps, mir$period_type_jumps, "type", info0)
})

test_that("rd_trendcell() mirror invariance: S0_3, scheme=pc, type_by=rd_side", {
  info0 <- "S0_3 t_rd=3 comparisons=1,2 scheme=pc type_by=rd_side"

  orig <- rd_trendcell(S0_3, y = "Y", x = "R", time = "t", id = "id",
                        t_rd = 3, comparisons = c(1, 2),
                        type_by = "rd_side", scheme = "pc")
  mir  <- rd_trendcell(mirror_x(S0_3), y = "Y", x = "R", time = "t", id = "id",
                        t_rd = 3, comparisons = c(1, 2),
                        type_by = "rd_side", scheme = "pc")

  expect_equal(orig$statistic, mir$statistic, tolerance = 1e-10, info = paste(info0, "statistic"))
  expect_equal(orig$df, mir$df, tolerance = 1e-10, info = paste(info0, "df"))
  expect_equal(orig$p_value, mir$p_value, tolerance = 1e-10, info = paste(info0, "p_value"))

  check_jump_swap(orig$cell_period_jumps, mir$cell_period_jumps, "cell", info0)
})
