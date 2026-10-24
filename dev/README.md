# dev/ — developer notes for rddid

**2026-10-08:** the composition-adjusted family (`rd_att`, `rd_sadjust`, `rd_c`, `rd_adjust`) was removed from the package (UX sweep, T1644); references to it in `appB_map.md`, `site_plan.md`, `atu_estimand_plan.md` are historical. Recover at tag `v0.4.0.9000-composition`.

Not part of the package build (`.Rbuildignore`).

**NEWS (2026-10-16).** The user-facing `NEWS.md` is `inst/NEWS.md`: it is installed with the package (`news(package = "rddid")`) and pkgdown builds no changelog page from it. It describes the package as released, not its development. The development log up to the first public version is `dev/news-history.md`; add development notes there.

**Numerical gate (UX sweep, 2026-10-08).** `dev/snapshot_all.R` records the full return object of every exported function on seeded data under every option (134 cells) to an `.rds`; `dev/snapshot_compare.R <baseline> <new>` requires every baseline leaf to be `identical()` (new leaves allowed, reworded error messages allowed). Current baseline: `dev/snapshots/baseline_bwrot.rds` (2026-10-24, 0.4.0.9002: the joint bandwidth's constants as in App. B.4 "Bandwidth estimation"; recorded on the branch and compared against a fresh snapshot of master 8300e37 — exactly the 31 `rddid(bwselect = "joint")` cells changed. Note: master 8300e37 itself no longer matched `baseline_tests_final.rds`, as PRs #2–#9 changed outputs and argument names without re-recording). Previous: `dev/snapshots/baseline_tests_final.rds` (2026-10-15, T1688: studentized sup-t default, `studentize`, `joint_tests` with normal and studentized references, `rd_checks()`, intercept weights kept by `rd_period()`; against `baseline_sizefix.rds` only the share tests' p-values and `joint_tests` changed, every Wald and estimation leaf identical). Previous: `dev/snapshots/baseline_sizefix.rds` (2026-10-15, T1678: the validation tests' size corrections; against `baseline_ux2.rds` only the four tests' default results moved — typecont/compstable `statistic`/`p_value` (sup-t), compstable `joint/ll_wald` (full covariance), homog/trendcell under the common-bandwidth default — and every estimation leaf and every explicit `bwselect = "cct"`/fixed-`h` test cell is identical). Earlier ones kept for history: `baseline_ux2.rds` (2026-10-08, after the review round: the no-switcher panel cells of the four tests now error); `baseline_bugfix.rds` (after the five bug fixes; differs from `baseline_step1.rds` only in the `rd_homog(p = 2)` cells). Run both before and after any change to `R/`.

Contents: `appB_map.md` and `tests_map.md` (code ↔ paper maps), `check_appB_labels.R`, `snapshot_rddid.R`, `site_plan.md` and `site_dgp_check.R` (pkgdown site rebuild, 2026-09), `atu_estimand_plan.md` (the `estimand` argument, 2026-09-10).

## Versions

Every pull request that changes the code (`R/`, `data/`, `src/`, `NAMESPACE`) raises the version
in `DESCRIPTION`. Between releases only the fourth component moves (0.4.0.9001 -> 0.4.0.9002):
run `usethis::use_version("dev")` or edit by hand. Change the first heading of `inst/NEWS.md`
to match (`# rddid 0.4.0.9002`). Documentation, tests, vignettes and CI need no bump.

The reason: outputs record the package version (the application scripts write it into their
`.tex` headers), so each version number should name one state of the code. The CI job
`version-bump` (`.github/workflows/version-bump.yaml`) enforces both rules and is a required
check on `master`. A release sets a three-part version (e.g. 0.5.0); the next change after it
moves to 0.5.0.9000.

## Pre-commit hook

After cloning, enable the doc-sync pre-commit hook (regenerates `man/*.Rd` from roxygen and blocks commits where the generated docs are stale — the mismatch that otherwise fails `R CMD check`):

``` sh
git config core.hooksPath .githooks
```
