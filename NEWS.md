# rddid 0.4.0.9000 (development)

* **Site rebuilt around three vignettes** (`vignettes/`, plan in `dev/site_plan.md`): *Estimating
  the ATT with rddid()* (per-period `rd_bw_cct()` + `rd_period()`, then `rddid()` under constant,
  linear and custom weights), *Bandwidth rules and sampling schemes* (`bwselect = "cct" / "joint" /
  "iter"`, `start`, fixed `h`; `scheme` detection and the three standard errors) and *Composition
  validation tests* (`rd_typecont()`, `rd_compstable()`, `rd_homog()`, `rd_trendcell()` on
  simulated scenarios that satisfy or violate each assumption). Every vignette hand-codes its DGP
  and prints the truth next to each estimate; `dev/site_dgp_check.R` verifies the scenario claims.
  `_pkgdown.yml` groups the reference into *Estimation* and *Validation tests*; README gains a
  quick start (now generated from `README.Rmd`).

* The composition-adjusted estimators (`rd_adjust()`, `rd_sadjust()`, `rd_c()`, `rd_att()`) are
  tagged `@keywords internal`: still exported and tested, but no longer listed on the site. They
  belong to a companion paper and are not part of the current manuscript.

* `Suggests` gains `knitr`, `rmarkdown`, `ggplot2`; `VignetteBuilder: knitr`.

# rddid 0.3.5.9000 (development)

* **Synced with Appendix B of the paper (rewritten 2026-09-08).** A code <-> equation map now
  lives in `dev/appB_map.md`: one row per object the package computes, with the paper's exact
  expression, the implementing symbol, conventions, audit status, and the test that pins it.
  `dev/check_appB_labels.R` verifies that every cited paper label still exists in `main.tex`;
  `dev/snapshot_rddid.R` is a numerical regression snapshot.

* Bandwidth selectors follow App. B.4 at general polynomial order `p` (previously the
  exponents and constants were hard-coded for `p = 1`): `rd_period()$b_const` is
  `(p+1)! (D - D_bc) / h^{p+1}`; the common `h*` (`eq:common_h_opt`) uses the constant
  `((p+1)!)^2 / (2(p+1))` and exponent `1/(2p+3)`; the period-specific objective and its
  regularization use `h^{p+1}/(p+1)!`. Numerically identical at `p = 1`.

* `bwselect = "iter"` under `scheme = "pc"`: the same-side cross-period covariance term of the
  aggregate AMSE now scales as `omega(h_t/h_s)/h_s` per side, with `omega(rho)` the kernel
  constant of Lemma `cov-pc` (new internal module `R/kernel_constants.R`, ported from the paper's
  simulation-verification code), instead of the previous `1/max(h_t, h_s)` approximation (exact
  only for the uniform kernel at `p = 0`). `.bw_joint_iter()` additionally returns the objective
  value and the objective function for diagnostics.

* `rd_period()`: the active set is the union of the pilot and main windows (was the pilot window
  alone, which silently truncated the main fit when `b < h`). Identical whenever `b >= h`.

* New tests: `test-appB-conformance.R` (from-the-equations reference implementation of the
  estimator, bias correction and the three sampling-scheme variances, 1e-10),
  `test-appB-bandwidth.R` (B.4 objectives and selectors), `test-kernel-constants.R`.

* Regularization of the bandwidth selectors: the curvature-variance estimate `Var(B-hat_t)`
  now pairs the pilot-window influence weights with the residuals of the order-`q` pilot fit at
  `b` (the same convention as the bias-corrected variance) instead of the order-`p` residuals at
  `h`. Moves `bwselect = "joint"`/`"iter"` bandwidths by a fraction of a percent at the CCT pilot
  ratio; at wide pilots (`b/h >= 3`) the old estimate over-stated the variance.

* Validation tests synced with Section 4.4 of the paper (`dev/tests_map.md`): documentation
  cites the assumptions by label (`ass:type-cont`, `ass:comp-stable`, `ass:homog`,
  `ass:trend-cell`) instead of numbers that changed in the paper; `rd_compstable()` drops one
  reference type instead of pseudo-inverting the structurally singular covariance (same
  statistic for binary types, exact df); a comparison-period unit exactly at the cutoff now
  stays on the reflected side; the joint-over-pairs result is documented as approximate. New
  from-the-text conformance tests `test-s44-conformance.R`.

* Guards: `bwselect = "cct"` and the joint pilot go through the guarded `rd_bw_cct()`
  (fallback + finiteness checks) instead of calling `rdrobust::rdbwselect()` directly;
  `.bw_joint()` warns when `h*` exceeds the running-variable radius; the coordinate descent
  warns when a bandwidth ends on the search boundary; `rddid(scheme = "pc"/"pv")` warns when
  no unit id repeats across periods (all cross-period covariances are then zero); the search
  cap is NA-safe.

* Stale references in code comments to the removed "Appendix C", `lem:coercive` and
  `eq:amse-ps` replaced by the current labels (`app:est-bw`, `eq:amse-att`, `eq:update`,
  `alg:coorddesc`).

# rddid 0.3.0.9000 (development)

* `rddid(..., bwselect = "iter")` gains a `start` argument controlling the
  coordinate-descent seed: `"hstar"` (default) seeds all periods at the common
  joint-optimal bandwidth h*; `"cct"` seeds each period at its own CCT/IK pilot
  h; or supply a named numeric/list of per-period bandwidths for a manual seed.
  Each run weakly improves the joint AMSE over its own start; seeding from h*
  therefore weakly dominates the common-h rule. (Seed from `"cct"` when the
  per-period biases nearly cancel, where h* is inflated.)

* Bug fix (reproducibility): `rd_typecont()` called `set.seed(NULL)` immediately
  before its Canay-Kamat permutation, which **re-initialised** the RNG from
  system entropy and discarded any seed the caller had set — so the CK p-value
  changed on every call. Removed; the permutation now inherits the caller's RNG
  state, so `set.seed()` before the call makes the CK p reproducible.
  (`rd_compstable()` never had this and was already reproducible.)

* `rd_typecont()`: the Canay-Kamat permutation test now chooses the number of
  nearest observations per side, `q`, by the Canay & Kamat (2018) rule of thumb
  **by default** (`q = NULL`), per period. A fixed `q` over-rejects in finite
  samples when the type distribution varies steeply in the running variable at
  the cutoff; the rule of thumb shrinks `q` as that association strengthens.
  Pass an integer `q` to force a fixed value; the per-period `q` used is
  returned in `meta$q_used`.

* `rd_compstable()`: same change — its Canay-Kamat permutation shares the
  identical fixed-`q` exposure, so `q` now defaults to the rule of thumb
  (`q = NULL`), chosen per `(t_RD, t_0)` pair on the pooled reflected sample.
  The per-pair `q` used is returned in `meta$q_used` (and echoed in each
  `pairs[[...]]$q`). Pass an integer `q` to force a fixed value.

* Bug fix (numerical robustness): the joint Wald pseudo-inverse (`.joint_wald()`,
  used by `rd_typecont()` and `rd_compstable()`) now uses the `MASS::ginv`
  relative tolerance `sqrt(eps)*max(sv)`. The previous, tighter tolerance could
  leave the structural-zero singular value (the per-period type indicators sum
  to 1) just above the cut on some LAPACK builds, inflating the statistic into a
  platform-dependent false rejection. Caught by the new CI on Ubuntu.

* Continuous integration: added an `R-CMD-check` GitHub Action (standard
  multi-OS matrix) and R-CMD-check / MIT-license badges to the README. Removed
  a placeholder ORCID from `DESCRIPTION`.
* `R CMD check` is now clean (0 errors / 0 warnings). Fixes: replaced the
  `\insertCite{}` macros in `rd_compstable()` (Rdpack was not a dependency)
  with the plain-text citations already in the References; documented the
  `regularize`/`reg_const` arguments of `rddid()`; and dropped the
  `VignetteBuilder: knitr` declaration (and the `knitr`/`rmarkdown` Suggests)
  since the package ships no vignettes.

* Internal: consolidated duplicated logic into shared helpers. A new
  `.cov_scheme()` (the `cs`/`pc`/`pv` scheme-combine of `.cross_cov()`) now
  backs the cross-period covariance in `rd_typecont()`, `rd_compstable()`, and
  `rd_homog()`, replacing three near-identical inline copies (including the
  former `.cross_cov_homog()`/`.match_sum_homog()`). A new `.scheme_from_long()`
  primitive backs both `.detect_scheme()` and `rd_typecont()`'s scheme
  detection. No change in results (verified to machine precision).

# rddid 0.2.1

* Internal: removed a duplicate `.build_types()` (a second copy lived in
  `test_homog.R` and shadowed the canonical one in `test_helpers.R` at load
  time). There is now a single shared implementation used by `rd_typecont()`
  and `rd_homog()`.
* Convention: units exactly at the cutoff are now treated as above it
  (`V_i = 1{R_i >= c}`) everywhere, including sampling-scheme detection. A unit
  sitting on the cutoff no longer registers as a separate "side" and so cannot
  be misread as a side-switch.
* Internal: dropped the unused single-cell `.ck_perm()` helper (and its test).
  `rd_typecont()`/`rd_compstable()` use an inlined *joint* Canay–Kamat
  permutation with one shared per-period shuffle; the standalone helper was a
  dead parallel path.
* Documentation: corrected the assumption numbering in the Section 3.4 test
  functions to match the manuscript's `\begin{assumption}` ordering —
  `rd_typecont()` is Assumption A7 (was mislabelled A6), `rd_compstable()` is
  A8 (was A7), and `rd_homog()` is A9 (was A8). Affects titles, `print` output,
  and cross-references only; no behaviour change.
* Documentation: fixed an unmatched apostrophe in the `rd_homog()` `@examples`
  comment that was silently dropping the entire example from the rendered help
  page.

# rddid 0.2.0

* Added tests for the Section 3.4 identifying assumptions of the
  time-varying-running-variable design:
  * `rd_typecont()` — continuity of the type distribution (LL-Wald and
    Canay–Kamat permutation, both necessary & sufficient; McCrary within-type
    sufficient-not-necessary; McCrary pooled neither).
  * `rd_compstable()` — composition stability across periods, via the reflection
    construction (LL-Wald and Canay–Kamat permutation, both necessary &
    sufficient). The permutation uses the partially-overlapping-samples scheme
    for units above the cutoff in more than one period.
  * `rd_homog()` — type-homogeneous confounding, tested in comparison periods
    (only suggestive of the assumption at the RD period: neither necessary nor
    sufficient there).

# rddid 0.1.0

* From-scratch rewrite: per-period local-linear RD engine (`rd_period()`,
  validated against `rdrobust` to machine precision), the `rddid()` aggregate
  estimator with constant/linear/custom weights, CS/PC/PV sampling-scheme
  variances, and joint / CCT / period-specific bandwidth selection.
