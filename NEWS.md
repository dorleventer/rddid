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
