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
