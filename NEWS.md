# rddid 0.2.0.9000 (development)

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
