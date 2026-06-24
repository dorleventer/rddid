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
