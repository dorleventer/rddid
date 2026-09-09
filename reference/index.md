# Package index

## Estimation

The RD-DID estimator, its per-period local-linear RD engine, and the
single-period CCT bandwidth.

- [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) :
  RD-DID estimation and inference
- [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  : Single-period local-linear RD discontinuity
- [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  : CCT (MSE-optimal) bandwidth for a single local-linear RD

## Validation tests

Tests of the identification assumptions for a time-varying running
variable (Section 4.4 of the paper).

- [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  : Test the continuity of the type distribution
- [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  : Test composition stability across periods
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  : Test of homogeneous confounding
- [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  : Test of a constant within-type confounding discontinuity
