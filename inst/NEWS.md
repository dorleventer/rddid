# didisc 0.4.0.9006

* The signed weight of each period's jump in the estimate (+1 for the RD period, minus the weight
  for a comparison period) is now `fit$period_weights` (was `fit$coef`), and
  `fit$period_weights_bc` (was `fit$coef_bc`). In `summary()` the column is `weight` (was `coef`)
  and `weight_bc` (was `coef_bc`). `coef()` still returns the estimate. No aliases for the old
  names. Estimates are unchanged.

# didisc 0.4.0.9005

* The documentation and printed output call the design and its estimate "diff-in-disc"
  (difference-in-discontinuities), as the paper does, in place of "RD-DID". `print()` of a fit
  now starts "Diff-in-disc estimate of the ATT ...". Estimates are unchanged.

# didisc 0.4.0.9004

* The package is renamed from rddid to didisc, after the name the paper gives the design
  (difference-in-discontinuities) and its estimator (DiDisc). The main function is now
  `didisc()`, its class `"didisc"` (with `summary.didisc()`, `plot.didisc()` and the other
  methods), and the simulated data `didisc_sim` and `didisc_sim_pv`. The old names are gone, with
  no aliases. Estimates are unchanged.

# rddid 0.4.0.9003

First public version.

* `rddid()` estimates the effect of a treatment that switches on at a cutoff of a running
  variable in one period, the RD period, when a confounding policy switches at the same cutoff
  in every period. It takes the jump in the outcome at the cutoff in the RD period and subtracts
  a weighted average of the jumps in comparison periods, where the treatment of interest is
  uniform at the cutoff. The weights follow the confounding-trend assumption (the confounding
  jump is constant, or moves linearly, over time), or are the minimum-variance weights within
  it. Several RD periods can be averaged, with weights the user sets.
* Standard errors for a repeated cross-section and for a panel, with a running variable fixed or
  moving over time. Each fit reports the local-linear estimate with its conventional standard
  error and the bias-corrected estimate with its robust standard error. The bandwidth is one
  common bandwidth chosen for the RD-DID estimate, or each period's own CCT bandwidth. The
  common bandwidth is chosen as rdrobust chooses one for a single period: the variance, in a
  panel the covariances between periods, and the design part of the bias at one rule-of-thumb
  bandwidth shared by every period, the curvature from each period's pilot fit, and a
  regularization that includes the covariances of the bias estimates between periods.
* Tests of the identification assumptions when the running variable moves over time:
  `rd_typecont()`, `rd_compstable()`, `rd_homog()` and `rd_trendcell()`, or all four at once
  with `rd_checks()`.
* `summary()`, `plot()`, `tidy()` and `glance()` methods; the building blocks `rd_period()` and
  `rd_bw_cct()`; simulated data `rddid_sim` and `rddid_sim_pv`.
