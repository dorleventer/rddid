# Package index

## Package overview

- [`rddid-package`](https://dorleventer.github.io/rddid/reference/rddid-package.md)
  : rddid: treatment effects at a cutoff shared with a confounding
  policy

## Estimate

The RD-DID estimator and what to do with a fit.

- [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) :
  Estimate the effect of a treatment at a cutoff shared with a
  confounding policy
- [`summary(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/summary.rddid.md)
  [`print(`*`<summary.rddid>`*`)`](https://dorleventer.github.io/rddid/reference/summary.rddid.md)
  : Summary of an RD-DID fit
- [`plot(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/plot.rddid.md)
  : Plot the per-period RD fits behind an RD-DID estimate
- [`coef(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-methods.md)
  [`confint(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-methods.md)
  [`nobs(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-methods.md)
  : Coefficients, confidence intervals and sample size of an RD-DID fit
- [`tidy(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  [`glance(`*`<rddid>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  [`tidy(`*`<rd_typecont>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  [`tidy(`*`<rd_compstable>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  [`tidy(`*`<rd_homog>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  [`tidy(`*`<rd_trendcell>`*`)`](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)
  : Tidy output for RD-DID fits and validation tests

## Check the assumptions

Tests of the identification assumptions when the running variable moves
over time, so that units can change side of the cutoff between periods.

- [`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)
  : Test of type continuity
- [`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)
  : Test of composition stability
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  : Test of homogeneous confounding
- [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  : Test of constant within-type confounding
- [`plot(`*`<rd_typecont>`*`)`](https://dorleventer.github.io/rddid/reference/plot.rd_typecont.md)
  : Plot a type-continuity test: one comparison period against one RD
  period
- [`plot(`*`<rd_compstable>`*`)`](https://dorleventer.github.io/rddid/reference/plot.rd_compstable.md)
  : Plot a composition-stability test: the reflected sample of one pair
  of periods
- [`plot(`*`<rd_homog>`*`)`](https://dorleventer.github.io/rddid/reference/plot.rd_homog.md)
  : Plot a homogeneous-confounding test: the confounding jump of each
  type, by comparison period
- [`plot(`*`<rd_trendcell>`*`)`](https://dorleventer.github.io/rddid/reference/plot.rd_trendcell.md)
  : Plot a constant-within-type-confounding test: each type's
  confounding jump over time
- [`plot_switchers()`](https://dorleventer.github.io/rddid/reference/plot_switchers.md)
  : Plot the switchers: the running variable in one period against
  another

## Building blocks

The local-linear RD in one period and its CCT bandwidths, for working
period by period.

- [`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
  : Local-linear RD in one period (building block)
- [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  : CCT bandwidths for one period (building block)

## Data

Simulated panels for trying things out.

- [`rddid_sim`](https://dorleventer.github.io/rddid/reference/rddid_sim.md)
  : Simulated RD-DID panel (running variable fixed over time)
- [`rddid_sim_pv`](https://dorleventer.github.io/rddid/reference/rddid_sim_pv.md)
  : Simulated RD-DID panel (running variable moves over time)
