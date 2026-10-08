# rddid: treatment effects at a cutoff shared with a confounding policy

A treatment of interest switches on at a cutoff of a running variable in
one period, the **RD period**. A **confounding policy** switches at the
same cutoff, in every period, so the jump in the outcome at the cutoff
in the RD period mixes the treatment effect with the **confounding
jump**. In the **comparison periods** the treatment of interest is
uniform at the cutoff (nobody treated, or everybody treated), so the
jump there *is* the confounding jump.
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
estimates the jump in every period by local-linear RD and subtracts a
weighted average of the comparison-period jumps from the RD-period jump.
How the weights are set is the **confounding-trend assumption**:
constant (equal weights) or linear in time.

## Workflow

1.  **Estimate.**
    `fit <- rddid(data, y = , x = , time = , id = , t_rd = )`; printing
    `fit` shows the estimate, its standard error and confidence
    interval, the comparison periods with their weights, the sampling
    scheme and the bandwidth. See
    [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md).

2.  **Look inside.** `summary(fit)` shows the jump in every period and
    the standard error under each sampling scheme
    ([`summary.rddid()`](https://dorleventer.github.io/rddid/reference/summary.rddid.md));
    [`coef()`](https://rdrr.io/r/stats/coef.html),
    [`confint()`](https://rdrr.io/r/stats/confint.html) and
    [`nobs()`](https://rdrr.io/r/stats/nobs.html) work as usual
    ([rddid-methods](https://dorleventer.github.io/rddid/reference/rddid-methods.md)),
    and `tidy()` and `glance()` feed table makers
    ([rddid-tidiers](https://dorleventer.github.io/rddid/reference/rddid-tidiers.md)).

3.  **Check the assumptions.** When the running variable moves over
    time, units can change side of the cutoff between periods, and four
    tests check the assumptions this adds: type continuity
    ([`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md)),
    composition stability
    ([`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)),
    homogeneous confounding
    ([`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md))
    and constant within-type confounding
    ([`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)).

The building blocks
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md)
(the local-linear RD in one period) and
[`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
(its CCT bandwidths) are exported for users who want to work period by
period.
[rddid_sim](https://dorleventer.github.io/rddid/reference/rddid_sim.md)
and
[rddid_sim_pv](https://dorleventer.github.io/rddid/reference/rddid_sim_pv.md)
are simulated panels for trying things out. `citation("rddid")` gives
the reference below.

## References

Leventer, D. and D. Nevo (2024). *Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data.*
arXiv:2408.05847. <https://arxiv.org/abs/2408.05847>

## See also

The package website, with a getting-started guide and articles on the
computation, the options and the tests:
<https://dorleventer.github.io/rddid/>

## Author

**Maintainer**: Dor Leventer <leventerdor@gmail.com>

Authors:

- Daniel Nevo
