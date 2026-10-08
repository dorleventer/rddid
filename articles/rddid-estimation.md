# Get started: estimate an RD-DID effect

## The setting

A treatment of interest switches on at a cutoff of a running variable in
one period, the **RD period**. A **confounding policy** switches at the
same cutoff, in every period, so the jump in the outcome at the cutoff
in the RD period mixes the treatment effect with the **confounding
jump**. In the **comparison periods** the treatment of interest is
uniform at the cutoff (nobody treated, or everybody treated), so,
provided the treatment of interest has no anticipation or carry-over
effects there (which the paper assumes), the jump there *is* the
confounding jump.
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
estimates the jump in every period by local-linear RD and subtracts a
weighted average of the comparison-period jumps from the RD-period jump:
``` math
\mathrm{ATT}(t_{\mathrm{RD}}) = D_{t_{\mathrm{RD}}} - \sum_t w_t D_t ,
```
where $`t_{\mathrm{RD}}`$ is the RD period, $`D_t`$ is the jump in the
outcome at the cutoff in period $`t`$, the sum runs over the comparison
periods, and $`w_t`$ is the weight of comparison period $`t`$. How the
weights are set is the **confounding-trend assumption**: equal weights
if the confounding jump is the same in every period, a straight-line
extrapolation if it moves linearly in time.

## The data

``` r

library(rddid)
head(rddid_sim)
#>   id year          R V W           Y
#> 1  1    1 -0.7785626 0 0 -0.38918276
#> 2  2    1 -0.9019249 0 0 -0.75738418
#> 3  3    1  0.1621249 1 0  1.10233495
#> 4  4    1  0.1695525 1 0  1.50736287
#> 5  5    1 -0.7321141 0 0 -0.04729411
#> 6  6    1  0.8406841 1 0  1.59100391
```

`rddid_sim` is a simulated panel of 1,000 units observed in three years,
in long format: one row per unit and year.
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md) uses
four columns: the outcome `Y`, the running variable `R` (cutoff 0), the
period `year` and the unit `id`. The other two record the design: the
confounding policy `V` is on above the cutoff in every year and raises
`Y` by 0.5; the treatment of interest `W` is on above the cutoff in year
3 only and raises `Y` by 1. So year 3 is the RD period, years 1 and 2
are comparison periods, and **the true effect is 1**.

With your own data,
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
needs the same shape: a data frame in long format, one row per unit and
period, with the outcome, the running variable, the period and, for a
panel, a unit identifier (the panel need not be balanced). The cutoff is
`c` (default 0), and a unit with `x >= c` is above it. The period column
is usually numeric, a year; a character column works with the default
constant trend, while `trend = "linear"` fits its line on the period
values and so needs a numeric column. By default the comparison periods
are all periods other than `t_rd`, in sorted order.

## Estimate

``` r

fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3868, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

Line by line:

- **`RD-DID estimate of the ATT in period 3`**: the estimand and the RD
  period (`t_rd = 3`). The ATT is the effect of the treatment on the
  units just above the cutoff, which are the ones treated in the RD
  period. If everybody is treated in the comparison periods,
  `estimand = "atu"` labels the same numbers as the ATU (see [the ATU
  section](#atu) below).
- **`Comparison periods: 1, 2 (constant confounding trend; weights 0.5, 0.5)`**:
  by default every period other than `t_rd` is a comparison period. The
  default `trend = "constant"` assumes the confounding jump is the same
  in every period, so each comparison period gets the same weight and
  their average stands in for the confounding jump in year 3.
- **`Sampling scheme`**: how the data were sampled, read off `id` and
  `R`: here a panel (the same units every year) in which no unit changes
  side of the cutoff, because the running variable does not move. The
  standard error takes into account that the same units enter every
  year’s fit. Under the default bandwidth rule the scheme can also move
  the bandwidth, and so the estimate; [Bandwidth rules and sampling
  schemes](https://dorleventer.github.io/rddid/articles/rddid-options.md)
  shows how.
- **`Bandwidth: common h`**: the default rule, `bwselect = "joint"`,
  uses one main bandwidth `h` in every period, chosen to minimize the
  asymptotic mean squared error of the RD-DID estimate.
  **`Pilot bandwidth b (period = value)`** lists the bandwidth of the
  bias correction in each period, in time order.
- **The two rows.** `Conventional` is the local-linear estimate with its
  conventional standard error. `Robust (bias-corrected)` is the
  bias-corrected estimate with its robust standard error (Calonico,
  Cattaneo and Titiunik, 2014). Each row has its z statistic, p-value
  and 95% confidence interval.

The true effect is 1; the two estimates are 1.09 and 1.14. A plain RD in
year 3 alone would estimate the effect plus the confounding jump, 1 +
0.5 = 1.5. The year-3 fit is stored in the object:

``` r

fit$fits[["3"]]
#> Single-period RD (p=1, h=0.2672, b=0.3951, kernel=triangular, n=1000)
#>   D (conventional)   = +1.664  (se 0.1641)
#>   D (bias-corrected) = +1.6854  (se 0.1931)
```

Its jump, 1.66, is that naive RD estimate, whose target is 1.5, not 1.
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
subtracts the average of the year-1 and year-2 jumps, (0.51 + 0.63) / 2
= 0.57, its estimate of the confounding jump of 0.5.

## Summary, coefficients and tables

``` r

summary(fit)
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3868, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   Per-period local-linear fits, in time order (estimate = sum of coef x jump):
#>   period   role          coef      n        h        b       jump      s.e.  jump (bc) s.e. (rb)
#>   1        comparison    -0.5   1000   0.2672   0.4102     0.5119    0.1688     0.4639    0.1948
#>   2        comparison    -0.5   1000   0.2672   0.3868     0.6306    0.1686     0.6240    0.2012
#>   3        RD               1   1000   0.2672   0.3951     1.6640    0.1641     1.6854    0.1931
#> 
#>   Robust s.e. under each sampling scheme (the printed one is for "pc"; the others are for comparison):
#>     cross-section 0.2385   panel, no unit changes side 0.1494   panel, some change side 0.1494
```

[`summary()`](https://rdrr.io/r/base/summary.html) adds the per-period
fits. Each row gives the period’s role, its coefficient in the estimate
(+1 for the RD period, minus its weight for a comparison period), its
number of observations, its bandwidths, and its jump with standard
error, conventional (`jump`, `s.e.`) and bias-corrected (`jump (bc)`,
`s.e. (rb)`). The estimate is the sum of coefficient times jump: 1.664 -
0.5 × 0.512 - 0.5 × 0.631 = 1.093.

The last lines give the robust standard error under each of the three
sampling schemes. The printout uses the one for the scheme detected in
the data; the other two show what the formulas for the other schemes
give on the same data. Here the cross-section standard error, 0.239,
ignores that the same units appear in every year; the panel one, which
accounts for it, is 0.149.

``` r

coef(fit)
#> Conventional       Robust 
#>     1.092709     1.141380
confint(fit)
#>                  2.5 %   97.5 %
#> Conventional 0.8449613 1.340457
#> Robust       0.8485735 1.434187
generics::tidy(fit)
#>           term estimate std.error statistic      p.value  conf.low conf.high
#> 1 Conventional 1.092709 0.1264043  8.644556 5.401507e-18 0.8449613  1.340457
#> 2       Robust 1.141380 0.1493939  7.640071 2.171020e-14 0.8485735  1.434187
```

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the two point
estimates and [`confint()`](https://rdrr.io/r/stats/confint.html) their
confidence intervals. `tidy()` (from the generics package) returns one
row per estimate, so table makers such as modelsummary work with
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
fits; [A table for a paper](#a-table-for-a-paper) below has an example.

## A confounding jump that moves linearly

``` r

fit_lin <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3,
                 trend = "linear")
fit_lin
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (linear confounding trend; weights -1, 2)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3869, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9147     0.2582    3.54   <0.001   [0.4087, 1.4207]
#>   Robust (bias-corrected)      0.9012     0.3068    2.94    0.003   [0.3000, 1.5025]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

With `trend = "linear"` the confounding jump may move linearly in time.
The line through the year-1 and year-2 jumps is extrapolated to year 3,
which gives weights -1 and 2: the confounding jump in year 3 is
estimated by $`2 D_2 - D_1`$. This needs at least two comparison
periods. The standard error grows (0.258 against 0.126) because
extrapolating a line amplifies noise. In `rddid_sim` the confounding
jump is constant, so both estimates are valid; the linear one gives up
precision to allow for a linear trend. A numeric vector instead of
`"constant"` or `"linear"` sets the weights directly, one per comparison
period.

## A table for a paper

`tidy()` gives the estimates and `glance()` the design (RD period,
comparison periods, trend, bandwidth rule, `h`, scheme), which is what
table makers read. With the modelsummary package installed, the two fits
side by side, with confidence intervals:

``` r

modelsummary::modelsummary(
  list(Constant = fit, Linear = fit_lin), statistic = "conf.int",
  gof_map = list(list(raw = "h", clean = "h", fmt = 3),
                 list(raw = "nobs", clean = "Num.Obs.", fmt = 0)))
```

|              | Constant         | Linear           |
|--------------|------------------|------------------|
| Conventional | 1.093            | 0.915            |
|              | \[0.845, 1.340\] | \[0.409, 1.421\] |
| Robust       | 1.141            | 0.901            |
|              | \[0.849, 1.434\] | \[0.300, 1.502\] |
| h            | 0.267            | 0.267            |
| Num.Obs.     | 3000             | 3000             |

The same call with `output = "rddid_table.tex"` (or `.html`, `.md`,
`.docx`) writes the table to a file;
`write.csv(generics::tidy(fit), "rddid_estimate.csv", row.names = FALSE)`
writes the estimate table to a plain CSV file.

`gof_map` sets which rows of `glance()` appear below the estimates and
their digits (`fmt`): here `h` with three digits and the number of
observations. Without it, every column of `glance()` is listed and `h`
prints with all its digits. The last line writes the estimate table to a
plain CSV file.

## When the running variable moves over time

In `rddid_sim` every unit stays on its side of the cutoff. In many
applications the running variable moves (a population count, a test
score, an income), and some units, the **switchers**, are above the
cutoff in one period and below it in another. The units at the cutoff in
the RD period can then be a different mix from the units at the cutoff
in a comparison period, and the difference of jumps can be biased. A
unit’s **type** records its side of the cutoff in periods other than the
one at hand. The tests use two versions:

| Test | A unit’s type in period *t* | Printed as |
|----|----|----|
| type continuity, composition stability | its sides of the cutoff in all the other periods, in time order | `++`, `+-`, `-+`, `--` (`+` above, `-` below) |
| homogeneous confounding, constant within-type confounding | its side of the cutoff in the RD period (the default, `type_by = "rd_side"`) | `+`, `-`; in plots “Above in 3”, “Below in 3” |

With two periods, both versions are the unit’s side of the cutoff in the
other period.

The RD-DID estimate identifies the ATT when the type distribution is
continuous at the cutoff and the confounding jump within each type is
constant over time, and in addition either the type composition is
stable across periods or the confounding jump is the same for every type
(Leventer and Nevo, 2024). Four tests check these assumptions.
`rddid_sim_pv` is the `rddid_sim` design with a running variable that
drifts between years. The estimate whose assumptions the four tests
below check:

``` r

fit_pv <- rddid(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
fit_pv
#> RD-DID estimate of the ATT in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Bandwidth: common h = 0.2367 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.3889, 2 = 0.3847, 3 = 0.3933
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 0.9654     0.1933    4.99   <0.001   [0.5866, 1.3442]
#>   Robust (bias-corrected)      0.9823     0.2292    4.29   <0.001   [0.5332, 1.4315]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The sampling scheme is now a panel in which some units change side of
the cutoff between years.

**Type continuity.** H0: the share of each type jumps by zero at the
cutoff, in every period. A rejection means units sort around the cutoff
by their side in other periods, so the jump in the outcome partly
reflects who the units are; the estimate can then be biased, whatever
the other tests say.

``` r

tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
tc
#> Test of a continuous type distribution  [rd_typecont()]
#>   H0: the share of each type jumps by zero at the cutoff, in every period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Periods: 1, 2, 3   Types: ++, +-, -+, --   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Joint Wald chi-squared(9) = 4.782,  p = 0.853
#>     Period 1: chi-squared(3) = 1.531,  p = 0.675
#>     Period 2: chi-squared(3) = 0.382,  p = 0.944
#>     Period 3: chi-squared(3) = 2.202,  p = 0.532
```

Every test prints its null (H0) and the sampling scheme read off the
data. `Types: ++, +-, -+, --` are the types of the first row of the
table above. `Bandwidth: CCT MSE-optimal, chosen per cell` means that
each local-linear fit the test runs (one per period and type, a *cell*)
gets its own MSE-optimal bandwidth from
[`rdrobust::rdbwselect`](https://rdrr.io/pkg/rdrobust/man/rdbwselect.html)
(Calonico, Cattaneo and Titiunik, 2014; hence “CCT”).
`Joint Wald chi-squared(9)` tests all the jumps in the type shares at
once: in each year the four shares sum to one, so three jumps are free,
and three years give 9. Each `Period` line tests one year alone.

**Composition stability.** H0: the share of each type among the units
just above the cutoff is the same in the RD period and in each
comparison period. A rejection means the units at the cutoff in the RD
period are a different mix of types from those in a comparison period;
the estimate is then valid only if the next test’s null holds.

``` r

cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
cs
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just above the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: chi-squared(3) = 21.982,  p = <0.001
#>     n above the cutoff: 490 (RD period), 508 (comparison), 430 in both
#>   Pair 3::2: chi-squared(3) = 6.137,  p = 0.105
#>     n above the cutoff: 490 (RD period), 497 (comparison), 418 in both
#> 
#>   Joint over pairs (sum of chi-squared): chi-squared(6) = 28.119,  p = <0.001
```

`Pair 3::1` compares year 3, the RD period, with year 1; its
`chi-squared(3)` tests the three free jumps in the type shares.
`n above the cutoff` counts the units above the cutoff in the RD year,
in the comparison year, and in both. `Joint over pairs` adds up the pair
statistics and their degrees of freedom; it treats the pairs as
independent although they share the year-3 units, so its p-value is
approximate (the paper’s test is per pair).

**Homogeneous confounding.** H0: in each comparison period the
confounding jump is the same for every type. Here, and in the next test,
the types are those of the second row of the table above: a unit’s side
of the cutoff in the RD period. A rejection means the confounding jump
differs by type; together with a rejection of composition stability, the
comparison-period jumps do not measure the confounding jump of the units
at the RD-period cutoff, and the estimate can be biased.

``` r

hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
hg
#> Test of homogeneous confounding  [rd_homog()]
#>   H0: in each comparison period the confounding jump is the same for every type
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Comparison periods: 1, 2
#> 
#>   Wald chi-squared(2) = 0.465,  p = 0.793
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Period     Type               jump       s.e.       n
#>     1          -                0.8985     0.2213     510  (reference)
#>     1          +                0.7374     0.3582     490
#>     2          -                1.0282     0.4007     510  (reference)
#>     2          +                0.7356     0.3588     490
```

The table gives the confounding jump of each type in each comparison
year, with its standard error and number of observations. The Wald
statistic has one contrast per year, type `+` against the type marked
`(reference)`, so two comparison years give `chi-squared(2)`.

**Constant within-type confounding.** H0: within each type, the
confounding jump is the same in every comparison period. This is a
pre-trends check: the assumption concerns the RD period, where the
confounding jump cannot be seen apart from the effect, so the test asks
whether it is stable across the comparison periods. A rejection means
the constant-trend weights do not cancel the confounding jump.
`trend = "linear"` allows a confounding jump that moves linearly in
time; `rd_trendcell(trend = "linear")` tests that version (three or more
comparison periods).

``` r

tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
tr
#> Test of a constant within-type confounding discontinuity  [rd_trendcell()]
#>   H0: within each type, the confounding jump is the same in every comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Comparison periods: 1, 2   Trend: constant
#> 
#>   Wald chi-squared(2) = 0.078,  p = 0.962
#> 
#>   Per-cell local-linear jumps (comparison periods):
#>     Type       Period             jump       s.e.       n
#>     +          1                0.7374     0.3582     490  (reference)
#>     +          2                0.7356     0.3588     490
#>     -          1                0.8985     0.2213     510  (reference)
#>     -          2                1.0282     0.4007     510
```

The jumps are those of
[`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
grouped by type. The Wald statistic has one contrast per type: year 2
against the year marked `(reference)`, the type’s first comparison year.

One table for all four:

``` r

do.call(rbind, lapply(list(tc, cs, hg, tr), generics::tidy))
#>                               test   statistic df      p.value
#> 1                  type continuity  4.78163259  9 8.529133e-01
#> 2            composition stability 28.11902600  6 8.923469e-05
#> 3          homogeneous confounding  0.46500200  2 7.925490e-01
#> 4 constant within-type confounding  0.07832992  2 9.615921e-01
```

The composition-stability row is the joint-over-pairs statistic, whose
p-value is approximate; the per-pair tests are in its printout above.

In `rddid_sim_pv` three of the four nulls are true by construction:
types are continuous at the cutoff, and the confounding jump is 0.5 for
every unit in every year. Composition stability is false: units at the
year-3 cutoff were mostly on the same side in years 1 and 2, while units
at the year-1 cutoff are spread evenly over the four types. In this
sample the four p-values are 0.853, below 0.001, 0.793 and 0.962 (in the
order of the table). A p-value below 0.05 for a true null is a false
rejection, which a test at the 5% level is designed to make in about one
sample in twenty. Because the other three assumptions hold, homogeneous
confounding among them, the failure of composition stability does not
bias `fit_pv`: its estimate is 0.97, against a true effect of 1. The
article on the identification assumptions shows how each statistic is
built and the options of the tests.

## Comparison periods where everybody is treated: the ATU

When everybody is treated in the comparison periods, rather than nobody,
the same difference of jumps estimates the ATU: the effect for the units
just below the cutoff, which are untreated in the RD period.
`estimand = "atu"` labels the output as the ATU: the estimate, standard
errors and bandwidths are the same as without it (among the tests, it
changes only
[`rd_compstable()`](https://dorleventer.github.io/rddid/reference/rd_compstable.md)).
For illustration, make everybody treated in years 1 and 2 of
`rddid_sim`, with the same effect of 1:

``` r

sim_atu <- transform(rddid_sim, W = ifelse(year < 3, 1, W), Y = Y + (year < 3))
rddid(sim_atu, y = "Y", x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu")
#> RD-DID estimate of the ATU in period 3
#>   Comparison periods: 1, 2   (constant confounding trend; weights 0.5, 0.5)
#>   Sampling scheme: panel, no unit changes side of the cutoff (detected)
#>   Bandwidth: common h = 0.2672 (rule "joint": one bandwidth, chosen for the RD-DID estimate)
#>   Pilot bandwidth b (period = value): 1 = 0.4102, 2 = 0.3868, 3 = 0.3951
#> 
#>                              Estimate  Std. err.       z  p-value   95% CI
#>   Conventional                 1.0927     0.1264    8.64   <0.001   [0.8450, 1.3405]
#>   Robust (bias-corrected)      1.1414     0.1494    7.64   <0.001   [0.8486, 1.4342]
#> 
#>   summary() shows the per-period fits and the s.e. under every sampling scheme.
```

The numbers are those of `fit`, now labeled the ATU; the true ATU here
is 1.

## What the paper reports

- **The estimate.** The paper reports both rows. It bases its confidence
  intervals and p-values on the `Robust (bias-corrected)` row, following
  Calonico, Cattaneo and Titiunik (2014): at a bandwidth chosen to
  minimize the mean squared error, the conventional interval is centered
  on a biased estimate and can under-cover.
- **The design.** The choices behind an estimate (the RD period, the
  comparison periods, the trend assumption with its weights, the
  bandwidth rule and `h`, the sampling scheme) are in the first lines of
  the printout; `generics::glance(fit)` returns them as a one-row table.
- **The tests.** When the running variable moves over time, the table
  above collects the four tests of the assumptions the estimate relies
  on.

## Next

- [Checking the identification
  assumptions](https://dorleventer.github.io/rddid/articles/rddid-validation-tests.md):
  how each test statistic is built, the tests under the ATU design, and
  their options.
- [Plots of the estimate and the
  checks](https://dorleventer.github.io/rddid/articles/rddid-plots.md):
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit and
  on each test, the switchers picture, and saving a figure.
- [Bandwidth rules and sampling
  schemes](https://dorleventer.github.io/rddid/articles/rddid-options.md):
  `bwselect`, fixed bandwidths, and the standard error under each
  scheme.
- [How rddid() computes the
  estimate](https://dorleventer.github.io/rddid/articles/rddid-how-it-works.md):
  the estimate and the tests rebuilt by hand, for referees and readers
  who want to see the pieces.

## References

`citation("rddid")` returns the reference for the package, the paper
below, with a BibTeX entry.

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
*Econometrica* 82(6), 2295–2326.
