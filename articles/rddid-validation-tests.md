# Checking the identification assumptions

[Get
started](https://dorleventer.github.io/rddid/articles/rddid-estimation.md)
runs the four tests and reads their printouts. This article adds how
each statistic is built, the tests under the ATU design, their options,
and how the four results read together.

## When the tests apply

The tests are for a panel in which the running variable moves over time,
so that some units, the **switchers**, are above the cutoff in one
period and below it in another. The units at the cutoff in the RD period
can then be a different mix from the units at the cutoff in a comparison
period, and the difference of the two jumps can mix the treatment effect
with that difference in who is at the cutoff. A unit’s **type** records
its side of the cutoff in other periods; the table in [Get
started](https://dorleventer.github.io/rddid/articles/rddid-estimation.html#when-the-running-variable-moves-over-time)
gives the two versions the tests use. The RD-DID estimate identifies the
ATT when

1.  **type continuity** holds,
2.  **constant within-type confounding** holds, and
3.  **composition stability** or **homogeneous confounding** holds (one
    of the two is enough).

This is the paper’s result for a time-varying running variable (Leventer
and Nevo, 2024). There is one test per assumption. When the running
variable does not move, as in `rddid_sim`, every unit has the same type
in every period and the tests are not informative; in a repeated
cross-section the types are not observed.

The examples use `rddid_sim_pv`: the `rddid_sim` design (a confounding
jump of 0.5 in every year, a treatment effect of 1 in year 3, the RD
period) with a running variable that drifts between years. Most units
stay on one side; some switch:

``` r

library(rddid)
above <- rddid_sim_pv$R >= 0
switcher <- tapply(above, rddid_sim_pv$id, function(a) length(unique(a)) > 1)
sum(switcher)
#> [1] 176
```

Every test takes the same long data frame and column names as
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md). By
default each test uses the bias-corrected jumps with their robust
variances (`bc = TRUE`), like the `Robust (bias-corrected)` row of
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md), and
picks a separate CCT bandwidth for every local-linear fit it runs
(`bwselect = "cct"`).

## Type continuity: `rd_typecont()`

**H0: the share of each type jumps by zero at the cutoff, in every
period.**

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

- `Types: ++, +-, -+, --`: a unit’s sides of the cutoff in the other two
  years, in time order (the first row of the type table in Get started).
- The test runs, in each year and for each type, a local-linear RD of
  the indicator “the unit is of this type” on the running variable; the
  jump is the difference between the type’s share just above and just
  below the cutoff.
- `Joint Wald chi-squared(9)`: in each year the four shares sum to one,
  so three jumps are free; three years give 9. The `Period` lines test
  each year alone.

Here p = 0.853. Type continuity holds in `rddid_sim_pv` by construction,
so a p-value below 0.05 would be a false rejection.

**If it rejects:** units sort around the cutoff according to their side
in other periods, so the jump in the outcome partly reflects a jump in
who the units are. The RD-DID estimate can then be biased, whatever the
other tests say.

## Composition stability: `rd_compstable()`

**H0: the share of each type among the units just above the cutoff is
the same in the RD period and in each comparison period.**

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

- One block per pair of the RD period and a comparison period:
  `Pair 3::1` compares year 3 with year 1. The types are those of type
  continuity (a unit’s sides in the other two years), so there are four
  types and three free shares (`chi-squared(3)`).
- The test stacks the units above the cutoff in the two years, with the
  comparison year’s running variable reflected below an artificial
  cutoff, and tests the jump in each type share there.
- `n above the cutoff`: the units above the cutoff in the RD year, in
  the comparison year, and in both. Units above the cutoff in both years
  appear on both sides of the stacked regression; the standard error
  accounts for them.
- `Joint over pairs`: the sum of the pair statistics and degrees of
  freedom. It treats the pairs as independent although they share the
  RD-period units, so its p-value is approximate; the paper’s test is
  per pair.

Here p \< 0.001. Composition stability is false in `rddid_sim_pv` by
construction: the units at the year-3 cutoff were mostly on the same
side in years 1 and 2, while the units at the year-1 cutoff are spread
evenly over the four types. A rejection here is a correct detection;
with 1,000 units the test does not always detect it.

**If it rejects:** the units at the cutoff in the RD period are a
different mix of types from those in a comparison period. If the
confounding jump also differs across types, the comparison-period jump
measures the confounding jump of a different mix, and the estimate can
be biased. If the confounding jump is the same for every type (the next
test), the mix does not matter and the estimate is still valid.

## Homogeneous confounding: `rd_homog()`

**H0: in each comparison period the confounding jump is the same for
every type.**

Here a unit’s type is its side of the cutoff in the RD period (the
second row of the type table in Get started; `type_by = "rd_side"`, the
default).

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

- The table gives, in each comparison year, the jump in the outcome at
  the cutoff among the units below (`-`) and above (`+`) the cutoff in
  year 3. In a comparison period the jump is the confounding jump, so
  these are the confounding jumps by type.
- `Wald chi-squared(2)`: one contrast (type `+` against the reference
  type `-`) per comparison year.

Here p = 0.793. The null is true in `rddid_sim_pv` (the confounding jump
is 0.5 for every unit), so a p-value below 0.05 would be a false
rejection. The per-type jumps are imprecise (standard errors 0.22 to
0.40).

**If it rejects:** the confounding jump differs across types. If
composition stability holds, a confounding jump that differs across
types does not bias the estimate. If both nulls are false, the estimate
can be biased.

## Constant within-type confounding: `rd_trendcell()`

**H0: within each type, the confounding jump is the same in every
comparison period.**

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

- The per-cell jumps are those of
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md),
  compared the other way:
  [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  compares types within a comparison period,
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md)
  compares comparison periods within a type.
- `Wald chi-squared(2)`: one contrast (year 2 against year 1) per type.

This is a pre-trends check. The assumption concerns the RD period, where
the confounding jump cannot be separated from the effect, so the test
asks whether the within-type confounding jump is stable across the
comparison periods. Not rejecting does not establish the assumption in
the RD period. Here p = 0.962.

With `trend = "linear"` the null is that the within-type confounding
jump moves linearly across the comparison periods. That needs at least
three comparison periods; with two, the function says so and returns no
statistic.

**If it rejects:** the within-type confounding jump changes over time,
and constant weights do not cancel it. `rddid(trend = "linear")` allows
a within-type confounding jump that moves linearly in time, and
`rd_trendcell(trend = "linear")` tests that version (three or more
comparison periods).

## ATU designs

When everybody is treated in the comparison periods
(`estimand = "atu"`), the estimate concerns the units just below the
cutoff, and composition stability concerns the below-cutoff shares.
`rd_compstable(estimand = "atu")` mirrors the running variable and runs
the same test on the units below the cutoff; the other three tests give
the same numbers with either estimand, so `estimand` only labels them.
`rddid_sim_pv` is an ATT design, so the call below only shows the
mechanics:

``` r

rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3, estimand = "atu")
#> Test of composition stability  [rd_compstable()]
#>   H0: the share of each type among the units just below the cutoff is the same in the RD period and in each comparison period
#>   Sampling scheme: panel, some units change side of the cutoff (detected)
#>   Estimand: ATU (the units below the cutoff are the ones untreated in the RD period, so the test is on their shares (mirrored design))
#>   RD period: 3   Comparison periods: 1, 2   Bandwidth: CCT MSE-optimal, chosen per cell
#> 
#>   Pair 3::1: chi-squared(3) = 6.206,  p = 0.102
#>     n below the cutoff: 510 (RD period), 492 (comparison), 432 in both
#>   Pair 3::2: chi-squared(3) = 3.164,  p = 0.367
#>     n below the cutoff: 510 (RD period), 503 (comparison), 431 in both
#> 
#>   Joint over pairs (sum of chi-squared): chi-squared(6) = 9.370,  p = 0.154
```

In this printout the H0 line and the counts refer to the units below the
cutoff: under the ATU design those are the units untreated in the RD
period, whose composition the estimate relies on.

## Run all four and tabulate

``` r

tests <- list(
  rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3),
  rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3),
  rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3),
  rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
)
do.call(rbind, lapply(tests, tidy))
#>                               test   statistic df      p.value
#> 1                  type continuity  4.78163259  9 8.529133e-01
#> 2            composition stability 28.11902600  6 8.923469e-05
#> 3          homogeneous confounding  0.46500200  2 7.925490e-01
#> 4 constant within-type confounding  0.07832992  2 9.615921e-01
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns one
row per test, so the four stack into one table for a paper or a table
maker such as modelsummary. The composition-stability row is the
joint-over-pairs statistic, whose p-value is approximate.

Together, the four results map onto the list at the top. A rejection of
type continuity, or of constant within-type confounding, would call the
estimate into question; so would a rejection of both composition
stability and homogeneous confounding. A rejection of only one of those
two would not on its own, because either one is enough for
identification. In `rddid_sim_pv` the truth is known: composition
stability fails and the other three assumptions hold, so the estimate is
valid:

``` r

rddid(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
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

The true effect is 1.

## Options

- `bwselect`: `"cct"` (default) chooses a CCT bandwidth for every
  local-linear fit the test runs; `"rot"` uses a rule of thumb. `h`
  fixes the bandwidth.
- `bc`: `TRUE` (default) tests the bias-corrected jumps with their
  robust variances; `FALSE` uses the conventional ones.
- `scheme`: detected from the data as in
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md);
  set it to override.
- [`rd_homog()`](https://dorleventer.github.io/rddid/reference/rd_homog.md)
  and
  [`rd_trendcell()`](https://dorleventer.github.io/rddid/reference/rd_trendcell.md):
  `type_by = "rd_side"` (default, the side in the RD period) or
  `"pattern"` (the sides in all other periods); `min_n` drops cells with
  too few observations on either side of the cutoff, with a message.

The article [How rddid() computes the
estimate](https://dorleventer.github.io/rddid/articles/rddid-how-it-works.md)
rebuilds the type-continuity and composition-stability tests by hand and
runs all four on simulated designs in which one assumption fails at a
time, with the resulting bias in the estimate.

## Reference

Leventer, D. and D. Nevo (2024). Correcting Invalid Regression
Discontinuity Designs Using Multiple Time-Period Data. arXiv:2408.05847.
<https://arxiv.org/abs/2408.05847>
