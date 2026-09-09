# Test composition stability across periods

Wald test of the composition-stability assumption (Section 4.4 of
Leventer and Nevo): for each RD-period / comparison-period pair
\\(t\_{\mathrm{RD}}, t_0)\\, the share of each type among the units just
above the cutoff is the same in the two periods,
\\\pi\_{t\_{\mathrm{RD}},(+)}(v) = \pi\_{t_0,(+)}(v)\\. With two periods
the type is binary (the unit's side in the other period) and this is the
single share-jump test of the paper's Section 4.4; with more periods the
type \\(\mathbf{u}, b)\\ collects the unit's sides in the other
comparison periods, \\\mathbf{u}\\, and its side \\b\\ in the partner
period of the pair.

## Usage

``` r
rd_compstable(
  data,
  x,
  time,
  id,
  t_rd,
  comparisons = NULL,
  c = 0,
  h = NULL,
  bwselect = c("cct", "rot"),
  kernel = "triangular",
  scheme = c("auto", "cs", "pc", "pv"),
  bc = TRUE,
  ...
)
```

## Arguments

- data:

  a long data frame, one row per unit-period. A unit's type in period
  \\t\\ is read from its running variable in the other period(s); units
  unobserved there are dropped from period \\t\\, so the panel need not
  be balanced.

- x:

  Column name (string) for the running variable.

- time:

  Column name (string) for the period indicator.

- id:

  Column name (string) for the unit identifier.

- t_rd:

  Value of `time` identifying the RD period.

- comparisons:

  Values of `time` to use as comparison periods. If `NULL` (default),
  all periods except `t_rd` are used.

- c:

  Cutoff for the running variable (default 0).

- h:

  Bandwidth. If `NULL` (default), the bandwidth is determined by
  `bwselect`; an explicit numeric value overrides `bwselect` and is used
  directly.

- bwselect:

  Bandwidth selection rule when `h = NULL`: `"cct"` (default) computes a
  per-cell CCT MSE-optimal bandwidth via
  [`rd_bw_cct()`](https://dorleventer.github.io/rddid/reference/rd_bw_cct.md)
  for each type indicator RD in the reflected space; `"rot"` uses \\0.5
  \times \mathrm{IQR}(x)\\ as a rule-of-thumb applied to the full
  sample. Ignored when `h` is supplied explicitly.

- kernel:

  Kernel for the local-linear RD: `"triangular"` (default),
  `"epanechnikov"`, or `"uniform"`.

- scheme:

  Covariance scheme for the Wald test:

  `"auto"`

  :   Detects whether any unit appears on both sides of the artificial
      cutoff (i.e., above the true cutoff in both periods). If yes, uses
      `"pv"` (time-varying panel); otherwise `"cs"`.

  `"cs"`

  :   Treats the two sides as independent.

  `"pc"`

  :   Includes same-side cross-period covariance only.

  `"pv"`

  :   Full panel with time-varying running variable: includes same-side
      minus opposite-side cross-period covariance.

- bc:

  Use robust bias-corrected jumps and variances in the LL-Wald
  (Calonico, Cattaneo and Titiunik 2014). `TRUE` (default) aligns the
  test with the bias-corrected
  [`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)
  estimator; `FALSE` uses the conventional local-linear jumps and
  variances.

- ...:

  Currently unused.

## Value

An object of class `"rd_compstable"`, a named list with:

- `pairs`:

  A list, one element per \\(t\_{\mathrm{RD}}, t_0)\\ pair (named
  `"trd::t0"`), each containing:

  `ll_wald`

  :   list with `stat`, `df`, `p`.

  `type_values`

  :   Character vector of \\(\mathbf{u},b)\\ type labels present in this
      pair.

  `scheme`

  :   Scheme actually used.

  `n_trd`

  :   Number of above-cutoff units from \\t_RD\\.

  `n_t0`

  :   Number of above-cutoff units from \\t_0\\.

  `n_both`

  :   Number of units above the cutoff in both periods.

- `joint`:

  Joint result over all pairs (stacked Wald):

  `ll_wald`

  :   list with `stat`, `df`, `p` (sum of the per-pair statistics and
      df, which assumes independent pairs; the pairs share the RD-period
      above-cutoff group, so treat the joint as approximate — the
      paper's test is per pair).

- `meta`:

  list with `t_rd`, `comparisons`, `h` (NA when `bwselect = "cct"`),
  `bwselect`, `c`, `bc`.

## Details

### Reflection construction

For each RD-period / comparison-period pair \\(t\_{\mathrm{RD}}, t_0)\\,
take the above-cutoff units of each period. For \\t_0\\-above units flip
the centred running variable: \\x' = -(R\_{i,t_0} - c)\\, placing them
just BELOW an artificial cutoff at 0. For \\t\_{\mathrm{RD}}\\-above
units set \\x' = R\_{i,t\_{\mathrm{RD}}} - c\\, keeping them just ABOVE
0. Stack the two groups into one artificial cross-section. At the
artificial cutoff the left/right limits of the \\(\mathbf{u},b)\\ type
share are then \\\pi\_{t_0,(+)}(\mathbf{u},b)\\ and
\\\pi\_{t\_{\mathrm{RD}},(+)}(\mathbf{u},b)\\, so a jump at 0 equals the
composition difference. A local-linear RD of each \\(\mathbf{u},b)\\
indicator on the reflected running variable, with a joint Wald that all
jumps are zero (dropping one reference type, since the type shares sum
to 1 and the full set of jumps is rank-deficient), is then this test.

### Unit-level wrinkle

A unit that is above the cutoff in BOTH periods \\t\_{\mathrm{RD}}\\ and
\\t_0\\ appears on BOTH sides of the artificial cutoff (as a
\\t\_{\mathrm{RD}}\\-above observation above the artificial 0 and a
\\t_0\\-above observation below it). The covariance matrix between the
left-side and right-side intercept estimates must therefore include the
id-matched cross-side covariance term (the two g-vectors can share unit
ids). The function computes \\(\text{cov}\_{++} + \text{cov}\_{--} -
\text{cov}\_{+-} - \text{cov}\_{-+})\\ — the same formula as the PV
scheme in the main estimator — rather than assuming the two sides are
independent.

## References

Leventer, D. and Nevo, D. "Correcting Invalid Regression Discontinuity
Designs." Working paper.

## See also

[`rd_typecont()`](https://dorleventer.github.io/rddid/reference/rd_typecont.md),
[`rd_period()`](https://dorleventer.github.io/rddid/reference/rd_period.md),
[`rddid()`](https://dorleventer.github.io/rddid/reference/rddid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-period panel with no composition shift (null DGP).
set.seed(1)
n <- 500
eta <- rnorm(n)
dat <- data.frame(
  id   = rep(seq_len(n), 2),
  time = rep(1:2, each = n),
  R    = c(eta + rnorm(n), eta + rnorm(n))
)
rd_compstable(dat, x = "R", time = "time", id = "id", t_rd = 2,
              comparisons = 1, h = 0.5)
} # }
```
