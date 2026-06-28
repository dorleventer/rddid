# rddid

<!-- badges: start -->
[![R-CMD-check](https://github.com/dorleventer/rddid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dorleventer/rddid/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Estimation and inference for regression-discontinuity
difference-in-discontinuities (RD-DID) designs, where a confounding policy
switches at the same cutoff as the treatment of interest. Implements the
framework of Leventer and Nevo.

Paper: <https://arxiv.org/abs/2408.05847>

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dorleventer/rddid")
```

## Development

After cloning, enable the doc-sync pre-commit hook (regenerates `man/*.Rd`
from roxygen and blocks commits where the generated docs are stale — the
mismatch that otherwise fails `R CMD check`):

``` sh
git config core.hooksPath .githooks
```
