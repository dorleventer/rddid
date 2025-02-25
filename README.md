
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rddid

<!-- badges: start -->
<!-- badges: end -->

The goal of rddid is to estimate treatment effects in RDDID settings,
implementing the estimation framework proposed in Leventer and Nevo
(2024). Functions for RD estimation in single time periods are based on
[rdrobust](https://github.com/rdpackages/rdrobust) code.

## Installation

You can install the development version of rddid from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("dorleventer/rddid")
```

## In Progress

This package in currently a work in progess. The following items are on
the to do list

- [ ] linear_in_time_period_est
  - [ ] Generalize function inputs
  - [ ] Generalize to more than two time periods (requires analytical
    derivation)

## Example

For examples see

- Estimation under constant PO discontinuity assumption
  [notebook](https://dorleventer.github.io/rddid/articles/rddid-constant.html).
- Estimation under linear in time PO discontinuity assumption
  [notebook](https://dorleventer.github.io/rddid/articles/rddid-linear-trend.html).

## References

Leventer, D., & Nevo, D. (2024). Correcting invalid regression
discontinuity designs with multiple time period data. arXiv preprint
[arXiv:2408.05847](https://arxiv.org/abs/2408.05847).
