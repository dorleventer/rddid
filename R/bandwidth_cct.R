# bandwidth_cct.R -- per-period CCT bandwidths: rd_bw_cct() (exported building block, wraps
# rdrobust::rdbwselect with an IQR fallback on failure) and .bw_cct() (one pair per period).
# Used by rddid(bwselect = "cct"), as the pilot fits of the joint rules, and per cell by the
# assumption tests.

#' Per-period CCT/IK MSE-optimal bandwidths
#'
#' @param plist named list of per-period data frames with columns `y`, `x`, `id`.
#' @keywords internal
#' @noRd
.bw_cct <- function(plist, c = 0, p = 1L, kernel = "triangular") {
  # rd_bw_cct() wraps rdrobust::rdbwselect(bwselect = "mserd") with the
  # finiteness / failure guards and the documented 0.5*IQR fallback.
  lapply(plist, function(d) rd_bw_cct(d$y, d$x, c = c, p = p, kernel = kernel))
}

#' CCT bandwidths for one period (building block)
#'
#' The building block behind `rddid(bwselect = "cct")`, the pilot fits of the
#' `"joint"` and `"iter"` rules, and the default bandwidths of the four tests
#' of the assumptions: the MSE-optimal main and pilot bandwidths of Calonico,
#' Cattaneo and Titiunik (2014) for a local-linear RD in one period, from
#' [rdrobust::rdbwselect()] with `bwselect = "mserd"`. Like [rd_period()] it
#' takes vectors, not a data frame and column names.
#'
#' If `rdbwselect()` fails, or returns a main bandwidth that is not a positive
#' number, `rd_bw_cct()` falls back to `h = b = 0.5 * IQR(x)` (`sd(x)` if the
#' interquartile range is zero) and says so in a message.
#'
#' @param y the outcome (a numeric vector).
#' @param x the running variable (a numeric vector, same length as `y`).
#' @param p order of the local polynomial (default 1, local linear).
#' @inheritParams rddid
#'
#' @return A named numeric vector `c(h = , b = )`: the main bandwidth (point
#'   estimate) and the pilot bandwidth (bias correction).
#'
#' @references
#' Calonico, S., M. D. Cattaneo and R. Titiunik (2014). Robust nonparametric
#' confidence intervals for regression-discontinuity designs. *Econometrica*
#' 82(6), 2295-2326.
#'
#' @family RD-DID estimation
#'
#' @examples
#' # each year of rddid_sim gets its own bandwidths (what bwselect = "cct" uses)
#' sapply(split(rddid_sim, rddid_sim$year), function(d) rd_bw_cct(y = d$Y, x = d$R))
#' @export
rd_bw_cct <- function(y, x, c = 0, p = 1L, kernel = "triangular") {
  fallback <- function(reason) {
    h0 <- 0.5 * stats::IQR(x)
    if (!is.finite(h0) || h0 <= 0) h0 <- stats::sd(x)
    message("rd_bw_cct: CCT unavailable (", reason,
            "); falling back to 0.5*IQR h=", round(h0, 4))
    c(h = h0, b = h0)
  }
  bw <- tryCatch(
    rdrobust::rdbwselect(y = y, x = x, c = c, p = p, kernel = kernel,
                         bwselect = "mserd"),
    error = function(e) e
  )
  if (inherits(bw, "error"))
    return(fallback(conditionMessage(bw)))
  h_val <- as.numeric(bw$bws[1, 1])
  b_val <- as.numeric(bw$bws[1, 3])
  if (!is.finite(h_val) || h_val <= 0)
    return(fallback(paste0("rdbwselect returned h=", h_val)))
  c(h = h_val, b = b_val)
}
