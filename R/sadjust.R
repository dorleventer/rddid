# Within-period composition term S_t of Proposition `sdisc` (Section S3.2).
#
# The observed period-t discontinuity decomposes (Prop sdisc) as
#   D_t = att(t) * 1{t in T_RD} + alpha_{t,0} + S_t,
# where the within-period composition term is
#   S_t = sum_a ( pi_{t,+}(a) - pi_{t,-}(a) ) * mu_{t,(0,0),-}(a).
# With BINARY types a in {above, below} = the unit's side of the cutoff in the
# PARTNER period, the type shares sum to one, so
#   S_t = dpi_above * ( mu_above^(-) - mu_below^(-) ),
# where dpi_above is the cutoff JUMP in the above-type share (a single RD on the
# type indicator; the Assumption-A7 quantity) and mu_a^(-) is the BELOW-cutoff
# (control-side) outcome limit among type-a units.
#
# KEY PROPERTY: mu_a^(-) uses the control (below-cutoff) side, observed in EVERY
# period including the RD period (treatment only switches on above the cutoff),
# so S_t is estimable in ALL periods -- unlike rd_adjust's cross-period C term.
# The composition-free discontinuity is D_t - S_t.
#
# Pairwise (P=2) design, matching the rest of S5: for each comparison period t0,
# form the pair {t0, t_rd} and report S_{t0} (type = t_rd side) together with
# S_{t_rd} (type = t0 side). Bandwidths default to per-(period, type, fit) CCT;
# in the bootstrap they are held FIXED at the point-estimate values (recomputing
# CCT on the thin switcher cells is unstable).

# ---- internal: each unit's side (1 above, 0 below) of the cutoff in `period` --
.rd_sadjust_side <- function(data, x, time, id, period, c) {
  keep <- data[[time]] == period
  stats::setNames(as.integer(data[[x]][keep] >= c),
                  as.character(data[[id]][keep]))
}

# ---- internal: focal-period frame with a partner-side type indicator ----------
# Reads the type (partner-period side) per unit from the running variable in the
# partner period; units absent from the partner period get type NA and drop out.
.rd_sadjust_frame <- function(data, y, x, time, id, focal, partner, c) {
  side_partner <- .rd_sadjust_side(data, x, time, id, partner, c)
  keep <- data[[time]] == focal
  ids  <- as.character(data[[id]][keep])
  data.frame(
    y          = as.numeric(data[[y]][keep]),
    x          = as.numeric(data[[x]][keep]),
    type_above = unname(side_partner[ids]),
    stringsAsFactors = FALSE
  )
}

# ---- internal: point-estimate CCT bandwidths for a frame's four fits ----------
.rd_sadjust_pointbw <- function(frame, c, kernel) {
  d  <- frame[stats::complete.cases(frame$y, frame$x, frame$type_above), ,
              drop = FALSE]
  ab <- d$type_above == 1L
  be <- d$type_above == 0L
  cct <- function(yy, xx) rd_bw_cct(yy, xx, c = c, kernel = kernel)
  list(pi = cct(d$type_above, d$x), ab = cct(d$y[ab], d$x[ab]),
       be = cct(d$y[be], d$x[be]), D = cct(d$y, d$x))
}

# ---- internal: S_t (conv + bc) and D_t (conv + bc) at FIXED bandwidths --------
# bws is list(pi=, ab=, be=, D=), each a named c(h=, b=).
.rd_sadjust_one <- function(frame, bws, c, kernel) {
  d  <- frame[stats::complete.cases(frame$y, frame$x, frame$type_above), ,
              drop = FALSE]
  ab <- d$type_above == 1L
  be <- d$type_above == 0L
  rp <- function(yy, xx, key)
    rd_period(y = yy, x = xx, h = unname(bws[[key]]["h"]),
              b = unname(bws[[key]]["b"]), c = c, kernel = kernel)
  fpi <- rp(d$type_above, d$x, "pi")
  fab <- rp(d$y[ab], d$x[ab], "ab")
  fbe <- rp(d$y[be], d$x[be], "be")
  fD  <- rp(d$y, d$x, "D")
  dpi    <- fpi$D
  dpi_bc <- fpi$D_bc
  mu_ab    <- fab$sides[["-"]]$beta0
  mu_ab_bc <- fab$sides[["-"]]$beta0_bc
  mu_be    <- fbe$sides[["-"]]$beta0
  mu_be_bc <- fbe$sides[["-"]]$beta0_bc
  c(S    = dpi    * (mu_ab    - mu_be),
    S_bc = dpi_bc * (mu_ab_bc - mu_be_bc),
    D    = fD$D,
    D_bc = fD$D_bc)
}

#' Within-period composition-adjusted RD estimator (Proposition `sdisc`)
#'
#' Estimates the within-period composition term \eqn{S_t} that contaminates the
#' observed period-\eqn{t} discontinuity \eqn{D_t} when units sort across the
#' cutoff between regimes, and reports the composition-free discontinuity
#' \eqn{D_t - S_t}. With binary types (the unit's side of the cutoff in the
#' partner period) the term collapses to
#' \eqn{S_t = \Delta\pi_{above}\,(\mu_{above}^{(-)} - \mu_{below}^{(-)})},
#' where \eqn{\Delta\pi_{above}} is the cutoff jump in the above-type share and
#' \eqn{\mu_a^{(-)}} is the below-cutoff (control-side) outcome limit among
#' type-\eqn{a} units. Because \eqn{\mu_a^{(-)}} uses the control side it is
#' observed in every period, so \eqn{S_t} is estimable even at the RD period.
#'
#' Uses the pairwise (P = 2) design: for each comparison period in
#' `comparisons`, the pair `{t0, t_rd}` yields \eqn{S_{t0}} (type = the unit's
#' `t_rd` side) and \eqn{S_{t\_rd}} (type = the unit's `t0` side). Standard
#' errors are a unit-level cluster bootstrap at bandwidths held fixed at the
#' point-estimate values.
#'
#' @param data long data frame, one row per unit-period.
#' @param y,x,time,id column names (strings) for outcome, running variable,
#'   period, and unit id.
#' @param comparisons comparison-period values of `time`; the pair `{t0, t_rd}`
#'   is formed for each.
#' @param t_rd RD-period value of `time`.
#' @param c cutoff (default 0).
#' @param h optional fixed numeric bandwidth used for every fit; if `NULL`
#'   (default) a per-(period, type, fit) CCT bandwidth is selected.
#' @param bwselect bandwidth selector label when `h` is `NULL` (default
#'   `"cct"`); retained for forward compatibility.
#' @param kernel `"triangular"` (default), `"epanechnikov"`, or `"uniform"`.
#' @param se `"bootstrap"` (default) or `"none"`.
#' @param B bootstrap replications (default 499).
#'
#' @return An object of class `"rd_sadjust"`: a list with a tidy per-period
#'   `table` (columns `comparison`, `period`, `role`, `S`, `S_bc`, `S_se`,
#'   `S_bc_se`, `D`, `D_bc`, `adjusted`, `adjusted_bc`, `n`,
#'   `n_control_switchers`) and meta fields.
#' @export
rd_sadjust <- function(data, y, x, time, id, comparisons, t_rd, c = 0,
                       h = NULL, bwselect = "cct", kernel = "triangular",
                       se = c("bootstrap", "none"), B = 499L) {
  cl <- match.call()
  se <- match.arg(se)
  for (nm in c(y, x, time, id))
    if (!nm %in% names(data)) stop("column '", nm, "' not found in `data`.")
  comparisons <- sort(unique(comparisons))
  if (length(comparisons) < 1L) stop("supply at least one comparison period.")
  if (t_rd %in% comparisons) stop("`t_rd` must not be in `comparisons`.")

  # one (focal, partner, role) spec per comparison pair x role
  specs <- list()
  for (t0 in comparisons) {
    specs[[length(specs) + 1L]] <-
      list(comparison = t0, focal = t0, partner = t_rd, role = "comparison")
    specs[[length(specs) + 1L]] <-
      list(comparison = t0, focal = t_rd, partner = t0, role = "rd")
  }
  K <- length(specs)

  # frames + bandwidths fixed at the point estimate
  frames <- lapply(specs, function(s)
    .rd_sadjust_frame(data, y, x, time, id, s$focal, s$partner, c))
  bws <- lapply(frames, function(fr) {
    if (!is.null(h)) {
      hp <- c(h = h, b = h)
      list(pi = hp, ab = hp, be = hp, D = hp)
    } else {
      .rd_sadjust_pointbw(fr, c, kernel)
    }
  })

  # point estimate: 4 x K matrix (rows S, S_bc, D, D_bc)
  pt <- vapply(seq_len(K), function(k)
    .rd_sadjust_one(frames[[k]], bws[[k]], c, kernel), numeric(4))
  rownames(pt) <- c("S", "S_bc", "D", "D_bc")

  # ---- bootstrap SE (resample units, rebuild frames, recompute at fixed bw) ---
  se_mat <- NULL
  n_boot_ok <- NA_integer_
  if (se == "bootstrap") {
    idx_by_id <- split(seq_len(nrow(data)), as.character(data[[id]]))
    U <- names(idx_by_id)
    arr <- array(NA_real_, dim = c(B, 4L, K))
    for (rb in seq_len(B)) {
      drawn <- sample(U, length(U), replace = TRUE)
      rows_list <- idx_by_id[drawn]
      bd <- data[unlist(rows_list, use.names = FALSE), , drop = FALSE]
      bd[[id]] <- rep(seq_along(drawn), lengths(rows_list))   # fresh cluster ids
      bframes <- lapply(specs, function(s)
        .rd_sadjust_frame(bd, y, x, time, id, s$focal, s$partner, c))
      for (k in seq_len(K))
        arr[rb, , k] <- tryCatch(
          .rd_sadjust_one(bframes[[k]], bws[[k]], c, kernel),
          error = function(e) rep(NA_real_, 4L))
    }
    se_mat <- apply(arr, c(2L, 3L), stats::sd, na.rm = TRUE)   # 4 x K
    rownames(se_mat) <- c("S", "S_bc", "D", "D_bc")
    n_boot_ok <- sum(stats::complete.cases(matrix(arr[, 1L, ], nrow = B)))
  }

  # per-frame counts
  ns <- vapply(frames, function(fr) {
    d <- fr[stats::complete.cases(fr$y, fr$x, fr$type_above), , drop = FALSE]
    c(n = nrow(d), nsw = sum(d$x < c & d$type_above == 1L))
  }, numeric(2))

  S    <- pt["S", ];    S_bc <- pt["S_bc", ]
  D    <- pt["D", ];    D_bc <- pt["D_bc", ]
  S_se    <- if (is.null(se_mat)) rep(NA_real_, K) else se_mat["S", ]
  S_bc_se <- if (is.null(se_mat)) rep(NA_real_, K) else se_mat["S_bc", ]

  tab <- data.frame(
    comparison          = vapply(specs, function(s) as.numeric(s$comparison), numeric(1)),
    period              = vapply(specs, function(s) as.numeric(s$focal), numeric(1)),
    role                = vapply(specs, function(s) s$role, character(1)),
    S                   = S,
    S_bc                = S_bc,
    S_se                = S_se,
    S_bc_se             = S_bc_se,
    D                   = D,
    D_bc                = D_bc,
    adjusted            = D - S,
    adjusted_bc         = D_bc - S_bc,
    n                   = ns["n", ],
    n_control_switchers = ns["nsw", ],
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )

  structure(list(
    table       = tab,
    t_rd        = t_rd,
    comparisons = comparisons,
    c           = c,
    h           = h,
    bwselect    = if (is.null(h)) bwselect else "fixed",
    kernel      = kernel,
    se          = se,
    B           = if (se == "bootstrap") B else 0L,
    n_boot_ok   = n_boot_ok,
    call        = cl
  ), class = "rd_sadjust")
}

#' @export
print.rd_sadjust <- function(x, ...) {
  cat("Within-period composition-adjusted RD estimate (Proposition sdisc)\n")
  cat(sprintf("  RD period: %s   comparison periods: %s   bw: %s\n",
              x$t_rd, paste(x$comparisons, collapse = ", "), x$bwselect))
  tab <- x$table
  sefmt <- function(v, s) if (is.na(s)) sprintf("%+.4g", v) else
    sprintf("%+.4g (%.3g)", v, s)
  for (i in seq_len(nrow(tab)))
    cat(sprintf(
      "  t=%g [%-10s | pair %g]: S=%s  D=%+.4g  D-S=%+.4g  (n=%d, sw=%d)\n",
      tab$period[i], tab$role[i], tab$comparison[i],
      sefmt(tab$S[i], tab$S_se[i]), tab$D[i], tab$adjusted[i],
      tab$n[i], tab$n_control_switchers[i]))
  if (x$se == "bootstrap")
    cat(sprintf("  SE: unit cluster bootstrap, B=%d (%d ok)\n",
                x$B, x$n_boot_ok))
  invisible(x)
}
