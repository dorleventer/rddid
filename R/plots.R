# plots.R -- ggplot2 pictures of a fit and of the four assumption tests, in the style of the
# validation figures of the paper's application: binned shares of a type by running variable
# with the local-linear fits on each side (type continuity, composition stability), point-range
# plots of the within-type confounding jumps (homogeneous confounding, constant within-type
# jump), the per-period RD plots behind an rddid() fit, and the switchers picture.
# House rules: no titles, subtitles or statistics in a figure (they are in print()/summary());
# short axis labels; a small bottom legend without a title; theme_bw; one comparison period
# against one RD period, one line per side; periods in time order; stable colours (green =
# comparison period, pink = RD period; blue = below the cutoff in the RD period, orange = above).
# ggplot2 is in Suggests: every function stops with a clear message when it is not installed.
# Nothing here re-estimates the objects' numbers; the share plots fit the pairwise indicator
# they draw, at the bandwidth rule the test used.

utils::globalVariables(c("x", "y", "side", "jump", "lo", "hi", "x_a", "x_b", "status", "role",
                         "xval", "grp", "panel", "period"))

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("the plot methods need the ggplot2 package: install.packages(\"ggplot2\")",
         call. = FALSE)
}

# ---- palette and theme (the paper's) --------------------------------------------------------
.col_pair  <- c("Comparison period" = "#009E73", "RD period" = "#CC79A7")
.col_below <- "#0072B2"
.col_above <- "#E69F00"
.col_more  <- c("#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999")   # 3rd+ types (Okabe-Ito)
.pt_alpha  <- 0.5
.pt_size   <- 1.8
.lwd       <- 1

.rddid_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 9),
                   axis.text = ggplot2::element_text(size = 10),
                   axis.title = ggplot2::element_text(size = 11),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0.8, "lines"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10))
}
.cutoff_line <- function(cutoff) {
  ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey40")
}

# ---- helpers ---------------------------------------------------------------------------------

#' Period labels in time order (numeric order when every label is a number, else as given)
#' @keywords internal
#' @noRd
.period_order <- function(labels) {
  num <- suppressWarnings(as.numeric(labels))
  if (!anyNA(num)) labels[order(num)] else labels
}

#' Binned means of y by x on each side of the cutoff: `bins` equal-width bins from lo to the
#' cutoff and `bins` from the cutoff to hi, so no bin straddles the cutoff; one point per bin
#' @keywords internal
#' @noRd
.bin_sides <- function(x, y, cutoff, bins = 20L, lo = NULL, hi = NULL) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (is.null(lo)) lo <- min(x)
  if (is.null(hi)) hi <- max(x)
  one <- function(xs, ys, from, to) {
    if (!length(xs) || to <= from) return(NULL)
    edges <- seq(from, to, length.out = bins + 1L)
    bin   <- pmin(pmax(findInterval(xs, edges, rightmost.closed = TRUE), 1L), bins)
    mids  <- (edges[-1L] + edges[-length(edges)]) / 2
    n_bin <- tabulate(bin, nbins = bins)
    keep  <- n_bin > 0
    data.frame(x = mids[keep],
               y = as.numeric(tapply(ys, factor(bin, levels = seq_len(bins)), mean))[keep],
               n = n_bin[keep])
  }
  left  <- x < cutoff & x >= lo
  right <- x >= cutoff & x <= hi
  rbind(one(x[left], y[left], lo, cutoff), one(x[right], y[right], cutoff, hi))
}

#' The two fitted lines of an rd_period fit, each over its own side's bandwidth window
#' @keywords internal
#' @noRd
.fit_lines <- function(fit, cutoff, n_points = 25L) {
  if (is.null(fit)) return(NULL)
  h <- unname(fit$h)
  one <- function(side, sign) {
    s  <- fit$sides[[side]]
    xs <- cutoff + sign * seq(0, h, length.out = n_points)
    data.frame(x = xs, y = s$beta0 + s$slope * (xs - cutoff), side = side)
  }
  rbind(one("-", -1), one("+", +1))
}

#' Local-linear fit of a 0/1 indicator on x at the cutoff, with the bandwidth rule of a test
#' (`h` if given, else "cct"/"rot" as the test used); NULL if the fit fails
#' @keywords internal
#' @noRd
.indicator_fit <- function(y, x, cutoff, kernel, h, bwselect) {
  ok <- is.finite(x) & is.finite(y)
  y <- y[ok]; x <- x[ok]
  if (is.na(h)) h <- NULL
  bw <- .cell_bandwidth(y, x, cutoff, kernel, h, bwselect)
  tryCatch(rd_period(y = y, x = x, h = bw[["h"]], b = bw[["b"]], c = cutoff, kernel = kernel),
           error = function(e) NULL)
}

#' Legend-friendly type names: "+" / "-" (side in the RD period) become "Above in t" / "Below in t"
#' @keywords internal
#' @noRd
.type_names <- function(type, t_rd) {
  if (is.null(t_rd)) return(type)
  out <- type
  out[type == "+"] <- sprintf("Above in %s", t_rd)
  out[type == "-"] <- sprintf("Below in %s", t_rd)
  out
}

#' Colours for type names: blue for "Below ...", orange for "Above ...", then the rest
#' @keywords internal
#' @noRd
.type_palette <- function(levels) {
  cols <- character(length(levels))
  cols[startsWith(levels, "Below")] <- .col_below
  cols[startsWith(levels, "Above")] <- .col_above
  rest <- which(cols == "")
  cols[rest] <- if (length(rest) <= length(.col_more)) .col_more[seq_along(rest)] else
    grDevices::hcl.colors(length(rest), "Dark 3")
  stats::setNames(cols, levels)
}

# ---- rddid fit: the per-period RD plots --------------------------------------------------

#' Plot the per-period RD fits behind an RD-DID estimate
#'
#' One panel per period, in time order: the outcome averaged within `bins` equal-width bins
#' of the running variable on each side of the cutoff, and the two local-linear fits of that
#' period drawn over their bandwidth on each side. The jump between the two lines at the cutoff
#' is the period's discontinuity \eqn{D_t}; `summary()` lists every \eqn{D_t} and the RD-DID
#' estimate they combine into. Green panels are comparison periods, pink the RD period.
#'
#' @param x an object returned by [rddid()].
#' @param bins number of equal-width bins on each side of the cutoff (default 20, over the
#'   running variable's range in each period).
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed), which can be changed with `+` and
#'   saved with `ggplot2::ggsave()`.
#' @examples
#' fit <- rddid(rddid_sim, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
#' @family RD-DID estimation
#' @export
plot.rddid <- function(x, bins = 20L, ...) {
  .need_ggplot2()
  cutoff  <- x$c
  periods <- .period_order(names(x$fits))      # panels in time order
  role_of <- function(k) if (k == as.character(x$t_rd)) "RD period" else "Comparison period"
  pts <- do.call(rbind, lapply(periods, function(k) {
    b <- .bin_sides(x$data[[k]]$x, x$data[[k]]$y, cutoff, bins = bins)
    if (!is.null(b) && nrow(b)) cbind(b, period = k, role = role_of(k)) else NULL
  }))
  lines <- do.call(rbind, lapply(periods, function(k) {
    l <- .fit_lines(x$fits[[k]], cutoff)
    if (!is.null(l)) cbind(l, period = k, role = role_of(k)) else NULL
  }))
  lab <- stats::setNames(paste(vapply(periods, role_of, character(1)), periods), periods)
  pts$period   <- factor(lab[pts$period],   levels = lab[periods])
  lines$period <- factor(lab[lines$period], levels = lab[periods])
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role),
                        alpha = .pt_alpha, size = .pt_size) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = .lwd) +
    .cutoff_line(cutoff) +
    ggplot2::facet_wrap(~ period, nrow = 1) +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::scale_colour_manual(values = .col_pair, name = NULL, guide = "none") +
    ggplot2::labs(x = "Running variable", y = "Outcome") +
    .rddid_theme()
}

# ---- rd_typecont: one comparison period against one RD period ------------------------------

#' Plot a type-continuity test: one comparison period against one RD period
#'
#' The paper's type-continuity figure. Two panels: in the comparison period, the share of
#' units that are above the cutoff in the RD period, against the comparison-period running
#' variable; in the RD period, the share that are above the cutoff in the comparison period,
#' against the RD-period running variable. In each panel the share is averaged within `bins`
#' equal-width bins on each side of the cutoff inside the bandwidth, and the local-linear fit
#' is drawn on each side, at the bandwidth rule the test used. Under the null the two lines of
#' a panel meet at the cutoff: who a unit is in the other period does not jump there.
#'
#' With more than two periods the test itself uses the full sign pattern of the other periods
#' as the type; the picture shows the pairwise version, one pair of periods at a time.
#'
#' @param x an object returned by [rd_typecont()].
#' @param t_rd,comparison the RD period and the comparison period to draw (values of the
#'   period variable). Defaults: the `t_rd` given to [rd_typecont()] (else the last period),
#'   and the first other period.
#' @param bins number of equal-width bins on each side of the cutoff (default 20).
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
#' @examples
#' tc <- rd_typecont(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(tc, comparison = 1)
#' @family tests of the assumptions
#' @export
plot.rd_typecont <- function(x, t_rd = NULL, comparison = NULL, bins = 20L, ...) {
  .need_ggplot2()
  periods <- x$meta$periods
  rd <- if (!is.null(t_rd)) as.character(t_rd) else
    if (!is.null(x$meta$t_rd)) as.character(x$meta$t_rd) else periods[length(periods)]
  t0 <- if (!is.null(comparison)) as.character(comparison) else setdiff(periods, rd)[1L]
  if (!rd %in% periods || !t0 %in% periods || rd == t0)
    stop("`t_rd` and `comparison` must be two different periods among: ",
         paste(periods, collapse = ", "))
  cutoff <- x$meta$c
  w <- x$sides
  panel <- function(this, other, role) {
    R <- w[[paste0("R_", this)]]
    y <- as.numeric(w[[paste0("side_", other)]] == "+")
    fit <- .indicator_fit(y, R, cutoff, x$meta$kernel, x$meta$h, x$meta$bwselect)
    if (is.null(fit)) return(NULL)
    h <- unname(fit$h)
    b <- .bin_sides(R, y, cutoff, bins = bins, lo = cutoff - h, hi = cutoff + h)
    lab <- sprintf("%s %s", role, this)
    list(pts = cbind(b, role = role, panel = lab),
         lines = cbind(.fit_lines(fit, cutoff), role = role, panel = lab))
  }
  left  <- panel(t0, rd, "Comparison period")
  right <- panel(rd, t0, "RD period")
  if (is.null(left) || is.null(right)) stop("the local-linear fit failed in one period.")
  pts   <- rbind(left$pts, right$pts)
  lines <- rbind(left$lines, right$lines)
  lv <- c(left$pts$panel[1], right$pts$panel[1])
  pts$panel   <- factor(pts$panel,   levels = lv)
  lines$panel <- factor(lines$panel, levels = lv)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role),
                        alpha = .pt_alpha, size = .pt_size) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = .lwd) +
    .cutoff_line(cutoff) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::scale_colour_manual(values = .col_pair, name = NULL, guide = "none") +
    ggplot2::labs(x = "Running variable", y = "Pr(above cutoff in other period)") +
    .rddid_theme()
}

# ---- rd_compstable: the reflected sample of one pair -------------------------------------

#' Plot a composition-stability test: the reflected sample of one pair of periods
#'
#' The paper's composition-stability figure. The units above the cutoff in the comparison
#' period are placed to the left of an artificial cutoff at their mirrored distance
#' \eqn{-(R_{t_0} - c)}, the units above the cutoff in the RD period to the right at
#' \eqn{R_{t_{RD}} - c}; the outcome is whether the unit is above the cutoff in the other
#' period of the pair. The axis shows the running variable itself (so both sides read outward
#' from the cutoff); the binned share and the local-linear fit on each side are drawn at the
#' bandwidth rule the test used. Under the null the two lines meet at the artificial cutoff:
#' the units just above the cutoff are the same mix in both periods. With
#' `estimand = "atu"` the picture is drawn for the units below the cutoff, with the RD period
#' mirrored instead.
#'
#' With more than two periods the test's types also record the unit's sides in the remaining
#' periods; the picture shows the pairwise version.
#'
#' @param x an object returned by [rd_compstable()].
#' @param pair which pair to draw: an index or a name of `x$pairs` (default the first).
#' @param bins number of equal-width bins on each side (default 20, over the bandwidth).
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
#' @examples
#' cs <- rd_compstable(rddid_sim_pv, x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(cs)
#' @family tests of the assumptions
#' @export
plot.rd_compstable <- function(x, pair = 1L, bins = 20L, ...) {
  .need_ggplot2()
  if (!length(x$pairs)) stop("no pair was tested (every pair was skipped).")
  ok_pair <- (is.character(pair) && pair %in% names(x$pairs)) ||
    (is.numeric(pair) && length(pair) == 1L && pair >= 1 && pair <= length(x$pairs))
  if (!ok_pair) stop("`pair` must be an index 1..", length(x$pairs), " or one of: ",
                     paste(names(x$pairs), collapse = ", "))
  pr <- x$pairs[[pair]]
  pair_name <- if (is.character(pair)) pair else names(x$pairs)[pair]
  t_rd <- sub("::.*$", "", pair_name); t0 <- sub("^.*::", "", pair_name)
  atu <- identical(x$estimand, "atu")
  s <- pr$sample
  # the partner-period side is the last character of the type string ("1" = above)
  partner <- function(type) as.numeric(substr(type, nchar(type), nchar(type)) == "1")
  xx <- c(s$x_t0, s$x_trd)
  yy <- c(partner(s$type_t0), partner(s$type_trd))
  fit <- .indicator_fit(yy, xx, 0, x$meta$kernel, x$meta$h, x$meta$bwselect)
  if (is.null(fit)) stop("the local-linear fit on the reflected sample failed.")
  h <- unname(fit$h)
  # under "att" the comparison period is the mirrored one (left of 0); under "atu" the whole
  # sample was mirrored first, so the RD period (right of 0) is the mirrored one
  roles <- if (atu) c(sprintf("Comparison period %s", t0), sprintf("RD period %s (mirrored)", t_rd))
           else     c(sprintf("Comparison period %s (mirrored)", t0), sprintf("RD period %s", t_rd))
  pts <- .bin_sides(xx, yy, 0, bins = bins, lo = -h, hi = h)
  pts$role <- factor(ifelse(pts$x < 0, roles[1], roles[2]), levels = roles)
  lines <- .fit_lines(fit, 0)
  lines$role <- factor(ifelse(lines$side == "-", roles[1], roles[2]), levels = roles)
  cutoff <- x$meta$c
  # tick labels: the running variable itself on both sides (above the cutoff under "att",
  # below it under "atu"), as in the paper's figure
  fold <- if (atu) function(b) format(cutoff - abs(b), trim = TRUE) else
                   function(b) format(cutoff + abs(b), trim = TRUE)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role),
                        alpha = .pt_alpha, size = .pt_size) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = .lwd) +
    .cutoff_line(0) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(labels = fold, n.breaks = 6) +
    ggplot2::scale_colour_manual(values = stats::setNames(unname(.col_pair), roles), name = NULL) +
    ggplot2::labs(x = "Running variable", y = "Pr(above cutoff in other period)") +
    .rddid_theme()
}

# ---- rd_homog / rd_trendcell: the within-type confounding jumps ---------------------------

.jump_pointrange <- function(tab, colour_var, facet = FALSE, ref_line = NULL) {
  tab$lo   <- tab$jump - stats::qnorm(0.975) * tab$se
  tab$hi   <- tab$jump + stats::qnorm(0.975) * tab$se
  tab$xval <- factor(tab$period, levels = .period_order(unique(tab$period)))
  lv <- unique(tab[[colour_var]])
  lv <- c(lv[startsWith(lv, "Below")], lv[startsWith(lv, "Above")],
          lv[!startsWith(lv, "Below") & !startsWith(lv, "Above")])
  tab$grp <- factor(tab[[colour_var]], levels = lv)
  dodge <- ggplot2::position_dodge(width = 0.35)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = xval, y = jump, colour = grp, group = grp)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.16, linewidth = 0.6,
                           position = dodge) +
    ggplot2::geom_line(position = dodge, linewidth = 0.7) +
    ggplot2::geom_point(position = dodge, size = 2.2) +
    ggplot2::scale_colour_manual(values = .type_palette(lv), name = NULL) +
    ggplot2::labs(x = "Comparison period", y = "Jump at cutoff (95% CI)") +
    .rddid_theme()
  if (!is.null(ref_line)) {
    ref_line$grp <- factor(ref_line$grp, levels = lv)
    ref_line$xval <- factor(ref_line$period, levels = levels(tab$xval))
    p <- p + ggplot2::geom_line(data = ref_line, ggplot2::aes(x = xval, y = jump, colour = grp,
                                                               group = grp),
                                linetype = "dashed", linewidth = 0.6, show.legend = FALSE)
  }
  if (facet) p <- p + ggplot2::facet_wrap(~ grp) + ggplot2::theme(legend.position = "none")
  p
}

#' Plot a homogeneous-confounding test: the confounding jump of each type, by comparison period
#'
#' The paper's homogeneous-confounding figure: for each comparison period, the local-linear jump
#' in the outcome at the cutoff within each type (point) with its 95% interval, types side by
#' side (blue: below the cutoff in the RD period; orange: above). Under the null the types'
#' jumps coincide within each period.
#'
#' @param x an object returned by [rd_homog()].
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
#' @examples
#' hg <- rd_homog(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(hg)
#' @family tests of the assumptions
#' @export
plot.rd_homog <- function(x, ...) {
  .need_ggplot2()
  tab <- x$period_type_jumps
  if (is.null(tab) || !nrow(tab)) stop("no fitted cells to plot.")
  tab$type <- .type_names(tab$type, x$t_rd)
  .jump_pointrange(tab, colour_var = "type")
}

#' Plot a constant-within-type-confounding test: each type's confounding jump over time
#'
#' One panel per type: its local-linear jump in each comparison period with its 95% interval,
#' and a dashed reference line, the type's average jump (`trend = "constant"`) or the
#' least-squares line through its jumps (`trend = "linear"`). Under the null the points sit on
#' the dashed line, up to sampling error.
#'
#' @param x an object returned by [rd_trendcell()].
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
#' @examples
#' tr <- rd_trendcell(rddid_sim_pv, y = "Y", x = "R", time = "year", id = "id", t_rd = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(tr)
#' @family tests of the assumptions
#' @export
plot.rd_trendcell <- function(x, ...) {
  .need_ggplot2()
  tab <- x$cell_period_jumps
  if (is.null(tab) || !nrow(tab)) stop("no fitted cells to plot.")
  tab$type <- .type_names(tab$cell, x$t_rd)
  ref <- do.call(rbind, lapply(split(tab, tab$type), function(d) {
    tv <- suppressWarnings(as.numeric(d$period))
    fitted <- if (identical(x$trend, "linear") && !anyNA(tv) && length(unique(tv)) >= 2L)
      stats::fitted(stats::lm(d$jump ~ tv)) else rep(mean(d$jump), nrow(d))
    data.frame(period = d$period, jump = fitted, grp = d$type[1L])
  }))
  .jump_pointrange(tab, colour_var = "type", facet = TRUE, ref_line = ref)
}

# ---- switchers ------------------------------------------------------------------------------

#' Plot the switchers: the running variable in one period against another
#'
#' Each unit observed in both periods is a point; the dashed lines are the cutoff. Units in the
#' off-diagonal quadrants changed side of the cutoff between the two periods (the "switchers"
#' that the composition tests are about): orange, below in the first period and above in the
#' second; blue, above and then below; grey, the same side in both. A message gives the share
#' of switchers.
#'
#' @inheritParams rddid
#' @param periods the two periods to compare (values of `time`); default the first two in
#'   `data`.
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   plot_switchers(rddid_sim_pv, x = "R", time = "year", id = "id", periods = c(1, 3))
#' @family tests of the assumptions
#' @export
plot_switchers <- function(data, x, time, id, periods = NULL, c = 0, ...) {
  .need_ggplot2()
  .check_columns(data, c(x, time, id))
  cutoff <- c
  all_periods <- sort(unique(data[[time]]))
  if (is.null(periods)) periods <- all_periods[1:2]
  if (length(periods) != 2L || !all(periods %in% all_periods))
    stop("`periods` must be two values of `", time, "`.")
  a <- data[data[[time]] == periods[1], c(id, x)]
  b <- data[data[[time]] == periods[2], c(id, x)]
  names(a) <- c("id", "x_a"); names(b) <- c("id", "x_b")
  m <- merge(a, b, by = "id")
  m <- m[is.finite(m$x_a) & is.finite(m$x_b), ]
  if (!nrow(m)) stop("no unit is observed in both periods (the ids never repeat).")
  above_a <- m$x_a >= cutoff; above_b <- m$x_b >= cutoff
  lab <- c(stay = "Same side", up = "Below, then above", down = "Above, then below")
  m$status <- factor(lab[ifelse(above_a == above_b, "stay", ifelse(above_a, "down", "up"))],
                     levels = lab)
  m <- m[order(m$status != lab[["stay"]]), ]        # switchers drawn on top of the stayers
  message(sprintf("plot_switchers: %.1f%% of the %d units observed in both periods change side",
                  100 * mean(m$status != lab[["stay"]]), nrow(m)))
  cols <- stats::setNames(c("grey75", .col_above, .col_below), lab)   # colour = side in period 2
  ggplot2::ggplot(m, ggplot2::aes(x = x_a, y = x_b, colour = status)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.4) +
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey40") +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::labs(x = sprintf("Running variable, period %s", periods[1]),
                  y = sprintf("Running variable, period %s", periods[2])) +
    .rddid_theme()
}
