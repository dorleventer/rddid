# plots.R -- ggplot2 pictures of a fit and of the four assumption tests, mirroring the
# validation figures of the paper's application: binned shares of a type by running variable
# with the local-linear fits on each side (type continuity, composition stability), point-range
# plots of the within-type confounding jumps (homogeneous confounding, constant within-type
# jump), the per-period RD plots behind an rddid() fit, and the switchers picture.
# ggplot2 is in Suggests: every function stops with a clear message when it is not installed.
# Nothing here re-estimates anything: the methods read the fits and data stored in the objects.

utils::globalVariables(c("x", "y", "period", "type", "side", "jump", "se", "lo", "hi",
                         "x_a", "x_b", "status", "role", "xval", "grp", "panel"))

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("the plot methods need the ggplot2 package: install.packages(\"ggplot2\")",
         call. = FALSE)
}

#' Binned means of y by x: `bins` equal-width bins over `range`, one point per non-empty bin
#' @keywords internal
#' @noRd
.bin_means <- function(x, y, bins = 20L, range = NULL) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (is.null(range)) range <- base::range(x)
  if (length(x) == 0L || diff(range) <= 0) return(data.frame(x = numeric(0), y = numeric(0), n = integer(0)))
  edges <- seq(range[1], range[2], length.out = bins + 1L)
  bin   <- pmin(pmax(findInterval(x, edges, rightmost.closed = TRUE), 1L), bins)
  mids  <- (edges[-1L] + edges[-length(edges)]) / 2
  n_bin <- tabulate(bin, nbins = bins)
  keep  <- n_bin > 0
  data.frame(x = mids[keep],
             y = as.numeric(tapply(y, factor(bin, levels = seq_len(bins)), mean))[keep],
             n = n_bin[keep])
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

.rddid_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 11))
}
.cutoff_line <- function(cutoff) {
  ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey40")
}
.pair_colours <- c("Comparison period" = "#1b9e77", "RD period" = "#c2558b")
.type_colours <- c("#1f78b4", "#e6a100", "#6a3d9a", "#33a02c", "#e31a1c", "#b15928")

# ---- rddid fit: the per-period RD plots --------------------------------------------------

#' Plot the per-period RD fits behind an RD-DID estimate
#'
#' One panel per period, in time order: the outcome averaged within `bins` equal-width
#' bins of the running variable, and the two local-linear fits of that period drawn over their
#' bandwidth on each side of the cutoff. The jump between the two lines at the cutoff is the
#' period's discontinuity \eqn{D_t}; the RD-DID estimate is the RD-period jump minus the
#' weighted comparison-period jumps (see `summary()`).
#'
#' @param x an object returned by [rddid()].
#' @param bins number of equal-width bins for the binned means (default 20, over the
#'   running variable's range in each period).
#' @param ... unused.
#' @return A ggplot object (ggplot2 must be installed).
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
    b <- .bin_means(x$data[[k]]$x, x$data[[k]]$y, bins = bins)
    if (nrow(b)) cbind(b, period = k, role = role_of(k)) else NULL
  }))
  lines <- do.call(rbind, lapply(periods, function(k) {
    l <- .fit_lines(x$fits[[k]], cutoff)
    if (!is.null(l)) cbind(l, period = k, role = role_of(k)) else NULL
  }))
  lab <- stats::setNames(paste(vapply(periods, role_of, character(1)), periods), periods)
  pts$period   <- factor(lab[pts$period],   levels = lab[periods])
  lines$period <- factor(lab[lines$period], levels = lab[periods])
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role), alpha = 0.8) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = 1) +
    .cutoff_line(cutoff) +
    ggplot2::facet_wrap(~ period, nrow = 1) +
    ggplot2::scale_colour_manual(values = .pair_colours, guide = "none") +
    ggplot2::labs(x = "Running variable", y = "Outcome") +
    .rddid_theme()
}

# ---- shared: a local-linear fit of a 0/1 indicator at the test's bandwidth rule -----------

#' Local-linear fit of an indicator on x at the cutoff, with the bandwidth rule of a test
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

# ---- rd_typecont: one comparison period against one RD period ------------------------------

#' Plot a type-continuity test: one comparison period against one RD period
#'
#' The paper's type-continuity figure. Two panels: in the comparison period, the share of
#' units that are above the cutoff in the RD period, against the comparison-period running
#' variable; in the RD period, the share that are above the cutoff in the comparison period,
#' against the RD-period running variable. In each panel the share is averaged within `bins`
#' equal-width bins inside the bandwidth and the local-linear fit is drawn on each side of the
#' cutoff, at the bandwidth rule the test used. Under the null the two lines of a panel meet at
#' the cutoff: who a unit is in the other period does not jump there.
#'
#' With more than two periods the test itself uses the full sign pattern of the other periods
#' as the type; the picture shows the pairwise version, one pair of periods at a time.
#'
#' @param x an object returned by [rd_typecont()].
#' @param t_rd,comparison the RD period and the comparison period to draw (values of the
#'   period variable). Defaults: the `t_rd` given to [rd_typecont()] (else the last period),
#'   and the first other period.
#' @param bins number of equal-width bins for the binned shares (default 20).
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
    b <- .bin_means(R, y, bins = bins, range = c(cutoff - h, cutoff + h))
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
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role), alpha = 0.85) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = 1.1) +
    .cutoff_line(cutoff) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = .pair_colours, guide = "none") +
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
#' period of the pair. The binned share is drawn with the local-linear fit on each side, at the
#' bandwidth rule the test used. Under the null the two lines meet at the artificial cutoff:
#' the units just above the cutoff are the same mix in both periods. With
#' `estimand = "atu"` the same picture is drawn for the units below the cutoff.
#'
#' With more than two periods the test's types also record the unit's sides in the remaining
#' periods; the picture shows the pairwise version.
#'
#' @param x an object returned by [rd_compstable()].
#' @param pair which pair to draw: an index or a name of `x$pairs` (default the first).
#' @param bins number of equal-width bins for the binned shares (default 20, over the
#'   bandwidth on each side).
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
  pr <- x$pairs[[pair]]
  if (is.null(pr)) stop("`pair` must be an index or a name of x$pairs: ",
                        paste(names(x$pairs), collapse = ", "))
  pair_name <- if (is.character(pair)) pair else names(x$pairs)[pair]
  t_rd <- sub("::.*$", "", pair_name); t0 <- sub("^.*::", "", pair_name)
  s <- pr$sample
  # the partner-period side is the last character of the type string ("1" = above)
  partner <- function(type) as.numeric(substr(type, nchar(type), nchar(type)) == "1")
  xx <- c(s$x_t0, s$x_trd)
  yy <- c(partner(s$type_t0), partner(s$type_trd))
  fit <- .indicator_fit(yy, xx, 0, x$meta$kernel, x$meta$h, x$meta$bwselect)
  if (is.null(fit)) stop("the local-linear fit on the reflected sample failed.")
  h <- unname(fit$h)
  roles <- c(sprintf("Comparison period %s (mirrored)", t0), sprintf("RD period %s", t_rd))
  pts <- .bin_means(xx, yy, bins = bins, range = c(-h, h))
  pts$role <- factor(ifelse(pts$x < 0, roles[1], roles[2]), levels = roles)
  lines <- .fit_lines(fit, 0)
  lines$role <- factor(ifelse(lines$side == "-", roles[1], roles[2]), levels = roles)
  side_word <- if (identical(x$estimand, "atu")) "below" else "above"
  ggplot2::ggplot() +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = role), alpha = 0.85) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = role, group = side),
                       linewidth = 1.1) +
    .cutoff_line(0) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = stats::setNames(unname(.pair_colours), roles)) +
    ggplot2::labs(x = "Distance to cutoff (comparison period mirrored)",
                  y = sprintf("Pr(%s cutoff in other period)", side_word)) +
    .rddid_theme()
}

# ---- rd_homog / rd_trendcell: the within-type confounding jumps ---------------------------

.jump_pointrange <- function(tab, xvar, colour_var, facet_var = NULL, xlab) {
  tab$lo <- tab$jump - stats::qnorm(0.975) * tab$se
  tab$hi <- tab$jump + stats::qnorm(0.975) * tab$se
  tab$xval <- factor(tab[[xvar]], levels = unique(tab[[xvar]]))
  lv <- unique(tab[[colour_var]])
  lv <- c(lv[startsWith(lv, "Below")], lv[startsWith(lv, "Above")],
          lv[!startsWith(lv, "Below") & !startsWith(lv, "Above")])
  tab$grp  <- factor(tab[[colour_var]], levels = lv)
  cols <- stats::setNames(.type_colours[seq_along(levels(tab$grp))], levels(tab$grp))
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = xval, y = jump, colour = grp, group = grp)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15,
                           position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35), linewidth = 0.6) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35), size = 2.4) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = xlab, y = "Jump at cutoff (95% CI)") +
    .rddid_theme()
  if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(~ grp) +
    ggplot2::theme(legend.position = "none")
  p
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

#' Plot a homogeneous-confounding test: the confounding jump of each type, by comparison period
#'
#' The paper's homogeneous-confounding figure: for each comparison period, the local-linear jump
#' in the outcome at the cutoff within each type (point) with its 95% interval, types side by
#' side. Under the null the types' jumps coincide within each period.
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
  tab$type <- .type_names(tab$type, x$call$t_rd)
  .jump_pointrange(tab, xvar = "period", colour_var = "type", xlab = "Comparison period")
}

#' Plot a constant-within-type-confounding test: each type's confounding jump over time
#'
#' One panel per type: its local-linear jump in each comparison period with its 95% interval,
#' and a dashed line at the type's average jump. Under the null (`trend = "constant"`) the
#' jumps of a type are the same in every comparison period; under `trend = "linear"` they lie
#' on a line.
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
  tab$type <- .type_names(tab$cell, x$call$t_rd)
  p <- .jump_pointrange(tab, xvar = "period", colour_var = "type", facet_var = "type",
                        xlab = "Comparison period")
  means <- stats::aggregate(jump ~ type, data = tab, FUN = mean)
  means$grp <- factor(means$type, levels = levels(p$data$grp))
  p + ggplot2::geom_hline(data = means, ggplot2::aes(yintercept = jump, colour = grp),
                          linetype = "dashed", show.legend = FALSE)
}

# ---- switchers ------------------------------------------------------------------------------

#' Plot the switchers: the running variable in one period against another
#'
#' Each unit observed in both periods is a point; the dashed lines are the cutoff. Units in the
#' off-diagonal quadrants changed side of the cutoff between the two periods (the "switchers"
#' that the composition tests are about); the title gives their share.
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
  above_a <- m$x_a >= cutoff; above_b <- m$x_b >= cutoff
  kind <- ifelse(above_a == above_b, "stay", ifelse(above_a, "down", "up"))
  share <- function(k) 100 * mean(kind == k)
  lab <- c(stay = sprintf("Stayer (%.0f%%)", share("stay")),
           down = sprintf("Above, then below (%.1f%%)", share("down")),
           up   = sprintf("Below, then above (%.1f%%)", share("up")))
  m$status <- factor(lab[kind], levels = lab)
  cols <- stats::setNames(c("grey75", "#d95f02", "#7570b3"), lab)
  ggplot2::ggplot(m, ggplot2::aes(x = x_a, y = x_b, colour = status)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.4) +
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = sprintf("Running variable, period %s", periods[1]),
                  y = sprintf("Running variable, period %s", periods[2])) +
    .rddid_theme()
}
