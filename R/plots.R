# ============================================================
# plot_profiles
# ============================================================

#' Plot raw individual trajectories
#'
#' Produces a spaghetti plot of the raw observed trajectories for all
#' subjects (or a selected subset), with an optional group-level mean
#' profile and subject highlighting.
#'
#' @param data A \code{data.frame} containing the data.
#' @param group Character. Name of the grouping variable.
#' @param time Character. Name of the time variable.
#' @param response Character. Name of the response variable.
#' @param subjects Optional character or integer vector of subject IDs to
#'   plot. If \code{NULL} (default) all subjects are plotted.
#' @param mean_profile Logical. Whether to superimpose the group-level mean
#'   trajectory as a bold line. Default \code{TRUE}.
#' @param highlight Optional character or integer vector of subject IDs to
#'   draw in a contrasting colour (e.g., subjects excluded by the
#'   condition-number filter).
#' @param title Character. Plot title. Default \code{"Individual profiles"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' d <- as.data.frame(Theoph)
#' plot_profiles(d, group = "Subject", time = "Time", response = "conc")
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_colour_manual labs theme_bw theme
#' @importFrom stats ave aggregate
#' @importFrom rlang .data
#' @export
plot_profiles <- function(data, group, time, response,
                          subjects      = NULL,
                          mean_profile  = TRUE,
                          highlight     = NULL,
                          title         = "Individual profiles") {



  grp_orig <- tail(make.unique(c(names(data), "grp_orig")), 1)
  data[[grp_orig]] <- data[[group]]
  data[[group]] <- as.numeric(as.factor(data[[group]]))
  if (!is.null(subjects)) data <- data[data[[group]] %in% subjects, , drop = FALSE]
  data$`.group` <- as.character(data[[group]])
  data$`.hl`    <- ifelse(data$`.group` %in% as.character(highlight), "Highlighted", "Normal")
  cols <- c("Normal" = "grey60", "Highlighted" = "#E74C3C")
  data$`.group` <- data[[grp_orig]]

  p <- ggplot2::ggplot(data,
                       ggplot2::aes(x = .data[[time]],
                                    y = .data[[response]],
                                    group = .data$`.group`,
                                    colour = .data$`.hl`)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_colour_manual(values = cols,
                                 name = "",
                                 breaks = c("Highlighted")) +
    ggplot2::labs(title = title,
                  x = time, y = response) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position =
                     if (is.null(highlight)) "none" else "bottom")

  if (mean_profile) {
    data$`.pos` <- ave(seq_len(nrow(data)), data[[group]], FUN = seq_along)
    mn <- aggregate(list(m  = data[[response]],
                         t  = data[[time]]),
                    list(pos = data$`.pos`), mean, na.rm = TRUE)
    mn <- mn[order(mn$pos), ]
    p <- p +
      ggplot2::geom_line(data = mn,
                         ggplot2::aes(x = t, y = m, group = 1),
                         colour = "black", linewidth = 1.2,
                         inherit.aes = FALSE)
  }



  p
}

# ============================================================
# plot_perm_hist
# ============================================================

#' Histogram of the permutation null distribution
#'
#' Regenerates or customises the permutation histogram stored in a
#' \code{"Dtest"} object. The observed statistic \eqn{T_\text{obs}} is
#' shown as a vertical dashed line, and the empirical \eqn{p}-value and
#' rejection decision are annotated on the plot.
#'
#' @param Htest An object of class \code{"Dtest"} returned by
#'   \code{\link{Dhypothesis_test}}.
#' @param bins Integer. Number of histogram bins. Default \code{30}.
#' @param fill Character. Fill colour for the histogram bars.
#'   Default \code{"steelblue"}.
#' @param line_col Character. Colour for the \eqn{T_\text{obs}} line.
#'   Default \code{"red"}.
#' @param title Character. Plot title. If \code{NULL} a default title is
#'   generated.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline coord_cartesian annotate labs theme_bw
#' @importFrom ggplot2 ggplot_build
#' @importFrom stats quantile
#' @export
plot_perm_hist <- function(Htest,
                           bins     = 30,
                           fill     = "steelblue",
                           line_col = "red",
                           title    = NULL) {
  stopifnot(inherits(Htest, "Dtest"))
  lbl <- if (is.null(Htest$bi_out)) "All RE"
  else paste(Htest$bi_out, collapse = ", ")
  if (is.null(title))
    title <- paste0("Permutation null distribution (", lbl, ")")

  Htest$Tperm <- Htest$Tperm[!is.na(Htest$Tperm)]
  df <- data.frame(T = Htest$Tperm)

  P1 <- ggplot2::ggplot(df, ggplot2::aes(x = T)) +
    ggplot2::geom_histogram(bins = bins, fill = fill,
                            colour = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = Htest$Tobs,
                        linetype = "dashed",
                        colour = line_col, linewidth = 1)
  ymax <- max(ggplot2::ggplot_build(P1)$data[[1]]$count)

  ## annotation
  pvalue_num <- attr(Htest, "internal")$pvalue_num
  Pval <- if (pvalue_num < 0.001) {"<0.001"} else {round(pvalue_num, 4)}
  labs <- c(paste0("Tobs = ", round(Htest$Tobs, 4)),
            paste0("p = ", Pval),
            Htest$Decision)
  labs    <- format(labs, justify = "left")
  overlap <- Htest$Tobs > stats::quantile(Htest$Tperm, 0.8)
  xpos    <- if (overlap) min(Htest$Tperm) else  Htest$Tobs

  ## Annotete P1
  P1 +
    ggplot2::coord_cartesian(ylim = c(0, ymax * 1.2)) +
    ggplot2::annotate(
      "text",
      x      = xpos,
      y      = Inf,
      hjust  = if (overlap) 0 else -0.1,
      vjust  = c(1.5, 3, 4.5),
      label  = labs,
      colour = c("red", "black", "black"),
      size   = 4) +
    ggplot2::labs(title = title, x = "Test statistic T", y = "Count") +
    ggplot2::theme_bw()
}

# ============================================================
# plot_fitted
# ============================================================

#' Plot fitted curves overlaid on observed data
#'
#' Draws the model-predicted trajectories on top of the observed data for
#' each subject. When a \code{"Dmethod"} object is supplied, subject-specific
#' EBLUP-adjusted predictions are used; otherwise only the population
#' fixed-effects curve is shown.
#'
#' @param Dobj An object of class \code{"Dmethod"} returned by
#'   \code{\link{Dmethod}}.
#' @param time Character. Name of the time variable in \code{Dobj$data}.
#' @param subjects Optional character or integer vector of subject IDs to
#'   plot. If \code{NULL} (default) all subjects are plotted.
#' @param overlay Logical. If \code{TRUE} (default) observed and fitted
#'   curves are drawn in the same panel per subject. If \code{FALSE},
#'   separate facets are used.
#' @param ncol Integer. Number of columns for the subject facets.
#'   Default \code{4}.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap scale_colour_manual scale_linetype_manual labs theme_bw theme
#' @importFrom MASS ginv
#' @export
plot_fitted <- function(Dobj, time,
                        subjects = NULL,
                        overlay  = TRUE,
                        ncol     = 4) {
  stopifnot(inherits(Dobj, "Dmethod"))

  Dobj     <- c(Dobj, attr(Dobj, "internal")) #internal + public objects
  data     <- Dobj$data
  response <- as.character(Dobj$Expr[[2]])
  group    <- Dobj$group#
  grp_orig <- Dobj$grp_orig
  ids      <- names(Dobj$R)

  if (!is.null(subjects)) ids <- ids[ids %in% as.character(subjects)]

  ## Subject-specific fitted values via EBLUP
  fitted_vals <- lapply(ids, function(id) {
    Ri   <- Dobj$R[[id]]
    ei   <- Dobj$residuals[[id]]
    ni   <- length(ei)
    Dhat <- Dobj$Dhat
    s2   <- Dobj$Sigma2

    V    <- Ri %*% Dhat %*% t(Ri) + s2 * diag(ni)
    Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
    bhat <- Dhat %*% t(Ri) %*% Vinv %*% ei
    as.numeric(Dobj$Ypred[[id]] + Ri %*% bhat)
  })
  names(fitted_vals) <- ids

  ## Assemble plot data
  plot_list <- lapply(ids, function(id) {
    di  <- data[data[[group]] == id, , drop = FALSE]
    ti  <- di[[time]]
    yi  <- di[[response]]
    fi  <- fitted_vals[[id]]
    data.frame(subject  = di[[grp_orig]],
               time     = ti,
               Observed = yi,
               Fitted   = fi)
  })
  pd <- do.call(rbind, plot_list)

  if (overlay) {
    pd_long <- rbind(
      data.frame(subject = pd$subject, time = pd$time,
                 value = pd$Observed, type = "Observed"),
      data.frame(subject = pd$subject, time = pd$time,
                 value = pd$Fitted,   type = "Fitted")
    )
    ggplot2::ggplot(pd_long,
                    ggplot2::aes(x = time, y = value,
                                 colour = type,
                                 linetype = type,
                                 group = interaction(subject, type))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ subject, ncol = ncol) +
      ggplot2::scale_colour_manual(
        values = c(Observed = "grey50", Fitted = "#2980B9"),
        name = "") +
      ggplot2::scale_linetype_manual(
        values = c(Observed = "solid", Fitted = "dashed"),
        name = "") +
      ggplot2::labs(title = "Fitted vs observed trajectories",
                    x = time, y = response) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  } else {
    ggplot2::ggplot(pd) +
      ggplot2::geom_line(ggplot2::aes(x = time, y = Observed),
                         colour = "grey50") +
      ggplot2::geom_line(ggplot2::aes(x = time, y = Fitted),
                         colour = "#2980B9", linetype = "dashed") +
      ggplot2::facet_wrap(~ subject, ncol = ncol) +
      ggplot2::labs(title = "Fitted vs observed trajectories",
                    x = time, y = response) +
      ggplot2::theme_bw()
  }
}

# ============================================================
# plot_residuals
# ============================================================

#' Plot residuals versus fitted values
#'
#' Produces a residuals-versus-fitted-values scatter plot with a horizontal
#' reference line at zero and a loess smoother to aid detection of
#' systematic misfit or heteroscedasticity.
#'
#' @param Dobj An object of class \code{"Dmethod"}.
#' @param time Character. Name of the time variable.
#' @param type Character. Type of residuals: \code{"response"} (default)
#'   for raw residuals \eqn{Y_{ij} - \hat Y_{ij}};
#'   \code{"standardised"} for residuals divided by \eqn{\hat\sigma};
#'   \code{"subject"} for subject-level mean residuals.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth labs theme_bw
#' @importFrom MASS ginv
#' @export
plot_residuals <- function(Dobj, time, type = c("response",
                                                "standardised",
                                                "subject")) {
  stopifnot(inherits(Dobj, "Dmethod"))
  Dobj     <- c(Dobj, attr(Dobj, "internal")) #both internal and public objects
  type     <- match.arg(type)
  data     <- Dobj$data
  response <- as.character(Dobj$Expr[[2]])
  group    <- Dobj$grp_orig
  ids      <- names(Dobj$R)
  s        <- sqrt(Dobj$Sigma2)

  res_list <- lapply(ids, function(id) {
    Ri   <- Dobj$R[[id]]
    ei   <- Dobj$residuals[[id]]
    ni   <- length(ei)
    fhat <- Dobj$Ypred[[id]]

    ## EBLUP subject-specific fitted
    Dhat <- Dobj$Dhat
    V    <- Ri %*% Dhat %*% t(Ri) + Dobj$Sigma2 * diag(ni)
    Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
    bhat <- Dhat %*% t(Ri) %*% Vinv %*% ei
    fhat_ss <- as.numeric(fhat + Ri %*% bhat)
    raw_res <- ei - as.numeric(Ri %*% bhat)

    r <- switch(type,
                response      = raw_res,
                standardised  = raw_res / s,
                subject       = rep(mean(raw_res), ni)
    )
    data.frame(subject = id, fitted = fhat_ss, residual = r)
  })
  pd <- do.call(rbind, res_list)

  ylab <- switch(type,
                 response     = "Residuals",
                 standardised = "Standardised residuals",
                 subject      = "Subject-level mean residuals"
  )

  ggplot2::ggplot(pd, ggplot2::aes(x = fitted, y = residual)) +
    ggplot2::geom_point(colour = "steelblue", alpha = 0.5, size = 1.5) +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                         se = FALSE, colour = "#E74C3C",
                         linewidth = 0.8) +
    ggplot2::labs(title = paste0("Residuals vs fitted (", type, ")"),
                  x = "Fitted values", y = ylab) +
    ggplot2::theme_bw()
}

# ============================================================
# plot_condition
# ============================================================

#' Plot per-subject condition numbers
#'
#' Displays the condition numbers \eqn{\kappa_i} from a \code{"MM_base"}
#' object or the \code{MM_base_obj} slot of a \code{"Dmethod"} object,
#' sorted in decreasing order. A horizontal reference line at
#' \code{kappa_max} marks the exclusion threshold.
#'
#' @param Dobj A \code{"MM_base_obj"} object from \code{\link{MM_base}}, or a
#'   \code{"Dmethod"} object (for MM/MMF methods only).
#' @param kappa_max Numeric. Exclusion threshold to highlight on the plot.
#'   Default \code{1e4}.
#' @param log_scale Logical. If \code{TRUE} (default), the \eqn{y}-axis
#'   is on the \eqn{\log_{10}} scale.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_colour_manual labs theme_bw theme element_text scale_y_log10
#' @export
plot_condition <- function(Dobj,  kappa_max = 1e4, log_scale = TRUE) {

  stopifnot(inherits(Dobj, "Dmethod"))
  stopifnot(Dobj$method %in% c("MM", "MMF") )

  Dobj      <- c(Dobj, attr(Dobj, "internal")) #both internal and public objects
  obj       <- Dobj$MM_base_obj
  group     <- Dobj$group
  grp_orig  <- Dobj$grp_orig
  data      <- Dobj$data[!duplicated(Dobj$data[[group]]), ]
  data$kapa <- obj$kappa[as.character(data[[group]])]
  data      <- data[ !is.na(data$kapa), ]
  data      <- data[order(data$kapa, decreasing = TRUE), ]


  pd <- data.frame(
    subject = factor(data[[grp_orig]], levels = data[[grp_orig]]),
    kappa_val     = data$kapa,
    status  = ifelse(data$kapa > kappa_max, "Excluded", "Retained")
  )

  p <- ggplot2::ggplot(pd,
                       ggplot2::aes(x = subject, y = kappa_val  ,
                                    colour = status, shape = status)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = kappa_max,
                        linetype = "dashed", colour = "red") +
    ggplot2::scale_colour_manual(
      values = c(Retained = "steelblue", Excluded = "#E74C3C")) +
    ggplot2::labs(title = "Per-subject condition numbers",
                  x = "Subject (sorted by kappa)",
                  y = expression(kappa[i]),
                  colour = "", shape = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       size  = 7),
                   legend.position = "bottom")

  if (log_scale) p <- p + ggplot2::scale_y_log10()
  p
}

# ============================================================
# plot_variance
# ============================================================

#' Compare variance component estimates across methods
#'
#' Produces a side-by-side dot plot (or bar chart) of the estimated
#' variance components \eqn{\hat d_{*jj}} for each random effect,
#' comparing results from multiple estimation methods.
#'
#' @param Dhatt_list A named list of \code{"Dmethod"} objects, one per
#'   method. For example,
#'   \code{list(VLS = DVLS, MM = DMM, MMF = DMMF)}.
#' @param component Character. \code{"diagonal"} (default) plots only
#'   the variance terms \eqn{\hat d_{*jj}}; \code{"full"} additionally
#'   includes the covariance terms.
#' @param title Character. Plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes scale_x_discrete geom_point position_dodge geom_hline labs theme_bw theme
#' @export
plot_variance <- function(Dhatt_list,
                          component = c("diagonal", "full"),
                          title     = "Variance component estimates") {
  component <- match.arg(component)
  stopifnot(is.list(Dhatt_list),
            all(sapply(Dhatt_list, inherits, "Dmethod")))

  rows <- lapply(names(Dhatt_list), function(nm) {
    Dhat <- Dhatt_list[[nm]]$Dhat
    k    <- nrow(Dhat)
    if (component == "diagonal") {
      data.frame(
        method    = nm,
        parameter = paste0("d['*", seq_len(k), seq_len(k), "']"),
        estimate  = diag(Dhat)
      )
    } else {
      idx <- which(lower.tri(Dhat, diag = TRUE), arr.ind = TRUE)
      data.frame(
        method    = nm,
        parameter = paste0("d['*", idx[,1], idx[,2], "']"),
        estimate  = Dhat[idx]
      )
    }
  })
  pd <- do.call(rbind, rows)
  pd$method <- factor(pd$method, levels = names(Dhatt_list))

  ggplot2::ggplot(pd,
                  ggplot2::aes(x = parameter, y = estimate,
                               colour = method, shape = method,
                               group  = method)) +
    ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    ggplot2::labs(title  = title,
                  x      = "Variance component",
                  y      = "Estimate",
                  colour = "Method",
                  shape  = "Method") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

# ============================================================
# S3 plot dispatch for Dmethod and Dtest
# ============================================================

#' Plot method for \code{Dmethod} objects
#'
#' By default all relevant diagnostic plots are produced and arranged in a
#' grid. Use the \code{which} argument to select specific plots.
#'
#' @param x A \code{"Dmethod"} object.
#' @param time Character. Name of the time variable (required).
#' @param which Integer vector specifying which plots to show:
#'   \code{1} = individual profiles,
#'   \code{2} = fitted vs observed,
#'   \code{3} = residuals vs fitted (response),
#'   \code{4} = standardised residuals,
#'   \code{5} = condition numbers (MM/MMF only).
#'   Default: all applicable plots.
#' @param ... Additional arguments passed to individual plot functions.
#'
#' @return Invisibly returns a list of \code{ggplot2} objects.
#'
#' @export
plot.Dmethod <- function(x, time, which = NULL, ...) {
  stopifnot(inherits(x, "Dmethod"))
  Dobj <- c(x, attr(x, "internal"))   ## merged public+internal objects
  if (missing(time))
    stop("plot.Dmethod: please supply the 'time' argument (name of the time variable).",
         call. = FALSE)

  plots <- list()
  avail <- c(1, 2, 3, 4)
  if (Dobj$method %in% c("MM", "MMF") && !is.null(Dobj$MM_base_obj))
    avail <- c(avail, 5)
  if (is.null(which)) which <- avail

  if (1 %in% which)
    plots[["profiles"]] <- plot_profiles(Dobj$data, Dobj$group, time,
                                         as.character(Dobj$Expr[[2]]), ...)
  if (2 %in% which)
    plots[["fitted"]] <- plot_fitted(x, time, ...)
  if (3 %in% which)
    plots[["residuals_raw"]] <- plot_residuals(x, time,
                                               type = "response", ...)
  if (4 %in% which)
    plots[["residuals_std"]] <- plot_residuals(x, time,
                                               type = "standardised", ...)
  if (5 %in% which && Dobj$method %in% c("MM", "MMF"))
    plots[["condition"]] <- plot_condition(x, ...)

  ## Display arranged grid
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    print(ggpubr::ggarrange(plotlist = plots,
                            ncol = 2,
                            nrow = ceiling(length(plots) / 2)))
  } else {
    for (p in plots) print(p)
  }
  invisible(plots)
}

#' Plot method for \code{Dtest} objects
#'
#' Displays the permutation null distribution histogram with the observed
#' test statistic annotated.
#'
#' @param x A \code{"Dtest"} object.
#' @param ... Additional arguments passed to \code{\link{plot_perm_hist}}.
#'
#' @return Invisibly returns the \code{ggplot2} object.
#'
#' @export
plot.Dtest <- function(x, ...) {
  p <- plot_perm_hist(x, ...)
  print(p)
  invisible(p)
}
