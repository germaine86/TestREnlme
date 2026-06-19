#' Compute the observed test statistic
#'
#' Computes \eqn{T_\text{obs} = N^{-1} \sum_i \text{tr}(R_i \hat D_* R_i^T)}
#' from a fitted \code{\link{Dmethod}} object.
#'
#' @param Dobj An object of class \code{"Dmethod"} returned by
#'   \code{\link{Dmethod}}.
#' @param bi_out Optional character vector of random-effect names to test.
#'   If \code{NULL} (default) all random effects are included.
#'
#' @return A scalar, the value of \eqn{T_\text{obs}}.
#'
#' @importFrom stats complete.cases
#' @export
Tstat <- function(Dobj, bi_out = NULL) {
  stopifnot(inherits(Dobj, "Dmethod"))
  stopifnot(is.null(bi_out) || all(bi_out %in% colnames(Dobj$Dhat)))
  R_list <- Dobj$R
  Dhat   <- Dobj$Dhat
  N      <- length(R_list)

  if (!is.null(bi_out)) {
    idx    <- which(colnames(Dhat) %in% bi_out)
    Dhat   <- Dhat[idx, idx, drop = FALSE]
    R_list <- lapply(R_list, function(Ri) Ri[, idx, drop = FALSE])
  }

  sum(vapply(R_list, function(Ri) {
    obs <- complete.cases(Ri)      ## remove NA rows
    Ri  <- Ri[obs, , drop = FALSE]
    sum((Ri %*% Dhat) * Ri)}, numeric(1))) / N
}

# ============================================================
#' Permutation test for random effects in linear and nonlinear mixed-effects models
#'
#' Performs a nonparametric permutation test for all random effects or any
#' user-specified subset, using the test statistic of Drikvandi et al.
#' (2013) adapted for nonlinear mixed-effects models in Uwimpuhwe (2026).
#'
#' @param data A \code{data.frame} containing all model variables.
#' @param Expr A two-sided formula specifying the nonlinear model
#'   \eqn{f_i(a_i, \gamma)}. The left-hand side is the response variable
#'   and the right-hand side defines the nonlinear function using the
#'   subject-specific parameter names (e.g., \code{ai1}, \code{ai2},
#'   \code{ai3}) that appear in \code{start} and \code{random}.
#' @param group Character. Name of the grouping variable.
#' @param random A character vector of two-sided formula strings, one per
#'   parameter in \eqn{a_i = A_i\beta + b_i}, mapping each subject-specific
#'   parameter (left-hand side, matching \code{Expr} and \code{start})
#'   to its fixed-effects expression plus random effect (right-hand side).
#'   For example, \code{c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2",
#'   "ai3 ~ B3 + bi3")} specifies that \code{ai1 = B1 + bi1},
#'   \code{ai2 = B2 + bi2}, \code{ai3 = B3 + bi3}.
#' @param start A named numeric vector of starting values for all parameters
#'   in \code{Expr}. Names must match those used in \code{Expr} (e.g.,
#'   \code{ai1}, \code{ai2}, \code{ai3}). If not supplied, \code{Dmethod}
#'   attempts to compute starting values automatically using
#'   \code{nls.multstart::nls_multstart()}, searching over multiple
#'   initial values within a specified range (e.g. \eqn{\pm 10}) for each
#'   parameter. If this step fails, provide \code{start} manually, or fit
#'   \code{nls_multstart()} or \code{nls()} separately with different
#'   starting values, search bounds, or optimisation settings, and use the
#'   resulting coefficients as starting values.
#' @param bi_out Optional character vector naming the random effects to
#'   test. If \code{NULL} (default) all random effects are tested jointly
#'   (Algorithm 1). Otherwise, only the named effects are tested
#'   conditionally on the rest being retained.
#' @param method Character. One of \code{"VLS"} (default),
#'   \code{"MM"}, or \code{"MMF"}.
#' @param Dhatt Optional pre-computed \code{"Dmethod"} object. When
#'   supplied, variance estimation is skipped. This is the recommended
#'   approach when running multiple tests on the same data, since
#'   \eqn{\hat D_*} only needs to be computed once and can be reused.
#' @param Thatt Optional pre-computed observed test statistic from
#'   \code{\link{Tstat}}.
#' @param MM_base_obj Optional pre-computed \code{"MM_base"} object. Only used
#'   when \code{method} is \code{"MM"} or \code{"MMF"} and \code{Dhatt}
#'   is \code{NULL}.
#' @param nperm Positive integer. Number of permutations \eqn{B}.
#'   Default \code{200}. Use \code{1000} or more for publication results.
#' @param seed Optional integer seed for reproducibility. Default \code{NULL}.
#' @param sig_alpha Significance level for the reject/do-not-reject
#'   decision. Default \code{0.05}.
#' @param kappa_max Condition-number threshold for MM/MMF. Default \code{1e4}.
#' @param RR_catof Exclusion criterion passed to \code{\link{MM_base}}.
#'   Either \code{"kappa"} (default), which uses the condition-number
#'   threshold \code{kappa_max}, or a user-specified numeric threshold
#'   applied directly.
#' @param verbose Integer (0, 1, or 2). Controls message output:
#'   \describe{
#'     \item{\code{0}}{Completely silent.}
#'     \item{\code{1}}{Prints summary: subjects used/excluded, observed
#'       statistic, \eqn{p}-value, and decision (default).}
#'     \item{\code{2}}{Also prints every \code{perm_freq} permutation counter.}
#'   }
#' @param perm_freq Integer. When \code{verbose = 2}, print a permutation
#'   progress message every \code{perm_freq} permutations. Default \code{10}.
#'
#' @return An object of class \code{"Dtest"}, a list with components:
#'   \describe{
#'     \item{\code{Decision}}{Character, \code{"Reject H0"} or
#'       \code{"Do not reject H0"}.}
#'     \item{\code{pvalue}}{The empirical permutation \eqn{p}-value.}
#'     \item{\code{Tobs}}{The observed test statistic \eqn{T_\text{obs}}.}
#'     \item{\code{Tperm}}{Numeric vector of length \code{nperm} containing
#'       the permutation statistics \eqn{T^{(1)}, \ldots, T^{(B)}}.}
#'     \item{\code{Dhatt}}{The \code{Dmethod} object used.}
#'     \item{\code{bi_out}}{The random effects tested.}
#'     \item{\code{plot}}{A \pkg{ggplot2} histogram of the permutation null
#'       distribution with \eqn{T_\text{obs}} annotated.}
#'   }
#'
#' @examples
#' \dontrun{
#' d      <- as.data.frame(Theoph)
#' Expr   <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
#'             (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
#'             (exp(ai2) - exp(ai3))
#' start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
#' random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")
#' DVLS  <- Dmethod(d, Expr, group = "Subject",
#'                    random = random, start = start)
#' H      <- Dhypothesis_test(d, Expr, group = "Subject",
#'                            random = random, start = start,
#'                            Dhatt = DVLS, nperm = 200, seed = 1)
#' H$pvalue
#' H$plot
#' }
#'
#' @references
#' Uwimpuhwe, G., Drikvandi, R., and Blozis, S.A. (2026).
#' Testing random effects in nonlinear mixed-effects models.
#' \emph{Statistics in Medicine}.
#' \doi{10.1002/sim.70605}
#'
#' Uwimpuhwe, G., Drikvandi, R. and Blozis, S. A. (in preparation).
#' TestREnlme: An R Package for Testing Random Effects in Nonlinear
#' Mixed-Effects Models. \emph{Journal of Statistical Software}.
#'
#' Demidenko, E. (2013).
#' \emph{Mixed Models: Theory and Applications with R} (2nd ed.). Wiley.
#'
#' Drikvandi, R., Verbeke, G., Khodadadi, A. and Nia, V. P. (2013).
#' Testing multiple variance components in linear mixed-effects models.
#' \emph{Biostatistics}, \strong{14 (1)}, 144--159.
#'
#' @importFrom stats setNames
#' @importFrom MASS ginv
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate labs theme_bw
#' @export
Dhypothesis_test <- function(data, Expr, group, random, start,
                             bi_out    = NULL,
                             method    = c("VLS", "MM", "MMF"),
                             Dhatt     = NULL,
                             Thatt     = NULL,
                             MM_base_obj= NULL,
                             nperm     = 200,
                             seed      = NULL,
                             sig_alpha = 0.05,
                             kappa_max = 1e4,
                             RR_catof  = "kappa",
                             verbose   = 1,
                             perm_freq = 10) {

  method   <- match.arg(method)

  ## --- Step 1: estimate D* if not supplied ---------------------------
  if (is.null(Dhatt)) {
    .vcat(verbose, 1, "\nDhypothesis_test: computing variance components ...")
    if(missing(start)){
    Dhatt <- Dmethod(data, Expr, group, random,
                     method = method, MM_base_obj = MM_base_obj,
                     kappa_max = kappa_max, RR_catof = RR_catof,
                     verbose = verbose)
    start <-  attr(Dhatt, "internal")$start_orig
    }else{
      Dhatt <- Dmethod(data, Expr, group, random, start,
                       method = method, MM_base_obj = MM_base_obj,
                       kappa_max = kappa_max, RR_catof = RR_catof,
                       verbose = verbose)
    }
  }

  ## --- Step 2: observed test statistic --------------------------------
  if (is.null(Thatt)) Thatt <- Tstat(Dhatt, bi_out = bi_out)
  .vcat(verbose, 1, "\n  T_obs = ", round(Thatt, 6))

  ## --- Retrieve stored objects ----------------------------------------
  data       <- attr(Dhatt, "internal")$data   ## already sorted by group inside Dmethod
  data       <- data[order(data[[group]]), ]
  R_list     <- Dhatt$R
  Ypred_list <- Dhatt$Ypred
  Dhat       <- Dhatt$Dhat
  ids        <- names(R_list)
  N          <- length(ids)
  re_names   <- colnames(Dhat)
  response      <- as.character(Expr[[2]])


  ## --- Step 3: compute adjusted residuals under H0 --------------------
  ## For subset test: residualise out retained RE via EBLUP
  idx_out  <- if (!is.null(bi_out)) which(re_names %in% bi_out) else seq_along(re_names)
  idx_keep <- setdiff(seq_along(re_names), idx_out)
  Newdata  <- .compute_Ystar_subset(data, group, Ypred_list, R_list, Dhat,
                                    Dhatt$Sigma2, ids, response, idx_keep)


  ## --- Step 4: permutation loop ----------------------------------------


  .vcat(verbose, 1, "\n  Running ", nperm, " permutations ...")
  ## Loop 1: generate all permutations quickly -- no progress message needed
  PYstar <- matrix(NA, nrow(data), nperm)
  max_t  <- max(Newdata$T_pos)

  if (!is.null(seed)) set.seed(seed)

  for (b in seq_len(nperm)) {
    Ystar <- Newdata$Ystar
    for (t in seq_len(max_t)) {
      idx <- which(Newdata$T_pos == t)
      if (length(idx) > 1)
        Ystar[idx] <- sample(Newdata$Ystar[idx], replace = FALSE)
    }
    PYstar[, b] <- Ystar

    ## Progress message every perm_freq permutations
    if (verbose >= 2 && b %% perm_freq == 0)
      .vcat(verbose, 2, "  Permutation ", b, " / ", nperm, " done.")
  }
  Yperm <- sweep(PYstar, 1, Newdata$Pred_H0, FUN = "+")



  ## Re-estimate D* on permuted data
  data_b        <- data
  .perm_counter <- 0L
  Tperm <- apply(Yperm, 2, function(ynew) {
    .perm_counter       <<- .perm_counter + 1L
    data_b[[response]]  <- c(ynew)
    Db <- tryCatch(
      suppressWarnings(
        Dmethod(data_b, Expr, group, random, start,
                method       = method,
                is_permuting = TRUE,
                kappa_max    = kappa_max,
                RR_catof     = RR_catof,
                verbose      = 0)
      ),
      error = function(e) NULL)
    if (is.null(Db) || all(diag(Db$Dhat) * c(Db$Sigma2) <= 0) ) NA_real_ else Tstat(Db, bi_out = bi_out)
  })



  ## --- Step 5: p-value and decision ------------------------------------
  pvalue_num <- mean(Tperm >= Thatt, na.rm = TRUE)
  decision   <- if (pvalue_num < sig_alpha) "Reject H0" else "Do not reject H0"
  pval       <- if (pvalue_num < 0.001) "< 0.001" else pvalue_num

  .vcat(verbose, 1,
        "\n  p-value  = ", pval,
        "\n  Decision : ", decision, "\n")

  ## --- Build permutation histogram plot --------------------------------
  plt <- .perm_hist_plot(Tperm, Thatt, pval, decision, bi_out)

  structure(
    list(
      Decision       = decision,
      pvalue         = pval,
      Tobs           = Thatt,
      Tperm          = Tperm,
      Dhatt          = Dhatt,
      bi_out         = bi_out,
      plot           = plt
    ),
    internal = list(pvalue_num = pvalue_num),
    class = "Dtest"
  )
}


# ============================================================
# Helper: compute Ystar for subset test (EBLUP of retained REs)
# ============================================================
#' @keywords internal

.compute_Ystar_subset <- function(data, group, Ypred_list, R_list,
                                  Dhat, sigma2, ids,
                                  response, idx_keep) {
  #Ystar   <- vector("list", length(ids)); names(Ystar) <- ids
  #Pred_H0 <- vector("list", length(ids)); names(keep_contrib) <- ids
  result <- vector("list", length(ids)); names(result) <- ids

  for (id in ids) {
    Ri   <- R_list[[id]]
    Yi   <- data[data[[group]] == id, response]
    fhat <- Ypred_list[[id]]
    ei   <- Yi - fhat
    ni   <- length(ei)

    ## Contribution of retained REs
    if (length(idx_keep) <= 0){Pred_H0_i <- fhat} else{
      V    <- Ri %*% Dhat %*% t(Ri) + sigma2 * diag(ni)
      Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
      bhat <- Dhat %*% t(Ri) %*% Vinv %*% ei   # k x 1:Full EBLUP

      bhat_keep        <- bhat[idx_keep, , drop = FALSE]
      Ri_keep          <- Ri[, idx_keep, drop = FALSE]
      Pred_H0_i <- fhat + as.numeric(Ri_keep %*% bhat_keep)
    }
    result[[id]] <- data.frame(
      id      = id,
      T_pos   = seq_len(ni),
      Pred_H0 = Pred_H0_i,
      Ystar   = Yi - Pred_H0_i
    )

  }
  do.call(rbind, result)
}


# ============================================================
# Permutation histogram ggplot
# ============================================================
#' @keywords internal
.perm_hist_plot <- function(Tperm, Tobs, pval, decision, bi_out) {
  df  <- data.frame(T = Tperm)

  pval_label <- if (pval == "< 0.001") "< 0.001" else paste0("= ", round(as.numeric(pval), 4))
  title <- if (is.null(bi_out)) { "Permutation test under H0: all random effects" } else {
    paste0("Permutation test under H0: ", paste(bi_out, collapse = ", "))  }

  ggplot2::ggplot(df, ggplot2::aes(x = T)) +
    ggplot2::geom_histogram(bins = 30,
                            fill = "steelblue", colour = "white",
                            alpha = 0.8) +
    ggplot2::geom_vline(xintercept = Tobs,
                        linetype = "dashed", colour = "red",
                        linewidth = 1) +
    ggplot2::annotate("text",
                      x = Tobs, y = Inf,
                      label = paste0(" T[obs]==", round(Tobs, 4)),
                      parse = TRUE,
                      hjust = -0.1, vjust = 1.5,
                      colour = "red", size = 3.5) +
    ggplot2::annotate("text",
                      x = -Inf, y = Inf,
                      label = paste0("p ", pval_label, "\n", decision),
                      hjust = -0.05, vjust = 1.5,
                      size = 3.5) +
    ggplot2::labs(
      title    = title,
      x        = "Test statistic T",
      y        = "Count"
    ) +
    ggplot2::theme_bw()
}

# ============================================================
# Print / summary for Dtest
# ============================================================

#' @export
print.Dtest <- function(x, ...) {
  cat("\nTestREnlme: Permutation test result\n")
  cat("  Tested RE   :",
      if (is.null(x$bi_out)) "All" else paste(x$bi_out, collapse = ", "), "\n")
  cat("  T_obs       :", round(x$Tobs, 6), "\n")
  cat("  p-value     :", x$pvalue, "\n")
  cat("  Decision    :", x$Decision, "\n")
  invisible(x)
}

#' @export
summary.Dtest <- function(object, ...) print(object, ...)
