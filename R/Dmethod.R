#' Non-parametric approach to estimate variance components in a linear
#' and nonlinear mixed-effects model
#'
#' Computes the scaled variance covariance matrix \eqn{\hat D_*} and the
#' error variance \eqn{\hat\sigma^2} using one of three nonparametric
#' estimators: Variance Least Squares (\code{"VLS"}), Method of Moments (\code{"MM"}),
#' or Method of Moments with First-Order Approximation (\code{"MMF"}).
#' The method also estimate weighted fixed effects
#'
#' @param data A \code{data.frame} containing all model variables.
#' @param Expr A two-sided formula specifying the nonlinear model
#'   \eqn{f_i(a_i, \gamma)}. The left-hand side is the response variable
#'   and the right-hand side defines the nonlinear function using the
#'   subject-specific parameter names (e.g., \code{ai1}, \code{ai2},
#'   \code{ai3}) that appear in \code{start} and \code{random}.
#' @param group Character. Name of the grouping variable in \code{data}.
#' @param random A character vector of two-sided formula strings, one per
#'   parameter in \eqn{a_i = A_i\beta + b_i}, mapping each subject-specific
#'   parameter (left-hand side, matching \code{Expr} and \code{start})
#'   to its fixed-effects expression plus random effect (right-hand side).
#'   For example, \code{c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2",
#'   "ai3 ~ B3 + bi3")} specifies that \code{ai1 = B1 + bi1},
#'   \code{ai2 = B2 + bi2}, \code{ai3 = B3 + bi3}.
#' @param start A named numeric vector of starting values for all parameters
#'   in \code{Expr}. Names must match those used in \code{Expr} (e.g.,
#'   \code{ai1}, \code{ai2}, \code{ai3}).  If \code{NULL} (the default),
#'    \code{Dmethod} attempts to compute starting values automatically using
#'   \code{nls.multstart::nls_multstart()}, searching over multiple
#'   initial values within a specified range (e.g. \eqn{\pm 10}) for each
#'   parameter. If this step fails, provide \code{start} manually, or fit
#'   \code{nls_multstart()} or \code{nls()} separately with different
#'   starting values, search bounds, or optimisation settings, and use the
#'   resulting coefficients as starting values.
#' @param method Character. One of \code{"VLS"} (default),
#'   \code{"MM"}, or \code{"MMF"} estimation method.
#' @param MM_base_obj  An optional pre-computed object of class \code{"MM_base"}
#'   returned by \code{\link{MM_base}}. When \code{NULL} (default) and
#'   \code{method} is \code{"MM"} or \code{"MMF"}, \code{MM_base()} is
#'   called internally. Supplying a pre-computed object avoids repeating the
#'   expensive first-stage NLS fits when calling both \code{"MM"} and
#'   \code{"MMF"} on the same data.
#' @param kappa_max Positive numeric. Condition-number threshold for
#'   excluding subjects in MM/MMF. Default \code{1e4}.
#' @param RR_catof Exclusion criterion passed to \code{\link{MM_base}}.
#'   Either \code{"kappa"} (default), which uses the condition-number
#'   threshold \code{kappa_max}, or a user-specified numeric threshold
#'   applied directly.
#' @param Beta_nls Optional named numeric vector of pre-computed fixed-effects
#'   estimates. If \code{NULL} (default), estimated internally via
#'   \code{\link{Beta_hat}}.
#' @param verbose Integer (0, 1, or 2). \code{0} = silent; \code{1} =
#'   summary messages (default); \code{2} = full progress.
#' @param is_permuting Logical. Internal flag set to \code{TRUE} during the
#'   permutation procedure, where the GLS-weighted refit of fixed effects
#'   is not needed and \code{Beta_nls} is used instead. Default \code{FALSE}.
#'
#' @return An object of class \code{"Dmethod"}, a list with components:
#'   \describe{
#'     \item{\code{Dhat}}{The estimated scaled covariance matrix of random effects
#'       \eqn{\hat D_*}.}
#'     \item{\code{Sigma2}}{The estimated error variance \eqn{\hat\sigma^2}.}
#'     \item{\code{Beta}}{The estimated fixed-effects vector \eqn{\hat\beta}.}
#'     \item{\code{method}}{The estimation method used.}
#'     \item{\code{R}}{Per-subject Jacobian matrices \eqn{R_i(\hat\theta_0)}.}
#'     \item{\code{Ypred}}{Per-subject fitted values under fixed effects only.}
#'     \item{\code{residuals}}{Per-subject residuals \eqn{\hat e_i}.}
#'     \item{\code{IDconvegence}}{Subject inclusion/exclusion summary (MM/MMF
#'       only; \code{NULL} otherwise).}
#'       }
#'   Additional internal objects are stored as an \code{"internal"}
#'   attribute (\code{attr(object, "internal")}) and used by package
#'   methods. They are not part of the user interface.
#'
#' @examples
#' d      <- as.data.frame(Theoph)
#' Expr   <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
#'             (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
#'             (exp(ai2) - exp(ai3))
#' start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
#' random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")
#' DVLS   <- Dmethod(d, Expr, group = "Subject",
#'                   random = random, start = start)
#' ## method defaults to VLS; explicit: method = "VLS"
#' DVLS[c("Dhat", "Sigma2", "Beta")]
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
#' \emph{Mixed Models: Theory and Applications with R} (2nd ed.). John Wiley & Sons
#'
#' Drikvandi, R., Verbeke, G., Khodadadi, A. and Nia, V. P. (2013).
#' Testing multiple variance components in linear mixed-effects models.
#' \emph{Biostatistics}, \strong{14 (1)}, 144--159.
#'
#' @importFrom MASS ginv
#' @importFrom utils tail
#' @importFrom Matrix bdiag
#' @importFrom nls.multstart nls_multstart
#' @importFrom stats as.formula coef setNames
#' @export
Dmethod <- function(data, Expr, group, random,
                    start       = NULL,
                    method      = c("VLS", "MM", "MMF"),
                    MM_base_obj = NULL,
                    kappa_max   = 1e4,
                    RR_catof    = "kappa",
                    Beta_nls    = NULL,
                    verbose     = 1,
                    is_permuting = FALSE) {



  method   <- match.arg(method)
  .vcat(verbose, 1, "\nDmethod: estimating variance components (method = ", method, ") ...")


  ## --- Automatic starting values -----------------------
  if (is.null(start)) start <- .auto_start(data, Expr )
  ## --- Reorder start to follow random's canonical order, then any
  ##     remaining (pure fixed-effect) names -----------------------
  random_names <- gsub("~.*| ", "", random)
  STnames      <- c(random_names, setdiff(names(start), random_names))
  start        <- start[STnames]
  Start_orig <-  start

  ## --- Map random, Expr and start from ai to B names -----------------------
  names0  <- strsplit(gsub(".*~| ", "", random), split = "\\+")
  names1  <- unlist(lapply(names0, function(x) gsub("\\*.*| ", "", x[1])))
  names(names1) <- gsub("~.*| ", "", random)
  ## Transform Expr: ai1->B1, ai2->B2, ai3->B3
  expr_chr <- paste(deparse(Expr), collapse = " ")
  for (i in seq_along(names1)) {
    expr_chr <- gsub(paste0("\\b", names(names1)[i], "\\b"), names1[i], expr_chr)
  }
  Expr <- as.formula(expr_chr)
  ## Transform random LHS: ai1->B1, ai2->B2, ai3->B3
  for (i in seq_along(names1)) {
    random[i] <- gsub(".*~", paste(names1[i], "~"), random[i])
  }
  ## Transform start names: ai1->B1, ai2->B2, ai3->B3
  names2 <- names(start) %in% names(names1)
  names(start)[names2] <- names1[names(start[names2])]

  ## --- Extract Q: subject level design --------------------------------
  Q0 <- all.vars(as.formula(random[[1]]))
  Q1 <- intersect(Q0, names(data))
  Q2 <- if (length(Q1)==0) {"~1"} else {paste("~", paste(Q1,collapse = "+"))}
  Q <- as.formula(Q2)


  ## --- organise data --------------------------------------------------
  #keep original group variable as it is
  grp_orig <- tail(make.unique(c(names(data), ".group_original")), 1)
  data[[grp_orig]] <- data[[group]]
  #transform to numeric
  data[[group]] <- as.numeric(as.factor(data[[group]]))
  data          <- data[order(data[[group]]), ]   ## subject-major order


  ## --- fixed effects --------------------------------------------------
  Formulars  <- Expressions(data = data, Expr = Expr, group = group,
                            random = random, start = start)
  Expr2 <- as.formula(Formulars$Expr2)
  start2 <- Formulars$start
  if (is.null(Beta_nls))Beta_nls <- Beta_hat(data, Expr2, start2, verbose = verbose)


  ## --- Jacobian and predicted values ---------------------------------
  ## Dnames: named character vector mapping each random-effect name
  ## (e.g. "bi1") to its corresponding subject-specific parameter name
  ## (e.g. "B1"), derived from 'random'.
  Dnames <- setNames(gsub("~.*| ","", random),gsub(".*\\+\\s*","", random))
  zp     <- ZandYPred(data, Expr2, group, random, Beta_nls)
  R_list     <- zp$R
  K_list     <- zp$K
  Ypred_list <- zp$Ypred
  resid_list <- zp$residuals
  ids        <- names(R_list)
  N          <- length(ids)
  k          <- ncol(R_list[[1]])
  re_names   <- colnames(R_list[[1]])
  ## NOTE: colnames(R_list[[1]]) holds the same values as Dnames (the
  ## subject-specific parameter names), so this line resolves to the
  ## full Dnames vector whenever every column of R has a matching entry
  ## in Dnames (the normal case). It is not a positional filter; the
  ## final re_names retains Dnames's names() (the bi*-style random
  ## effect labels) used below for labelling Dhat.
  re_names <- Dnames[re_names %in% Dnames]

  ## --- dispatch -------------------------------------------------------.
  result <- switch(method,
                   VLS = .vls_estimator(data, Expr2,R_list, K_list, resid_list,
                                        ids, k, re_names,  start = start2, group = group,
                                        verbose = verbose),
                   MM   = {
                     if (is.null(MM_base_obj ))
                       MM_base_obj  <- MM_base(data, Expr, group, random, start,
                                               kappa_max = kappa_max,
                                               RR_catof = RR_catof, verbose = verbose)
                     .mm_estimator(MM_base_obj , re_names, verbose = verbose)
                   },
                   MMF  = {
                     if (is.null(MM_base_obj ))
                       MM_base_obj  <- MM_base(data, Expr, group, random, start,
                                               kappa_max = kappa_max,
                                               RR_catof = RR_catof, verbose = verbose)
                     .mmf_estimator(MM_base_obj , re_names,  verbose = verbose)
                   }
  )

  .vcat(verbose, 1, "  Done. sigma^2 = ", round(result$Sigma2, 6))

  # GLS Beta hat and Ypred
  Dhat  = result$Dhat
  Sigma2= result$Sigma2

  if(!is_permuting){
    weights_list <- lapply(R_list, function(Ri) {
      V_i <- c(Sigma2) * diag(nrow(Ri)) + Ri %*% Dhat %*% t(Ri)
      solve(V_i) })
    Bhat_GLS <-  Beta_hat(data=data, Expr2, start2, weights = weights_list, group = group)
  }else {Bhat_GLS <-  Beta_nls}# This is not required during permutation procedure
  Beta <- Bhat_GLS[names(start2)]

  # Ypred
  zp   <- ZandYPred(data, Expr2, group, random, Beta)
  Ypred_list <- zp$Ypred
  resid_list <- zp$residuals

  #Objects for public
  Final <-  list(
    Dhat         = Dhat,
    Sigma2       = Sigma2,
    Beta         = Beta,
    method       = method,
    R            = R_list,
    Ypred        = Ypred_list,
    residuals    = resid_list,
    IDconvegence = if (method %in% c("MM", "MMF"))
      MM_base_obj$IDconvegence else NULL,
    MM_base_obj  = if (method %in% c("MM", "MMF")) MM_base_obj  else NULL)

  #Objects for internal use
  structure(
    Final,
    internal = list(
      data         = data,
      Expr         = Expr,
      group        = group,
      grp_orig     = grp_orig,
      random       = random,
      start        = start,
      start_orig   = Start_orig,
      Q            = Q,
      re_names     = re_names),
    class = "Dmethod")
}


# ============================================================;
# VLS internal estimator
# ============================================================;
#' @keywords internal
.vls_estimator <- function(data, Expr, R_list, K_list, resid_list,
                           ids, k, re_names,start, group,
                           verbose = 1) {


  id  <- data[[group]]
  Nt  <- nrow(data)
  m   <- ncol(K_list[[1]])
  N   <- length(ids)
  Y   <- unlist(lapply(ids, function(i) data[[as.character(Expr[[2]])]][id == i]))
  iN  <- MASS::ginv(Reduce("+", lapply(K_list, crossprod)))

  ## Initialise accumulators
  Hmat <- 0; Gmat <- 0; Cvec <- 0; Ze <- 0; ee <- 0

  ## Loop over subjects
  for (i in seq_along(ids)) {
    R_i <- R_list[[ids[i]]]
    X_i <- K_list[[ids[i]]]
    e_i <- resid_list[[ids[i]]]

    P_i  <- X_i %*% iN %*% t(X_i)
    tZPZ <- crossprod(R_i, P_i) %*% R_i
    zei  <- crossprod(R_i, e_i)
    tZX  <- crossprod(R_i, X_i)
    ZtZ  <- crossprod(R_i)

    Cvec <- Cvec + as.vector(ZtZ - tZPZ)
    Gmat <- Gmat + kronecker(tZX, tZX)
    Hmat <- Hmat + kronecker(ZtZ, ZtZ) -
      kronecker(ZtZ, tZPZ) -
      kronecker(tZPZ, ZtZ)
    Ze   <- Ze   + kronecker(zei, zei)
    ee   <- ee   + sum(e_i^2)
  }

  ## VLS system
  HGmat <- Hmat + Gmat %*% kronecker(iN, iN) %*% t(Gmat)
  Mat1  <- cbind(c(Nt - m, Cvec), rbind(t(Cvec), HGmat))
  Mat2  <- rbind(ee, Ze)

  vls <- MASS::ginv(Mat1) %*% Mat2

  ## Sigma: robust for small N
  sigma2 <- if (N < 30) {
    Z_star <- as.matrix(Matrix::bdiag(R_list[ids]))
    X_full <- Reduce("rbind", K_list[ids])
    S      <- cbind(X_full, Z_star)
    PS     <- S %*% MASS::ginv(crossprod(S)) %*% t(S)
    ev     <- eigen(diag(Nt) - PS, only.values = TRUE)$values
    Rank   <- sum(abs(ev) > 1e-10)
    as.numeric(crossprod(Y, (diag(Nt) - PS)) %*% Y / Rank)
  } else {
    max(as.numeric(vls[1]), 0)
  }

  ## Dhat and positive definite correction
  Dhat <- matrix(vls[-1], k, k)
  Dhat <- (Dhat + t(Dhat)) / 2
  Dhat <- .nearest_pd(Dhat, warn = FALSE, name = "D (VLS)")
  rownames(Dhat) <- colnames(Dhat) <- names(re_names)

  .vcat(verbose, 1, "  VLS: done. Sigma2 = ", round(sigma2, 6))

  list(Dhat = Dhat, Sigma2 = sigma2)
}


# ============================================================;
# MM internal estimator
# ============================================================;
#' @keywords internal
.mm_estimator <- function(mb,re_names,  verbose = 1) { #R_list, k, re_names,  data,

  IDkept  <- mb$IDconvegence
  IDworks <- IDkept$ID_used
  if (length(IDworks) == 0)
    stop("MM: no subjects retained after filtering. ",
         "Consider reducing the model complexity.", call. = FALSE)

  data    <- mb$data
  group   <- mb$group
  Ti0     <- Filter(function(x) any(!is.na(x)), mb$Ti)
  ddi     <- mb$ddi
  nud     <- length(ddi)
  k       <- mb$k



  # MM estimate: eq. 8.20 / 8.21
  ddTi <- lapply(seq_along(ddi),
                 function(i) ddi[[i]] - mb$Lamda * Ti0[[i]])
  Dmm  <- Reduce("+", ddTi) / nud
  Dhat <- (Dmm + t(Dmm)) / 2# symetric
  Dhat <- .nearest_pd(Dhat, warn = FALSE, name = "D (MM)")
  rownames(Dhat) <- colnames(Dhat) <- names(re_names)
  #Return
  list(Dhat = Dhat, Sigma2 = mb$sigma2)
}

# ============================================================;
# MMF internal estimator
# ============================================================;
#' @keywords internal
.mmf_estimator <- function(mb, re_names, verbose = 1) {

  IDkept  <- mb$IDconvegence
  IDworks <- IDkept$ID_used
  if (length(IDworks) == 0)
    stop("MMF: no subjects retained after filtering. ",
         "Consider increasing kappa_max or using VLS.", call. = FALSE)

  IDworks <- as.numeric(as.character(IDworks))
  data    <- mb$data
  group   <- mb$group
  Ti0     <- Filter(function(x) any(!is.na(x)), mb$Ti)
  ddi     <- mb$ddi
  nud     <- length(ddi)
  k       <- mb$k
  m       <- mb$m
  Q1      <- mb$Q
  Ti      <- mb$Ti


  # MMF estimate: eq. 8.28
  qtq   <- lapply(IDworks, function(i) Q1[i, ] %*% t(Q1[i, ]))
  QQTi  <- lapply(IDworks, function(i)
    drop(1 - t(Q1[i, ]) %*% solve(Reduce("+", qtq)) %*% Q1[i, ]) * Ti[[i]])
  dQQTi <- lapply(seq_along(QQTi),
                  function(i) ddi[[i]] - mb$Lamda * QQTi[[i]])
  Dmm_FA <- Reduce("+", dQQTi) / (nud - m)

  # Positive definite
  Dhat <- (Dmm_FA + t(Dmm_FA)) / 2
  Dhat <- .nearest_pd(Dhat, warn = FALSE, name = "D (MMF)")
  rownames(Dhat) <- colnames(Dhat) <- names(re_names)

  #Return
  list(Dhat = Dhat, Sigma2 = mb$sigma2)
}

# ============================================================;
# Print and summary methods
# ============================================================;

#' @export
print.Dmethod <- function(x, ...) {
  cat("\nTestREnlme: Variance component estimates\n")
  cat("Method :", x$method, "\n\n")
  cat("Fixed effects (Beta):\n")
  print(round(x$Beta, 6))
  cat("\nError variance (Sigma2):\n")
  cat(" ", round(x$Sigma2, 6), "\n")
  cat("\nVariance-covariance matrix (Dhat):\n")
  print(round(x$Dhat, 6))
  if (!is.null(x$IDconvegence)) {
    cat("\nSubject summary:\n")
    cat("  Total           :", length(x$IDconvegence$ID_all), "\n")
    cat("  Used            :", length(x$IDconvegence$ID_used), "\n")
    cat("  High-kappa excl.:", length(x$IDconvegence$ID_supressed$highRR), "\n")
    cat("  Non-converged   :", length(x$IDconvegence$ID_supressed$Unconverged), "\n")
  }
  invisible(x)
}

#' @export
summary.Dmethod <- function(object, ...) {
  print(object, ...)
}
