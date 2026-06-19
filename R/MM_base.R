#' Compute first-stage quantities for MM and MMF estimators
#'
#' Performs the first-stage nonlinear least squares (NLS) fit for each
#' subject individually, and computes the per-subject Jacobian matrices
#' \eqn{R_i^T R_i} and condition numbers \eqn{\kappa_i}. The returned
#' object can be passed to \code{\link{Dmethod}} via the \code{MM_base_obj}
#' argument to avoid recomputing the first stage when calling both
#' \code{method = "MM"} and \code{method = "MMF"} on the same data.
#'
#' @param data A \code{data.frame} containing all variables.
#' @param Expr A two-sided formula specifying the nonlinear model
#'   \eqn{f_i(a_i, \gamma)}. The left-hand side is the response variable
#'   and the right-hand side defines the nonlinear function using the
#'   subject-specific parameter names (e.g., \code{ai1}, \code{ai2},
#'   \code{ai3}) that appear in \code{start} and \code{random}.
#' @param group Character. Name of the grouping (subject) variable.
#' @param random A character vector of two-sided formula strings, one per
#'   parameter in \eqn{a_i = A_i\beta + b_i}, mapping each subject-specific
#'   parameter (left-hand side, matching \code{Expr} and \code{start})
#'   to its fixed-effects expression plus random effect (right-hand side).
#'   For example, \code{c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2",
#'   "ai3 ~ B3 + bi3")} specifies that \code{ai1 = B1 + bi1},
#'   \code{ai2 = B2 + bi2}, \code{ai3 = B3 + bi3}.
#' @param start A named numeric vector of starting values for all \code{Expr}
#'   parameters. Names must match those used in \code{Expr} (the
#'   subject-specific parameter names, e.g., \code{ai1}, \code{ai2},
#'   \code{ai3}). Good starting values can be obtained from a
#'   preliminary call to \code{nls()} on the pooled data.
#' @param kappa_max Positive numeric. Subjects whose per-subject Jacobian
#'   condition number exceeds this threshold are excluded from MM/MMF
#'   second-stage estimation. Default \code{1e4}.
#' @param RR_catof Character. Exclusion criterion: \code{"kappa"} (default)
#'   uses the condition-number threshold \code{kappa_max};
#'   or a user-specified numeric threshold on the \eqn{R_i^T R_i}.
#' @param verbose Integer (0, 1, or 2). \code{0} = silent; \code{1} =
#'   prints a summary of subjects retained and excluded (default);
#'   \code{2} = additionally prints per-subject convergence messages.
#'
#' @return A list of class \code{"MM_base"} with components:
#'   \describe{
#'     \item{\code{data}}{The (re-ordered, group-renumbered) \code{data.frame}
#'       used internally.}
#'     \item{\code{group}}{Name of the grouping variable.}
#'     \item{\code{start}}{The starting values used.}
#'     \item{\code{Beta_nls}}{Named numeric vector of pooled (first-pass)
#'       fixed-effects estimates from \code{\link{Beta_hat}}.}
#'     \item{\code{random}}{The \code{random} formulas, aligned with
#'       \code{start}.}
#'     \item{\code{re_names}}{Named character vector mapping random-effect
#'       names to subject-specific parameter names.}
#'     \item{\code{Tmatrix}}{Matrix of per-subject scaled covariance terms
#'       (one row per subject, \code{k^2} columns).}
#'     \item{\code{Ti}}{List of per-subject \eqn{k \times k} covariance
#'       matrices derived from \code{Tmatrix}.}
#'     \item{\code{ddi}}{List of per-subject \eqn{k \times k} second-stage
#'       deviation outer-product matrices (eq.~8.20).}
#'     \item{\code{Lamda}}{Minimum eigenvalue of the standardised scatter
#'       matrix, used to bias-correct the MM/MMF variance estimates.}
#'     \item{\code{sigma2}}{Pooled residual variance across retained
#'       subjects.}
#'     \item{\code{Beta_GB}}{Named numeric vector of non-random-effect
#'       fixed-effects estimates.}
#'     \item{\code{m}}{Total number of second-stage fixed-effects
#'       parameters (\code{k * ncol(Q)}).}
#'     \item{\code{k}}{Number of random-effects parameters.}
#'     \item{\code{id}, \code{uid}, \code{nud}}{Subject index vector, unique
#'       subject IDs, and number of unique subjects.}
#'     \item{\code{IDworks}}{Integer vector of subject IDs retained after
#'       filtering.}
#'     \item{\code{Q}}{Subject-level design matrix.}
#'     \item{\code{Aai}}{List with elements \code{Ai} (per-subject design
#'       matrices) and \code{a0i} (per-subject parameter estimates).}
#'     \item{\code{kappa}}{Named numeric vector of per-subject condition
#'       numbers.}
#'     \item{\code{IDconvegence}}{List with elements \code{ID_all},
#'       \code{ID_used}, and \code{ID_supressed} (itself a list with
#'       \code{highRR} and \code{Unconverged}).}
#'   }
#'
#' @examples
#' d      <- as.data.frame(Theoph)
#' Expr   <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
#'             (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
#'             (exp(ai2) - exp(ai3))
#' start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
#' random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")
#' mb <- MM_base(d, Expr, group = "Subject", random = random, start = start)
#'
#' @importFrom nlme nlsList
#' @importFrom MASS ginv
#' @importFrom expm sqrtm
#' @importFrom stats nls nls.control coef model.matrix as.formula
#' @export
MM_base <- function(data, Expr, group, random, start,
                    kappa_max = 1e4,
                    RR_catof = "kappa",
                    verbose = 1) {


  # RR_catof validation
  if (is.numeric(RR_catof)) {if (length(RR_catof) != 1 || RR_catof <= 0)
    stop("RR_catof must be a positive numeric value.")
  } else {RR_catof <- match.arg(RR_catof, "kappa")}

  ## 0. Align random with start
  names(random) <- gsub("~.*| ", "", random)
  random <- random[intersect(names(random), names(start))]
  re_names <- setNames(gsub("~.*| ","", random),gsub(".*\\+| ","", random))


  ## 1. Extract Q: subject level design
  Q0 <- all.vars(as.formula(random[[1]]))
  Q1 <- intersect(Q0, names(data))
  Q2 <- if (length(Q1)==0) {"~1"} else {paste("~", paste(Q1,collapse = "+"))}
  Q <- as.formula(Q2)


  ## 1. Data preparation
  data[[group]] <- as.numeric(as.factor(data[[group]]))
  data          <- data[order(data[[group]]), ]
  rownames(data)<- seq_len(nrow(data))
  k       <- length(random)
  Nt      <- nrow(data)
  m0      <- length(start)
  id      <- data[[group]]
  uid     <- unique(id)
  nud     <- length(uid)
  outcome <- all.vars(Expr)[1]
  Idcum   <- cumsum(table(id))
  Q1      <- as.matrix(model.matrix(Q, data = data)[Idcum, , drop = FALSE])
  rownames(Q1) <- seq_len(nrow(Q1))
  m <- k * ncol(Q1)

  ## 2. Overall Bhat, add non-RE fixed effects to data
  Bhat     <- Beta_hat(data = data, Expr = Expr, start = start)
  Beta_all <- setdiff(names(Bhat), names(random))
  for (ii in Beta_all) data[[ii]] <- Bhat[[ii]]

  ## 3. Build per-subject NLS expressions and fit
  Formulars  <- Expressions(data = data, Expr = Expr, group = group,
                            random = random, start = start)
  Expr_nls   <- as.formula(paste(Formulars$Expr_MM_all0, collapse = " "))
  start_RE   <- start[names(random)]

  .vcat(verbose, 2, "MM_base: fitting per-subject NLS (", nud, " subjects) ...")

  oi       <- try(nlme::nlsList(Expr_nls, data,
                                start   = start_RE, warn.nls=FALSE,
                                control = nls.control(maxiter = 200000)),  silent = TRUE)

  if (inherits(oi, "try-error")) {
    stop("MM_base: first-stage per-subject NLS fitting failed entirely. ",
         "Check starting values or model specification.", call. = FALSE)
  }

  covpar_all <- summary(oi)$cov.unscaled
  COEF_all   <- coef(oi)

  ## 4. Initial sigma on converged subjects
  COEF_all0   <- COEF_all[rowSums(is.na(COEF_all)) != ncol(COEF_all), ]
  IDRR        <- as.integer(rownames(COEF_all0))
  unconverged <- setdiff(uid, IDRR)
  Newdata     <- data[data[[group]] %in% IDRR, ]
  s2pool_init <- .Sigma2hat_MM(data = Newdata, COEFs = COEF_all,
                               Expr = Expr, group = group)



  ## 5. Loop over converged subjects - collect a_i, T_i, A_i
  a0_i    <- matrix(NA_real_, nrow = nud, ncol = k)
  uT_i    <- matrix(NA_real_, nrow = nud, ncol = k^2)
  IDworks <- NULL
  AA      <- NULL
  Aa      <- NULL
  Ai      <- vector("list", nud)
  Scov    <- 0
  high_kappa <- integer(0)
  kappa_vec <- setNames(rep(NA_real_, length(IDRR)), as.character(IDRR))

  for (i in IDRR) {
    covpar <- covpar_all[i, , , drop = FALSE][1, , ]

    ## Filtering
    Include_i <- TRUE
    if (identical(RR_catof, "kappa")) {
      FIT       <- oi[[as.character(i)]]
      kappa_val <- kappa(crossprod(FIT$m$gradient()))
      kappa_vec[as.character(i)] <- kappa_val
      if (!is.null(kappa_max) && kappa_val > kappa_max) {
        Include_i <- FALSE
        .vcat(verbose, 2, "  Subject ", i,
              ": excluded (kappa = ", round(kappa_val, 1), ").")
      }
    } else if (is.numeric(RR_catof)) {
      HighRR_cov <- sum(abs(covpar)) * s2pool_init
      if (HighRR_cov >= RR_catof) {
        Include_i <- FALSE
        .vcat(verbose, 2, "  Subject ", i,
              ": excluded (HighRR = ", round(HighRR_cov, 4), ").")
      }
    }

    if (Include_i) {
      a           <- c(unlist(COEF_all[i, ]))
      a0_i[i, ]   <- a
      uT_i[i, ]   <- as.vector(covpar)
      IDworks     <- c(IDworks, i)
      Qi          <- Q1[i, ]
      Ai[[i]]     <- kronecker(diag(1, k, k), t(Qi))
      AA          <- cbind(AA, c(t(Ai[[i]]) %*% Ai[[i]]))
      Aa          <- cbind(Aa, t(Ai[[i]]) %*% matrix(a, k, 1))
      Scov        <- Scov + covpar
    } else {
      high_kappa <- c(high_kappa, i)
    }
  }

  ## 6. Bookkeeping
  IDkept <- list(
    ID_all      = uid,
    ID_used     = IDworks,
    ID_supressed = list(
      highRR      = high_kappa,
      Unconverged = unconverged
    )
  )

  ## 7. Refined sigma and T matrices on retained subjects only
  data_used <- data[data[[group]] %in% IDworks, ]
  s2pool    <- .Sigma2hat_MM(data = data_used, COEFs = COEF_all,
                             Expr = Expr, group = group)
  Tmatrix   <- uT_i * s2pool
  Ti        <- lapply(seq_len(nrow(Tmatrix)),
                      function(x) matrix(Tmatrix[x, ], k, k))

  ## 8. Second-stage deviations di (eq. 8.20)
  invSAA  <- solve(matrix(rowSums(AA), m, m))
  SAa     <- matrix(rowSums(Aa), m, 1)
  a_names <- names(COEF_all[IDworks[1], ])

  di <- sapply(IDworks, function(x)
    matrix(a0_i[x, ], k, 1) - Ai[[x]] %*% invSAA %*% SAa)
  rownames(di) <- a_names

  ddi0 <- lapply(seq_len(length(IDworks)),
                 function(x) matrix(c(di[, x] %*% t(di[, x])), k, k))
  ddi  <- lapply(ddi0, function(x) {dimnames(x) <- list(a_names, a_names); x })


  ## 9. Lambda (minimum eigenvalue of standardised scatter)
  Sddi  <- Reduce("+", ddi)
  dimnames(Sddi) <- list(a_names, a_names)
  Tmat  <- matrix(colSums(Tmatrix, na.rm = TRUE), k, k)
  TddT  <- solve(expm::sqrtm(t(Tmat))) %*% Sddi %*% solve(expm::sqrtm(Tmat))
  Lamda <- min(eigen(TddT)$values)


  #10 other object to return
  a0_list        <- list(Ai=Ai, a0i=a0_i)
  Beta_GB        <- Bhat[Beta_all]
  id_used        <- IDworks
  ## Print summary
  .vcat(verbose, 1,
        "\nMM_base summary:",
        "\n  Total Observations : ", Nt,
        "\n  Total subjects     : ", nud,
        "\n  Subjects used      : ", length(id_used),
        "\n  High-kappa excluded: ", length(high_kappa),
        "\n  Non-converged      : ", length(unconverged))

  if (!all(is.na(kappa_vec)))
    .vcat(verbose, 1,
          "\n  Condition number range: [",
          round(min(kappa_vec, na.rm = TRUE), 1), ", ",
          round(max(kappa_vec, na.rm = TRUE), 1), "]")


  structure(
    list(
      data         = data,
      group        = group,
      start        = start,
      Beta_nls     = Bhat,
      random       = random,
      re_names     = re_names,
      Tmatrix      = Tmatrix,
      Ti           = Ti,
      ddi          = ddi,
      Lamda        = Lamda,
      sigma2       = s2pool,
      Beta_GB      = Beta_GB,
      m            = m,
      k            = k,
      id           = id,
      uid          = uid,
      nud          = nud,
      IDworks      = IDworks,
      Q            = Q1,
      Aai          = a0_list,
      kappa        = kappa_vec,
      IDconvegence = list(
        ID_all    = as.character(uid),
        ID_used   = as.character(id_used),
        ID_supressed = list(
          highRR      = high_kappa,
          Unconverged = unconverged
        )
      )
    ),
    class = "MM_base"
  )
}

#' @export
print.MM_base <- function(x, ...) {
  cat("MM_base object\n")
  cat("  Total subjects     :", length(x$IDconvegence$ID_all), "\n")
  cat("  Subjects used      :", length(x$IDconvegence$ID_used), "\n")
  cat("  High-kappa excluded:", length(x$IDconvegence$ID_supressed$highRR), "\n")
  cat("  Non-converged      :", length(x$IDconvegence$ID_supressed$Unconverged), "\n")
  cat("  Pooled sigma^2     :", round(x$sigma2, 6), "\n")
  invisible(x)
}
