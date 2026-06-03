#' @keywords internal
.pkg_env <- new.env(parent = emptyenv())

# ============================================================:
# Verbose helper
# ============================================================:

#' @keywords internal
.vcat <- function(verbose, level, ...) {
  if (verbose >= level) message(...)
}

#' @keywords internal
.vperm <- function(verbose, b, nperm, freq = 10) {
  if (verbose >= 2 && (b %% freq == 0 || b == nperm)) {
    message("  Permutation ", b, " / ", nperm, " ...")
  }
}

# ============================================================:
# Nearest positive-definite correction
# ============================================================:

#' @keywords internal
.nearest_pd <- function(M, warn = TRUE, name = "D") {
  if (!matrixcalc::is.positive.semi.definite(M + diag(1e-12, nrow(M)))) {
    if (warn)
      warning("Estimated ", name,
              " is not positive definite; ",
              "replacing with nearest positive-definite matrix.",
              call. = FALSE)
    M <- as.matrix(Matrix::nearPD(M, corr = FALSE,
                                   keepDiag = FALSE,
                                   maxit = 1000)$mat)
  }
  M
}

# ============================================================:
#' Estimate fixed effects by nonlinear least squares
#'
#' Estimates the fixed-effects parameter vector by unweighted NLS (first
#' attempt via \code{\link[stats]{nls}}), or by weighted NLS via
#' \code{\link[stats]{nlminb}} when \code{weights} is supplied or when the
#' unweighted NLS fails. Used both for initial unweighted estimation and for
#' GLS re-estimation at the end of VLS using the estimated variance
#' components to construct weight.
#'
#' @param data A \code{data.frame} containing all model variables.
#' @param Expr A two-sided formula specifying the nonlinear model.
#' @param start A named numeric vector of starting values.
#' @param weights Either \code{NULL} (default, unweighted), a numeric vector
#'   of weights, or a named list of per-subject inverse covariance matrices
#'   (used for GLS re-estimation after variance components are estimated).
#' @param group Character. Name of the grouping variable, required when
#'   \code{weights} is a list.
#' @param verbose Integer (0, 1, or 2). Default \code{1}.
#'
#' @return A named numeric vector of estimated fixed effects.
#'
#' @examples
#' data("Theoph", package = "nlme")
#' d      <- as.data.frame(Theoph)
#' Expr   <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
#'             (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
#'             (exp(ai2) - exp(ai3))
#' start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
#' Beta_hat(d, Expr, start)
#'
#' @export
Beta_hat <- function(data, Expr, start, weights = NULL, group = NULL,
                     verbose = 1) {
  
  ## Unweighted NLS as first try (skip if weights supplied)
  BOLS <- NULL
  if (is.null(weights)) {
    .vcat(verbose, 2, "  Beta_hat: trying unweighted NLS ...")
    Beta_hat_fct <- function(data,Expr, start){
      Expr0 <- as.character(Expr)
      Expr1 <- Expr0[grep(paste(all.vars(Expr)[-1],collapse="|"), Expr0)]
      Expr2 <- parse(text = Expr1)
      Xs <- all.vars(Expr)[-1][!all.vars(Expr)[-1]%in% names(start)]
      outcome <-   all.vars(Expr)[1]
      Y_all <- data[,outcome]
      YPred <- eval(parse(text = Expr0), c(as.list(data),start))
      sum(abs(Y_all-YPred)^2)
    }
    BOLS <- tryCatch(nlminb(start,  objective=Beta_hat_fct, data=data, Expr=Expr, control=list(iter.max=200000)),error = function(e) NULL)$par
    BOLS <- if(is.null(BOLS)) 
      tryCatch( coef(nls(Expr, data = data, start = as.list(start), control = nls.control(maxiter = 200000))), error = function(e) NULL) else BOLS
  }
  
  ## Fall back to nlminb (also used for weighted estimation)
  if (is.null(BOLS)) {
    .vcat(verbose, 2, "  Beta_hat: using nlminb ...")
    
    Beta_hat_fct <- function(par, data, Expr, weights, group) {
      Expr0 <- as.character(Expr)[3]
      Yvar  <- all.vars(Expr)[1]
      Y_all <- data[[Yvar]]
      for (ii in names(par)) data[[ii]] <- par[ii]
      preds  <- eval(parse(text = Expr0), data)
      resids <- Y_all - preds
      
      if (is.null(weights)) {
        return(sum(resids^2))
      } else if (!is.list(weights)) {
        return(sum(weights * resids^2))
      } else if (is.list(weights) && !is.null(group)) {
        ids   <- unique(data[[group]])
        total <- 0
        for (i in seq_along(ids)) {
          idx   <- which(data[[group]] == ids[i])
          r_i   <- matrix(resids[idx], ncol = 1)
          total <- total + drop(t(r_i) %*% weights[[i]] %*% r_i)
        }
        return(total)
      } else {
        stop("weights must be NULL, a vector, or a list with group specified.",
             call. = FALSE)
      }
    }
    
    OptBeta <- try(nlminb(start,
                          objective = Beta_hat_fct,
                          data      = data,
                          Expr      = Expr,
                          weights   = weights,
                          group     = group,
                          control   = list(iter.max = 200000)),
                   silent = TRUE)
    
    BOLS <- if (!inherits(OptBeta, "try-error")) {
      OptBeta$par
    } else {
      .vcat(verbose, 1, "  Beta_hat: nlminb failed, returning NA.")
      rep(NA_real_, length(start))
    }
    names(BOLS) <- names(start)
  }
  
  BOLS
}
# ============================================================================:
#' Second-stage weighted GLS beta estimate for MM methods
#'
#' Computes \eqn{\hat\beta_{\text{MM}}} using the two-stage generalised least
#' squares formula (Davidian & Giltinan 1995, eq. 8.22).
#'
#' @param Aai       List with `Ai` (per-subject design matrices) and `a0i`
#'   (per-subject coefficient matrix), as returned by [Dhat_MM_base()].
#' @param k         Integer. Number of random-effect parameters.
#' @param Ti        List of per-subject scaled covariance matrices.
#' @param Dhat      Estimated D matrix.
#' @param start     Named numeric starting-value vector.
#' @param MM_IDworks Integer vector of retained subject indices.
#' @param random    Named vector of random-effects formulas.
#'
#' @return Named numeric vector of fixed-effect estimates.
#' @keywords internal
.Beta_MM_TS <- function(Aai, k, Ti, Dhat, start, MM_IDworks = NULL, random) {
  Ai   <- Aai$Ai
  a0i  <- Aai$a0i
  if (is.null(MM_IDworks)) MM_IDworks <- seq_len(nrow(a0i))
  
  ATDA_i  <- lapply(MM_IDworks, function(x) t(Ai[[x]]) %*% solve(Ti[[x]] + Dhat) %*% Ai[[x]])
  ATDa_i  <- lapply(MM_IDworks, function(x) t(Ai[[x]]) %*% solve(Ti[[x]] + Dhat) %*% a0i[x, ])
  SATDA_i <- Reduce("+", ATDA_i)
  SATDa_i <- Reduce("+", ATDa_i)
  Beta    <- c(solve(SATDA_i) %*% SATDa_i)
  
  # Recover parameter names from random formulas
  names(random) <- gsub("~.*| ", "", random)
  random <- random[names(start)]
  names0 <- strsplit(gsub(".*~| ", "", random), split = "\\+")
  names1 <- unlist(lapply(names0, function(x) gsub("\\*.*| ", "", x[-length(x)])))
  names(Beta) <- names1
  Beta
}

# ============================================================================:
# Within-cluster variance estimators

#' Pooled within-cluster variance for MM methods
#'
#' Computes the pooled within-cluster residual variance estimate
#' \eqn{\hat\sigma^2 = \text{SSE} / \sum_i (n_i - k)}, where *k* is the
#' number of parameters in `start`.
#'
#' @param data  A `data.frame`.
#' @param start Named numeric starting-value vector (length = *k*).
#' @param Expr  First-stage formula.
#' @param group Character. Grouping column name.
#'
#' @return A scalar estimate of \eqn{\sigma^2}.
#' @keywords internal


.Sigma2hat_MM <- function(data, start, Expr, group){
  #1. data and Expression
  #---~---~---~---~---~-
  outcome <-   all.vars(Expr)[1]
  lhs <- as.character(Expr)[3]
  k <- length(start)#ncol(start)
  Nt <- nrow(data)
  id <- data[,group]#id
  uid <- unique(id)
  nud <- length(uid)
  SSE <-  dfpool <- 0
  SSE1 <- NULL
  
  #2.Loop through clusters
  #---~---~---~---~---
  data[,names(start)] <-NULL
  for(i in 1:nud) {
    data_i <- data[id==uid[i],]
    Bhat_i <- Beta_hat(data_i, Expr, start)
    ypred <- eval(parse(text=lhs),c(as.list(data_i),Bhat_i) )
    y_i <- data_i[,outcome]
    SSE1 <- c(SSE1, sum((y_i-ypred)^2))
    SSE <- SSE+sum((y_i-ypred)^2)
    dfpool <- dfpool+length(y_i)-k #sum n_i-k equ:8.17 page443
  }
  #3. Estimate within-cluster variance (sigma^2)
  #---~---~---~---~---~---~---~---~---~---~---
  Sigma2hat <- SSE/dfpool
  return(Sigma2hat)
}



# ============================================================:
#' Compute Jacobian matrix and fixed-effects predicted values
#'
#' Evaluates the Jacobian \eqn{R_i(\hat\theta_0)} of the nonlinear model
#' function with respect to the random-effects parameters, and computes
#' fitted values under the fixed-effects-only model.
#'
#' @param data A \code{data.frame}.
#' @param Expr A two-sided formula for the nonlinear model.
#' @param group Character name of the grouping variable.
#' @param random A named list or vector of one-sided formulas mapping each
#'   parameter to its random-effect expression, e.g.
#'   \code{c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2")}.
#' @param Bhat Named numeric vector of fixed-effects estimates from
#'   \code{\link{Beta_hat}}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{R}}{Named list of per-subject Jacobian matrices
#'       \eqn{R_i(\hat\theta_0)}.}
#'     \item{\code{Ypred}}{Named list of per-subject predicted vectors
#'       \eqn{f_i(A_i\hat\beta, \hat\gamma)}.}
#'     \item{\code{residuals}}{Named list of per-subject residual vectors
#'       \eqn{\hat e_i = Y_i - f_i}.}
#'   }
#'
#' @export
ZandYPred <- function(data, Expr0, Expr, group, random, Bhat0 = NULL,Bhat = NULL, start = NULL,Dnames) {
  if (is.null(Bhat0))  Bhat0 <- Beta_hat(data = data, Expr = Expr0, start = start[Dnames])
  if (is.null(Bhat))  Bhat <- Beta_hat(data = data, Expr = Expr, start = start)
  id    <- data[[group]]
  Yvar  <- as.character(Expr[[2]])
  rhs   <- as.character(Expr)[3]
  rhs0  <- as.character(Expr0)[3]
  ldata <- c(as.list(data), as.list(Bhat))
  ldata0<- c(as.list(data), as.list(Bhat0))

  ## Design matrices for FE (K) and RE (R)
  K0 <- lapply(names(Bhat), function(x) eval(D(parse(text = rhs), x), ldata))
  R0 <- lapply(Dnames,      function(x) eval(D(parse(text = rhs0), x), ldata0))
  K  <- as.matrix(setNames(data.frame(Reduce("cbind", K0)), names(Bhat)))
  R  <- as.matrix(setNames(data.frame(Reduce("cbind", R0)), Dnames))
  
  ## Fitted values and residuals
  data$yp <- eval(parse(text = rhs), ldata)
  data$y  <- data[[Yvar]]
  data$e  <- data$y - data$yp
  
  ## Split into per-subject lists
  ids        <- unique(id)
  K_list     <- lapply(ids, function(i) K[id == i, , drop = FALSE])
  R_list     <- lapply(ids, function(i) R[id == i, , drop = FALSE])
  Ypred_list <- lapply(ids, function(i) data$yp[id == i])
  resid_list <- lapply(ids, function(i) data$e[id == i])
  
  names(K_list) <- names(R_list) <-
    names(Ypred_list) <- names(resid_list) <- ids
  
  list(R = R_list, K = K_list, Ypred = Ypred_list, residuals = resid_list)
}

# ============================================================:
#' Build model expressions for MM/MMF estimation
#'
#' Constructs the set of model expressions needed by \code{\link{MM_base}}:
#' the full expression with random effects substituted, the second-stage
#' expression incorporating subject-level covariates from \code{Q}, and
#' the \code{nlsList}-style formula for per-subject NLS fitting.
#'
#' @param data A \code{data.frame}.
#' @param group Character. Name of the grouping variable.
#' @param random Named character vector of random-effects formulas,
#'   e.g. \code{c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2")}.
#' @param Expr A two-sided formula for the nonlinear model.
#' @param start Named numeric vector of starting values.
#' @param Q One-sided formula for the second-stage design. Default \code{~1}.
#'
#' @return A list with components \code{Expr}, \code{Expr1}, \code{Expr2},
#'   \code{Expr_MM_all0}, \code{random0}, \code{start}, and
#'   \code{start_MM_all}.
#'
#' @keywords internal
Expressions <- function(data, group, random, Expr, start) {
  
  ## 1. Extract Q: subject level design 
  Q0 <- all.vars(as.formula(random[[1]]))
  Q1 <- intersect(Q0, names(data))
  Q2 <- if (length(Q1)==0) {"~1"} else {paste("~", paste(Q1,collapse = "+"))}
  Q <- as.formula(Q2)
  
  #Fixed
  Fixed  <- names(start)
  Random <- trimws(sub("~.*", "", random))
  QRHS   <- setNames(trimws(sub(".*~", "", random)), Random)
  QRHS1  <- sub("\\+[^+]*$", "", QRHS)   # drop last term (random part)

  ## Full expression with random effects substituted
  Expr1      <- as.character(Expr)
  Expr1[[3]] <- mgsub::mgsub(Expr1[[3]], names(QRHS1), paste0("(", QRHS, ")"))
  Expr1      <- Reduce("paste", Expr1[c(2, 1, 3)])

  ## Reduced expression (first-stage fixed only, or second-stage combined)
  if (identical(deparse(Q), "~1")) {
    Expr2     <- Expr
    New_start <- start
  } else {
    QRHS2      <- paste0("(", QRHS1, ")")
    Expr2      <- as.character(Expr)
    Expr2[[3]] <- mgsub::mgsub(Expr2[[3]], names(QRHS1), QRHS2)
    Expr2      <- Reduce("paste", Expr2[c(2, 1, 3)])
    newPar     <- setdiff(all.vars(as.formula(Expr2)), c(Fixed, names(data)))
    start0     <- setNames(rnorm(length(newPar)), newPar)
    New_start  <- c(start, start0)
    New_start  <- New_start[names(New_start) %in% all.vars(as.formula(Expr2))]
  }

  ## nlsList-style formula: outcome ~ rhs | group
  Expr00       <- as.character(Expr)
  Expr_MM_all0 <- paste(Expr00[[3]], "|", group)
  Expr00[[3]]  <- Expr_MM_all0
  Expr_MM_all0 <- Reduce("paste", Expr00[c(2, 1, 3)])

  N            <- length(unique(data[, group]))
  NStart       <- setNames(names(start), names(start))
  start_MM_all <- lapply(NStart,
    function(x) if (x %in% Random) rep(start[x], N) else start[x])

  list(
    Expr         = Expr,
    Expr1        = Expr1,
    Expr2        = Expr2,
    Expr_MM_all0 = Expr_MM_all0,
    random0      = Random,
    start        = New_start,
    start_MM_all = start_MM_all
  )
}
