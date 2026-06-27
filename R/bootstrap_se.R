#' Bootstrap standard errors for variance component estimates
#'
#' Computes bootstrap standard errors for the fixed effects \eqn{\hat\beta},
#' random effect variance components \eqn{\hat D_*}, and error variance \eqn{\hat\sigma^2}
#' using one of two bootstrap strategies documented in the literature for
#' mixed-effects models.
#'
#' @param Dobj An object of class \code{"Dmethod"} returned by
#'   \code{\link{Dmethod}}. The estimation method (default \code{"VLS"})
#'   is inherited from \code{Dobj}.
#' @param nboot Positive integer. Number of bootstrap samples. Default
#'   \code{200}.
#' @param type Character. Bootstrap strategy:
#'   \describe{
#'     \item{\code{"case"}}{Resample subjects with replacement (default).
#'       Valid under minimal assumptions. Recommended for general use
#'       Thai et al. (2013).}
#'     \item{\code{"residual"}}{Keep subjects fixed, resample marginal
#'       residuals \eqn{\hat e_i = Y_i - \hat f_i} with replacement.
#'       Assumes i.i.d. residuals.
#'       See Carpenter et al. (2003), Thai et al. (2013).}
#'   }
#' @param seed Optional integer seed for reproducibility. Default \code{NULL}.
#' @param verbose Integer (0, 1, or 2). Default \code{1}.
#'
#' @return A list of class \code{"Dboot"} with components:
#'   \describe{
#'     \item{\code{Beta}}{Matrix with columns \code{Estimate} and \code{SE}
#'       for the fixed effects.}
#'     \item{\code{Dhat}}{Matrix with columns \code{Estimate} and \code{SE}
#'       for each variance/covariance component (flattened via
#'       row:column naming).}
#'     \item{\code{Sigma2}}{Matrix with columns \code{Estimate} and \code{SE}
#'       for the error variance.}
#'     \item{\code{Boot_Beta}}{Matrix of bootstrap fixed-effects estimates
#'       (\code{nboot x length(Beta)}).}
#'     \item{\code{Boot_Dhat}}{Array of bootstrap variance component matrices
#'       (\code{k x k x nboot}).}
#'     \item{\code{Boot_Sigma2}}{Numeric vector of bootstrap error variance
#'       estimates (length \code{nboot}).}
#'     \item{\code{type}}{The bootstrap strategy used.}
#'     \item{\code{nboot}}{Number of bootstrap samples requested.}
#'     \item{\code{nfail}}{Number of bootstrap samples that failed.}
#'   }
#'
#' @references
#' Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003).
#' A novel bootstrap procedure for assessing the relationship between
#' class size and achievement.
#' \emph{Journal of the Royal Statistical Society C},
#' \strong{52}, 431--443.
#'
#' Thai, H.T., Mentre, F., Holford, N.H.G., Veyrat-Follet, C. and
#' Comets, E. (2013).
#' A comparison of bootstrap approaches for estimating uncertainty of
#' parameters in linear mixed-effects models.
#' \emph{Pharmaceutical Statistics}, \strong{12}, 129--140.
#'
#' @examples
#' d      <- as.data.frame(Theoph)
#' Expr   <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
#'             (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
#'             (exp(ai2) - exp(ai3))
#' start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
#' random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")
#' DVLS   <- Dmethod(d, Expr, group = "Subject",
#'                   random = random, start = start, method = "VLS")
#' DVLS[c("Dhat", "Sigma2", "Beta")]
#' \donttest{
#' BootSE <- bootstrap_se(DVLS, nboot = 20, type = "case", seed = NULL, verbose = 1)
#' }
#'
#' @importFrom stats sd
#' @export
bootstrap_se <- function(Dobj, nboot = 200,
                         type  = c("case", "residual"),
                         seed  = NULL, verbose = 1) {

  stopifnot(inherits(Dobj, "Dmethod"))
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  #Original data with no NA added
  Dobj   <- c(Dobj, attr(Dobj, "internal")) #internal + public objects
  data   <- Dobj$data
  group  <- Dobj$group
  Ypred_list <- Dobj$Ypred
  resid_list <- Dobj$residuals


  ## --- Extract stored objects ------------------------------------------

  Expr     <- Dobj$Expr
  random   <- Dobj$random
  start    <- Dobj$start
  method   <- Dobj$method
  k        <- nrow(Dobj$Dhat)
  p        <- length(Dobj$Beta)
  response <- as.character(Expr[[2]])

  ids <- unique(data[[group]])
  N   <- length(ids)




  ## For residual pool marginal residuals
  if (type %in% c("residual")) {
    all_resids <- unlist(lapply(ids, function(id) {
      ei <- resid_list[[id]]
      ei[!is.na(ei)]
    }))
  }

  ## --- Containers ------------------------------------------------------
  Boot_Beta   <- matrix(NA, nboot, p, dimnames = list(NULL, names(Dobj$Beta)))
  Boot_Dhat   <- array(NA, dim = c(k, k, nboot))
  Boot_Sigma2 <- numeric(nboot)
  nfail       <- 0L

  .vcat(verbose, 1, "\nBootstrap_se: running ", nboot, " bootstrap samples (type = '", type, "') ...")

  for (b in seq_len(nboot)) {
    .vperm(verbose, b, nboot, freq = 10)

    ## --- Build bootstrap dataset per strategy --------------------------
    data_b <- switch(type,

                     ## Case bootstrap: resample subjects with replacement
                     case = {
                       ids_b <- sample(ids, N, replace = TRUE)
                       do.call(rbind, lapply(seq_along(ids_b), function(j) {
                         di          <- data[data[[group]] == ids_b[j], , drop = FALSE]
                         di[[group]] <- j
                         di
                       }))
                     },

                     ## Residual bootstrap: resample marginal residuals
                     residual = {
                       do.call(rbind, lapply(ids, function(id) {
                         di     <- data[data[[group]] == id, , drop = FALSE]
                         fhat   <- Ypred_list[[id]]
                         obs    <- !is.na(resid_list[[id]])
                         ni_obs <- sum(obs)
                         ei_star <- sample(all_resids, ni_obs, replace = TRUE)
                         di[[response]][obs] <- fhat[obs] + ei_star
                         di
                       }))
                     }
    )

    ## --- Refit Dmethod on bootstrap sample -----------------------------
    Db <- tryCatch(
      suppressWarnings(
        Dmethod(data_b, Expr, group, random, start,
                method       = method,
                is_permuting = FALSE,
                verbose      = 0)
      ),
      error = function(e) NULL
    )

    if (is.null(Db)) {
      nfail <- nfail + 1L
      next
    }

    Boot_Beta[b, ]   <- Db$Beta
    Boot_Dhat[, , b] <- Db$Dhat
    Boot_Sigma2[b]   <- Db$Sigma2
  }

  ## --- Remove failed samples -------------------------------------------
  ok          <- complete.cases(Boot_Beta) & !is.na(Boot_Sigma2)
  Boot_Beta   <- Boot_Beta[ok, , drop = FALSE]
  Boot_Dhat   <- Boot_Dhat[, , ok, drop = FALSE]
  Boot_Sigma2 <- Boot_Sigma2[ok]

  .vcat(verbose, 1, "\n  Completed : ", sum(ok), " / ", nboot, " samples (", nfail, " failed)\n")
  if (sum(ok) < 2) {stop("Fewer than 2 successful bootstrap samples; cannot compute SEs.")}

  ## --- Compute SEs -----------------------------------------------------
  SE_Beta   <- apply(Boot_Beta,   2,      sd)
  SE_Dhat   <- apply(Boot_Dhat,   c(1,2), sd)
  SE_Sigma2 <- sd(Boot_Sigma2)
  dimnames(SE_Dhat) <- dimnames(Dobj$Dhat)

  #Final
  Beta_table   <- cbind(  Estimate = Dobj$Beta,        SE = SE_Beta )
  Sigma2_table <- cbind( Estimate = Dobj$Sigma2,     SE = SE_Sigma2 )
  Dhat_table   <- cbind(Estimate = .unmatrix(Dobj$Dhat),SE = .unmatrix(SE_Dhat))


  structure(
    list(
      Beta     = Beta_table,
      Dhat     = Dhat_table,
      Sigma2   = Sigma2_table,
      Boot_Beta   = Boot_Beta,
      Boot_Dhat   = Boot_Dhat,
      Boot_Sigma2 = Boot_Sigma2,
      type        = type,
      nboot       = nboot,
      nfail       = nfail
    ),
    class = "Dboot"
  )
}

# ============================================================
#  to unmatrix
# ============================================================
#' @keywords internal
.unmatrix <- function(x, byrow = FALSE, sep = ":") {
  stopifnot(is.matrix(x))

  rn <- rownames(x)
  cn <- colnames(x)

  if (is.null(rn)) rn <- seq_len(nrow(x))
  if (is.null(cn)) cn <- seq_len(ncol(x))

  if (byrow) {
    vals <- as.vector(t(x))
    nms <- as.vector(t(outer(rn, cn, paste, sep = sep)))
  } else {
    vals <- as.vector(x)
    nms <- as.vector(outer(rn, cn, paste, sep = sep))
  }

  stats::setNames(vals, nms)
}


# ============================================================
# Print and summary methods
# ============================================================
#' @export
print.Dboot <- function(x, ...) {
  cat("\nTestREnlme: Bootstrap standard errors\n")
  cat("  Type              :", x$type, "\n")
  cat("  Bootstrap samples :", x$nboot,
      "(", x$nfail, "failed)\n\n")
  cat("Fixed effects SE:\n")
  print(round(x$Beta, 6))
  cat("\nError variance SE:\n")
  cat(" ", round(x$Sigma2, 6), "\n")
  cat("\nVariance components SE:\n")
  print(round(x$Dhat, 6))
  invisible(x)
}

#' @export
summary.Dboot <- function(object, ...) print(object, ...)
