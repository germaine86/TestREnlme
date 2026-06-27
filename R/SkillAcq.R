#' Skill acquisition data
#'
#' Response latency data from a study of quantitative skill acquisition
#' on a learning task (Blozis 2004). Log-transformed response latencies
#' are recorded for \eqn{N = 204} subjects across \eqn{J = 11} trial
#' blocks, stored in wide format. Used to illustrate fitting a nonlinear
#' mixed-effects model with a subject-level covariate (working memory)
#' incorporated through the second-stage model, and to demonstrate
#' reshaping wide-format longitudinal data into the long format required
#' by \code{\link{Dmethod}} and \code{\link{Dhypothesis_test}}.
#'
#' @format A data frame with 204 rows (one per subject) and 13 columns:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{ly1}{Log-transformed response latency at trial block 1.}
#'   \item{ly2}{Log-transformed response latency at trial block 2.}
#'   \item{ly3}{Log-transformed response latency at trial block 3.}
#'   \item{ly4}{Log-transformed response latency at trial block 4.}
#'   \item{ly5}{Log-transformed response latency at trial block 5.}
#'   \item{ly6}{Log-transformed response latency at trial block 6.}
#'   \item{ly7}{Log-transformed response latency at trial block 7.}
#'   \item{ly8}{Log-transformed response latency at trial block 8.}
#'   \item{ly9}{Log-transformed response latency at trial block 9.}
#'   \item{ly10}{Log-transformed response latency at trial block 10.}
#'   \item{ly11}{Log-transformed response latency at trial block 11.}
#'   \item{wm2}{Subject-level working-memory covariate.}
#' }
#'
#' @details
#' \code{SkillAcq} is stored in wide format, as commonly encountered in
#' practice. Before use with \code{\link{Dmethod}} it must be reshaped to
#' long format, with one row per subject-trial observation; see
#' \strong{Examples} below.
#'
#' The nonlinear mixed-effects model used for this dataset is
#' \deqn{Y_{ij} = a_{i1} - (a_{i1} + a_{i0}) \exp(a_{i2} T_{ij}) + \varepsilon_{ij},}
#' with
#' \deqn{a_{ik} = \beta_{0k} + \beta_{1k} wm2_i + b_{ik}, \quad k \in \{0, 1, 2\},}
#' where \eqn{a_{i0}}, \eqn{a_{i1}}, and \eqn{a_{i2}} represent,
#' respectively, the subject-specific initial performance offset, lower
#' asymptote, and learning rate, each with a regression component on
#' \code{wm2} and a subject-specific random effect. The questions of
#' interest are whether all three random effects are necessary, and
#' whether one or more can be removed to obtain a more parsimonious
#' model; see \code{\link{Dhypothesis_test}}.
#'
#' @source Blozis, S. A. (2004). Structured latent curve models for the
#'  study of change in multivariate repeated measures.
#' \emph{Psychological Methods}, \strong{9(3)}, 334--353.
#'  https://doi.org/10.1037/1082-989X.9.3.334
#'
#' Used as a worked example in Uwimpuhwe, G., Drikvandi, R. and Blozis,
#' S. A. (in preparation). TestREnlme: An R Package for Testing Random
#' Effects in Nonlinear Mixed-Effects Models. \emph{Journal of
#' Statistical Software}.
#'
#' @examples
#' \donttest{
#' ## Reshape from wide (ly1..ly11) to long format
#' qrt  <- data.frame(SkillAcq)
#' qrt1 <- reshape2::melt(qrt, id.vars = c("id", "wm2"),
#'                        variable.name = "ly", value.name = "Y")
#' qrt1$t <- as.numeric(sub("ly", "", qrt1$ly))
#'
#' ## Model: Y_ij = ai1 - (ai1 + ai0) * exp(ai2 * t_ij) + e_ij
#' ## with aik = B0k + B1k * wm2_i + bik, k in {0, 1, 2}
#' Expr_learn <- Y ~ ai1 - (ai1 + ai0) * exp(ai2 * t)
#' random_learn <- c("ai0 ~ B00 + B10 * wm2 + bi0",
#'                   "ai1 ~ B01 + B11 * wm2 + bi1",
#'                   "ai2 ~ B02 + B12 * wm2 + bi2")
#'
#'
#' ## Estimate variance components (VLS) and test all three random effects
#' DVLS_learn <- Dmethod(qrt1, Expr_learn, group = "id",
#'                       random = random_learn, start = NULL)
#' DVLS_learn[c("Dhat", "Sigma2", "Beta")]
#'
#' H_learn <- Dhypothesis_test(qrt1, Expr_learn, group = "id",
#'                             random = random_learn, start = NULL,
#'                             Dhatt = DVLS_learn, nperm = 200, seed = 1)
#' H_learn$pvalue
#' }
"SkillAcq"
