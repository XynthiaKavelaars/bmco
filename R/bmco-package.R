#' bmco: Bayesian Analysis for Multivariate Categorical Outcomes
#'
#' Provides Bayesian methods for comparing groups on multiple binary outcomes,
#' including basic tests, regression adjustment, and multilevel models.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{bmvb}}: Bayesian test using multivariate Bernoulli
#'   \item \code{\link{bglm}}: Subgroup analysis using Bayesian logistic regression analysis
#'   \item \code{\link{bglmm}}: Multilevel data using Bayesian multilevel logistic regression analysis
#' }
#'
#' @docType package
#' @name bmco-package
#'
#' @importFrom graphics abline axis barplot boxplot hist legend par points segments
#' @importFrom stats complete.cases cor integrate model.matrix quantile rnorm sd var
#' @importFrom coda autocorr.plot densplot traceplot
#' @references{
#' \insertRef{Kavelaars2020}{bmco}
#'
#' \insertRef{Kavelaars2024}{bmco}
#'
#' \insertRef{Kavelaars2023}{bmco}
#' }
#' @importFrom Rdpack reprompt
"_PACKAGE"
