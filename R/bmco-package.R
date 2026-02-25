#' bmco: Bayesian Analysis for Multivariate Categorical Outcomes
#'
#' Provides Bayesian methods for comparing groups on multiple binary outcomes,
#' including basic tests, regression adjustment, and multilevel models.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{bmvb}}: Basic Bayesian test using multivariate Bernoulli
#'   \item \code{\link{bglm}}: Test with covariate adjustment
#'   \item \code{\link{bglmm}}: Multilevel/clustered data analysis
#' }
#'
#' @docType package
#' @name bmco-package
#'
#' @importFrom graphics abline axis barplot boxplot hist legend par points segments
#' @importFrom stats complete.cases cor integrate model.matrix quantile rnorm sd var
#' @importFrom coda autocorr.plot densplot traceplot
"_PACKAGE"
