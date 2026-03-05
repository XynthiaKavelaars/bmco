#### bmvb.R - Multivariate Bernoulli Functions ####

#' Estimate Theta with Multivariate Bernoulli Distribution
#'
#' Estimate success probabilities with a Multivariate Bernoulli distribution
#' (reference approach).
#'
#' @param Y n x K matrix with bivariate binary responses.
#' @param prior_alpha Numeric vector of length Q with prior hyperparameters,
#'   or scalar if all are equal. Default is 0.5.
#' @param n_it Scalar. Number of draws from posterior distributions. Default is 10000.
#'
#' @return nIt x K matrix with bivariate Bernoulli probabilities. Currently
#'   supported for K=2 only.
#' @keywords internal
estimate_theta_mvb <- function(Y, prior_alpha = 0.5, n_it = 1e4) {
  y_mult <- multivariate2multinomial(Y)

  # Set prior
  if (length(prior_alpha) == 1) {
    prior <- rep(prior_alpha, ncol(y_mult))
  } else if (length(prior_alpha) == ncol(y_mult)) {
    prior <- prior_alpha
  } else {
    stop("prior_alpha must have length 1 or length Q")
  }

  # Draw from Dirichlet posterior
  m_phi <- MCMCpack::rdirichlet(n_it, colSums(y_mult) + prior)

  # Transform to marginal probabilities
  m_theta <- cbind(rowSums(m_phi[, c(1, 2)]), rowSums(m_phi[, c(1, 3)]))

  return(m_theta)
}

#' Print Method for bmvb Objects
#'
#' @param x A bmvb object.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' set.seed(2024)
#' trial_data <- data.frame(
#'   treatment = rep(c("placebo", "drug"), each = 50),
#'   y1 = rbinom(100, 1, rep(c(0.40, 0.60), each = 50)),
#'   y2 = rbinom(100, 1, rep(c(0.50, 0.70), each = 50))
#' )
#' fit <- bmvb(
#'   data = trial_data, grp = "treatment",
#'   grp_a = "placebo", grp_b = "drug",
#'   y_vars = c("y1", "y2"), n_it = 1000
#' )
#' print(fit)
#'
#' @seealso \code{\link{summary.bmvb}}
#' @export
print.bmvb <- function(x, digits = 3, ...) {

  cat("\nMultivariate Bernoulli Analysis\n")
  cat("=================================\n\n")

  est_tab <- data.frame(
    Group   = c(x$info$grp_a, x$info$grp_b),
    "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
    "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
    check.names = FALSE
  )
  print(est_tab, row.names = FALSE, digits = digits)
  cat(sprintf("  n(%s) = %d    n(%s) = %d\n\n",
              x$info$grp_a, x$sample_sizes$n_a,
              x$info$grp_b, x$sample_sizes$n_b))

  if (x$info$rule == "Comp" && !is.null(x$delta$w_delta)) {
    cat(sprintf("Posterior probability %s [%s, w = %s]: %.3f\n",
                x$info$test_label, x$info$rule,
                paste(round(x$info$w, digits), collapse = ", "),
                x$delta$pop))
  } else {
    cat(sprintf("Posterior probability %s [%s rule]: %.3f\n",
                x$info$test_label, x$info$rule, x$delta$pop))
  }

  cat("\nUse summary() for credible intervals and ESS.\n\n")

  invisible(x)
}
