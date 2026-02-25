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
#' @export
print.bmvb <- function(x, digits = 3, ...) {
  cat("\nMultivariate Bernoulli Analysis\n")
  cat("=================================\n\n")
  
  cat("Group Estimates:\n")
  est_tab <- rbind.data.frame(
    c(x$info$grp_a, x$estimates$mean_a, x$estimates$sd_a),
    c(x$info$grp_b, x$estimates$mean_b, x$estimates$sd_b)
  )
  colnames(est_tab) <- c("Group", "mean y1", "mean y2", "sd y1", "sd y2")
  print(est_tab, row.names = FALSE)
  cat("\n")
  
  cat("Sample Sizes:\n")
  cat(sprintf("  %s: n = %d\n", x$info$grp_a, x$sample_sizes$n_a))
  cat(sprintf("  %s: n = %d\n\n", x$info$grp_b, x$sample_sizes$n_b))
  
  cat("Treatment Effect:\n")
  cat(sprintf("  Delta mean: %s\n", paste(x$delta$mean_delta, collapse = ", ")))
  cat(sprintf("  Delta SE:   %s\n", paste(x$delta$se_delta, collapse = ", ")))
  
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weighted delta mean: %.3f\n", x$delta$w_delta))
  }
  
  cat(sprintf("  Posterior probability: %.3f\n\n", x$delta$pop))
  
  cat("Test Information:\n")
  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(x$info$w, collapse = ", ")))
  }
  cat(sprintf("  Hypothesis: %s\n", x$info$test_label))
  
  invisible(x)
}
