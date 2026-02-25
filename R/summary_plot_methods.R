#### summary.R - Summary and Plot Methods for bayesmulticat ####

#' Summary Method for bmvb Objects
#'
#' Provides detailed summary statistics for bmvb analysis results.
#'
#' @param object A bmvb object from bmvb() function.
#' @param prob Numeric. Probability for credible intervals. Default is 0.95.
#' @param ... Additional arguments (not used).
#'
#' @return A list of class "summary.bmvb" with detailed statistics.
#' @export
summary.bmvb <- function(object, prob = 0.95, ...) {

  # Create summary structure
  summary_obj <- list(
    call = object$call,
    estimates = object$estimates,
    sample_sizes = object$sample_sizes,
    delta = object$delta,
    info = object$info
  )

  # Add credible intervals if samples available
  if (!is.null(object$samples)) {
    alpha <- 1 - prob

    # CIs for group A
    ci_a <- apply(object$samples$theta_a, 2, function(x) {
      quantile(x, probs = c(alpha/2, 1 - alpha/2))
    })

    # CIs for group B
    ci_b <- apply(object$samples$theta_b, 2, function(x) {
      quantile(x, probs = c(alpha/2, 1 - alpha/2))
    })

    # CIs for delta
    ci_delta <- apply(object$samples$delta, 2, function(x) {
      quantile(x, probs = c(alpha/2, 1 - alpha/2))
    })

    summary_obj$credible_intervals <- list(
      prob = prob,
      theta_a = ci_a,
      theta_b = ci_b,
      delta = ci_delta
    )

    # Effective sample sizes
    if (requireNamespace("coda", quietly = TRUE)) {
      ess_a <- coda::effectiveSize(object$samples$theta_a)
      ess_b <- coda::effectiveSize(object$samples$theta_b)
      ess_delta <- coda::effectiveSize(object$samples$delta)

      summary_obj$effective_n <- list(
        theta_a = ess_a,
        theta_b = ess_b,
        delta = ess_delta
      )
    }
  }

  class(summary_obj) <- "summary.bmvb"
  return(summary_obj)
}

#' Print Method for summary.bmvb Objects
#'
#' @param x A summary.bmvb object.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#' @export
print.summary.bmvb <- function(x, digits = 3, ...) {

  cat("\nMultivariate Bernoulli Analysis - Summary\n")
  cat("===========================================\n\n")

  cat("Group Estimates:\n")
  est_tab <- rbind.data.frame(
    c(x$info$grp_a, x$estimates$mean_a, x$estimates$sd_a),
    c(x$info$grp_b, x$estimates$mean_b, x$estimates$sd_b)
  )
  colnames(est_tab) <- c("Group", "mean y1", "mean y2", "sd y1", "sd y2")
  print(est_tab, row.names = FALSE, digits = digits)
  cat("\n")

  # Credible intervals if available
  if (!is.null(x$credible_intervals)) {
    cat(sprintf("%d%% Credible Intervals:\n", x$credible_intervals$prob * 100))
    cat(sprintf("  %s, y1: [%.3f, %.3f]\n", x$info$grp_a,
                x$credible_intervals$theta_a[1, 1],
                x$credible_intervals$theta_a[2, 1]))
    cat(sprintf("  %s, y2: [%.3f, %.3f]\n", x$info$grp_a,
                x$credible_intervals$theta_a[1, 2],
                x$credible_intervals$theta_a[2, 2]))
    cat(sprintf("  %s, y1: [%.3f, %.3f]\n", x$info$grp_b,
                x$credible_intervals$theta_b[1, 1],
                x$credible_intervals$theta_b[2, 1]))
    cat(sprintf("  %s, y2: [%.3f, %.3f]\n\n", x$info$grp_b,
                x$credible_intervals$theta_b[1, 2],
                x$credible_intervals$theta_b[2, 2]))
  }

  cat("Sample Sizes:\n")
  cat(sprintf("  %s: n = %d\n", x$info$grp_a, x$sample_sizes$n_a))
  cat(sprintf("  %s: n = %d\n\n", x$info$grp_b, x$sample_sizes$n_b))

  cat("Treatment Effect:\n")
  cat(sprintf("  Delta mean: %s\n", paste(round(x$delta$mean_delta, digits), collapse = ", ")))
  cat(sprintf("  Delta SE:   %s\n", paste(round(x$delta$se_delta, digits), collapse = ", ")))

  if (!is.null(x$credible_intervals)) {
    cat(sprintf("  Delta %d%% CI (y1): [%.3f, %.3f]\n",
                x$credible_intervals$prob * 100,
                x$credible_intervals$delta[1, 1],
                x$credible_intervals$delta[2, 1]))
    cat(sprintf("  Delta %d%% CI (y2): [%.3f, %.3f]\n",
                x$credible_intervals$prob * 100,
                x$credible_intervals$delta[1, 2],
                x$credible_intervals$delta[2, 2]))
  }

  if (x$info$rule == "Comp" && !is.null(x$delta$w_delta)) {
    cat(sprintf("  Weighted delta: %.3f\n", x$delta$w_delta))
  }

  cat(sprintf("  Posterior probability: %.3f\n\n", x$delta$pop))

  cat("Test Information:\n")
  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(round(x$info$w, 3), collapse = ", ")))
  }
  cat(sprintf("  Hypothesis: %s\n", x$info$test_label))

  # Effective sample sizes if available
  if (!is.null(x$effective_n)) {
    cat("\nMCMC Diagnostics:\n")
    cat("  Effective sample sizes:\n")
    cat(sprintf("    theta_%s: %s\n", x$info$grp_a,
                paste(round(x$effective_n$theta_a, 0), collapse = ", ")))
    cat(sprintf("    theta_%s: %s\n", x$info$grp_b,
                paste(round(x$effective_n$theta_b, 0), collapse = ", ")))
    cat(sprintf("    delta: %s\n",
                paste(round(x$effective_n$delta, 0), collapse = ", ")))
  }

  invisible(x)
}

#' Summary Method for bglm Objects
#'
#' @param object A bglm object.
#' @param prob Numeric. Probability for credible intervals. Default is 0.95.
#' @param ... Additional arguments.
#'
#' @return A summary.bglm object.
#' @export
summary.bglm <- function(object, prob = 0.95, ...) {

  summary_obj <- summary.bmvb(object, prob = prob, ...)

  # Add regression coefficient information
  summary_obj$regression <- list(
    coefficients = object$estimates$b,
    se = object$estimates$b_se,
    marg_method = object$info$marg_method,
    sub_pop = object$info$sub_pop
  )

  class(summary_obj) <- c("summary.bglm", "summary.bmvb")
  return(summary_obj)
}

#' Print Method for summary.bglm Objects
#'
#' @param x A summary.bglm object.
#' @param digits Number of digits. Default is 3.
#' @param ... Additional arguments.
#'
#' @export
print.summary.bglm <- function(x, digits = 3, ...) {

  cat("\nBayesian Multivariate Logistic Regression - Summary\n")
  cat("====================================================\n\n")

  # Call parent print method for common elements
  NextMethod("print")

  # Add regression-specific information
  cat("\nRegression Coefficients:\n")
  cat(sprintf("  Marginalization: %s\n", x$regression$marg_method))
  cat(sprintf("  Subpopulation: [%s]\n",
              paste(x$regression$sub_pop, collapse = ", ")))

  invisible(x)
}


#' Summary Method for bglmm Objects
#'
#' @param object A bglmm object.
#' @param prob Numeric. Probability for credible intervals.
#' @param ... Additional arguments.
#'
#' @export
summary.bglmm <- function(object, prob = 0.95, ...) {

  summary_obj <- summary.bglm(object, prob = prob, ...)

  # Add multilevel-specific information
  summary_obj$multilevel <- list(
    n_clusters = object$sample_sizes$J,
    diagnostics = object$diags
  )

  class(summary_obj) <- c("summary.bglmm", "summary.bglm", "summary.bmvb")
  return(summary_obj)
}

#' Print Method for summary.bglmm Objects
#'
#' @param x A summary.bglmm object.
#' @param digits Number of digits.
#' @param ... Additional arguments.
#'
#' @export
print.summary.bglmm <- function(x, digits = 3, ...) {

  cat("\nBayesian Multilevel Multivariate Logistic Regression - Summary\n")
  cat("===============================================================\n\n")

  # Print most elements from parent
  NextMethod("print")

  # Add multilevel-specific information
  cat("\nMultilevel Structure:\n")
  cat(sprintf("  Number of clusters: %d\n", x$multilevel$n_clusters))

  if (!is.null(x$multilevel$diagnostics)) {
    cat("\nMCMC Convergence:\n")
    if (!is.null(x$multilevel$diagnostics$g)) {
      cat(sprintf("  Random effects MPSRF: %.4f\n",
                  x$multilevel$diagnostics$g$convergence$mpsrf))
    }
    if (!is.null(x$multilevel$diagnostics$tau)) {
      cat(sprintf("  Variance components MPSRF: %.4f\n",
                  x$multilevel$diagnostics$tau$convergence$mpsrf))
    }
  }

  invisible(x)
}
