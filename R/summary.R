#### summary.R - Summary Methods for bmco Package ####
##
## Design philosophy
## -----------------
##   print.*         : one screenful -- group means, sample sizes, posterior
##                     probability, decision rule.
##   summary.*       : full analysis -- group estimates with SD, credible
##                     intervals (when posterior samples are available),
##                     regression coefficient table, prior specification,
##                     marginalization details, MCMC diagnostics.
##

# -- summary.bmvb / print.summary.bmvb ----------------------------------------

#' Summary Method for bmvb Objects
#'
#' Provides a comprehensive summary of a \code{\link{bmvb}} analysis.
#' When the model was fitted with \code{return_samples = TRUE}, credible
#' intervals and effective sample sizes are included.
#'
#' @param object A \code{bmvb} object returned by \code{\link{bmvb}}.
#' @param prob Numeric. Coverage probability for credible intervals.
#'   Default is \code{0.95}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.bmvb}, a list containing all
#'   fields of \code{object} plus:
#' \describe{
#'   \item{credible_intervals}{If posterior samples are available: a list
#'     with \code{prob} and credible interval matrices for \code{theta_a},
#'     \code{theta_b}, and \code{delta}.}
#'   \item{effective_n}{If posterior samples are available: a list with
#'     effective sample sizes for \code{theta_a}, \code{theta_b}, and
#'     \code{delta}.}
#' }
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
#'   y_vars = c("y1", "y2"), n_it = 1000,
#'   return_samples = TRUE
#' )
#' summary(fit)
#'
#' @seealso \code{\link{bmvb}}, \code{\link{print.bmvb}}
#' @export
summary.bmvb <- function(object, prob = 0.95, ...) {

  summary_obj <- list(
    estimates    = object$estimates,
    sample_sizes = object$sample_sizes,
    delta        = object$delta,
    info         = object$info
  )

  if (!is.null(object$samples)) {
    alpha  <- 1 - prob
    probs  <- c(alpha / 2, 1 - alpha / 2)

    summary_obj$credible_intervals <- list(
      prob    = prob,
      theta_a = apply(object$samples$theta_a, 2, quantile, probs = probs),
      theta_b = apply(object$samples$theta_b, 2, quantile, probs = probs),
      delta   = apply(object$samples$delta,   2, quantile, probs = probs)
    )

    summary_obj$effective_n <- list(
      theta_a = coda::effectiveSize(object$samples$theta_a),
      theta_b = coda::effectiveSize(object$samples$theta_b),
      delta   = coda::effectiveSize(object$samples$delta)
    )
  }

  class(summary_obj) <- "summary.bmvb"
  return(summary_obj)
}

#' Print Method for summary.bmvb Objects
#'
#' @param x A \code{summary.bmvb} object returned by \code{\link{summary.bmvb}}.
#' @param digits Number of digits to display. Default is \code{3}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns \code{x}.
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
#'   y_vars = c("y1", "y2"), n_it = 1000,
#'   return_samples = TRUE
#' )
#' print(summary(fit))
#'
#' @seealso \code{\link{summary.bmvb}}, \code{\link{print.bmvb}}
#' @export
print.summary.bmvb <- function(x, digits = 3, ...) {

  has_ci  <- !is.null(x$credible_intervals)
  pct_lbl <- if (has_ci) paste0(x$credible_intervals$prob * 100, "%") else NULL

  cat("\nMultivariate Bernoulli Analysis\n")
  cat("=================================\n\n")

  ## Group estimates ----
  cat("Group Estimates:\n")
  grps <- c(x$info$grp_a, x$info$grp_b)

  if (has_ci) {
    ci_a <- x$credible_intervals$theta_a
    ci_b <- x$credible_intervals$theta_b
    tab <- data.frame(
      Group = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "CI y1"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 1], ci_a[2, 1]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 1], ci_b[2, 1])),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      "CI y2"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 2], ci_a[2, 2]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 2], ci_b[2, 2])),
      check.names = FALSE
    )
    names(tab)[4] <- paste0(pct_lbl, " CI y1")
    names(tab)[7] <- paste0(pct_lbl, " CI y2")
  } else {
    tab <- data.frame(
      Group   = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      check.names = FALSE
    )
  }
  print(tab, row.names = FALSE, digits = digits)
  cat(sprintf("\n  n(%s) = %d    n(%s) = %d\n\n",
              x$info$grp_a, x$sample_sizes$n_a,
              x$info$grp_b, x$sample_sizes$n_b))

  ## Treatment effect ----
  cat("Treatment Effect:\n")
  cat(sprintf("  Delta mean (y1, y2): %s\n",
              paste(round(x$delta$mean_delta, digits), collapse = ", ")))
  cat(sprintf("  Delta SE   (y1, y2): %s\n",
              paste(round(x$delta$se_delta,   digits), collapse = ", ")))

  if (has_ci) {
    ci_d <- x$credible_intervals$delta
    cat(sprintf("  %s CI delta: y1 [%.3f, %.3f]   y2 [%.3f, %.3f]\n",
                pct_lbl,
                ci_d[1, 1], ci_d[2, 1],
                ci_d[1, 2], ci_d[2, 2]))
  }

  if (x$info$rule == "Comp" && !is.null(x$delta$w_delta)) {
    cat(sprintf("  Weighted delta: %.3f\n", x$delta$w_delta))
  }

  cat(sprintf("  Posterior probability %s: %.3f\n\n", x$info$test_label,
              x$delta$pop))

  ## Test information ----
  cat("Test Information:\n")
  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(round(x$info$w, digits), collapse = ", ")))
  }
  cat(sprintf("  Hypothesis: %s\n", x$info$test_label))

  ## Effective sample sizes ----
  if (!is.null(x$effective_n)) {
    cat("\nEffective Sample Sizes:\n")
    cat(sprintf("  theta (%s): %s\n", x$info$grp_a,
                paste(round(x$effective_n$theta_a, 0), collapse = ", ")))
    cat(sprintf("  theta (%s): %s\n", x$info$grp_b,
                paste(round(x$effective_n$theta_b, 0), collapse = ", ")))
    cat(sprintf("  delta:      %s\n",
                paste(round(x$effective_n$delta, 0), collapse = ", ")))
  } else {
    cat("\n  (Run bmvb() with return_samples = TRUE for credible intervals",
        "and ESS.)\n")
  }
  cat("\n")

  invisible(x)
}


# -- summary.bglm / print.summary.bglm ----------------------------------------

#' Summary Method for bglm Objects
#'
#' Provides a comprehensive summary of a \code{\link{bglm}} analysis,
#' including the regression coefficient table, prior specification,
#' MCMC diagnostics (effective sample sizes and \eqn{\hat{R}} per
#' parameter), and, when the model was fitted with
#' \code{return_samples = TRUE}, credible intervals.
#'
#' @param object A \code{bglm} object returned by \code{\link{bglm}}.
#' @param prob Numeric. Coverage probability for credible intervals.
#'   Default is \code{0.95}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.bglm}, a list containing:
#' \describe{
#'   \item{estimates}{Posterior means and SDs of group probabilities and
#'     regression coefficients.}
#'   \item{sample_sizes}{Group sample sizes.}
#'   \item{delta}{Posterior mean differences, SEs, and posterior probability.}
#'   \item{info}{Prior specification, test settings, and marginalization
#'     details.}
#'   \item{credible_intervals}{If posterior samples are available: credible
#'     interval matrices for \code{theta_a}, \code{theta_b}, and
#'     \code{delta}.}
#'   \item{effective_n}{If posterior samples are available: effective sample
#'     sizes for \code{theta_a}, \code{theta_b}, and \code{delta}.}
#'   \item{mcmc_diags}{MCMC diagnostics for the regression coefficients
#'     (effective sample sizes and \eqn{\hat{R}}).}
#' }
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package:
#' summary(bglm_fit)
#'
#' @seealso \code{\link{bglm}}, \code{\link{print.bglm}}
#' @export
summary.bglm <- function(object, prob = 0.95, ...) {

  summary_obj <- list(
    estimates    = object$estimates,
    sample_sizes = object$sample_sizes,
    delta        = object$delta,
    info         = object$info
  )

  ## Credible intervals and ESS from posterior samples (optional) ----
  if (!is.null(object$samples)) {
    alpha <- 1 - prob
    probs <- c(alpha / 2, 1 - alpha / 2)

    summary_obj$credible_intervals <- list(
      prob    = prob,
      theta_a = apply(object$samples$theta_a, 2, quantile, probs = probs),
      theta_b = apply(object$samples$theta_b, 2, quantile, probs = probs),
      delta   = apply(object$samples$delta,   2, quantile, probs = probs)
    )

    summary_obj$effective_n <- list(
      theta_a = coda::effectiveSize(object$samples$theta_a),
      theta_b = coda::effectiveSize(object$samples$theta_b),
      delta   = coda::effectiveSize(object$samples$delta)
    )
  }

  ## MCMC diagnostics from stored diagnostics (always available by default) ----
  if (!is.null(object$diags$b)) {
    summary_obj$mcmc_diags <- list(
      n_eff = object$diags$b$n_eff,
      rhat  = object$diags$b$rhat,
      mpsrf = object$diags$b$convergence$mpsrf
    )
  }

  class(summary_obj) <- "summary.bglm"
  return(summary_obj)
}

#' Print Method for summary.bglm Objects
#'
#' @param x A \code{summary.bglm} object returned by
#'   \code{\link{summary.bglm}}.
#' @param digits Number of digits to display. Default is \code{3}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package:
#' print(summary(bglm_fit))
#'
#' @seealso \code{\link{summary.bglm}}, \code{\link{print.bglm}}
#' @export
print.summary.bglm <- function(x, digits = 3, ...) {

  has_ci  <- !is.null(x$credible_intervals)
  pct_lbl <- if (has_ci) paste0(x$credible_intervals$prob * 100, "%") else NULL

  cat("\nBayesian Multivariate Logistic Regression\n")
  cat("==========================================\n\n")

  ## Group estimates ----
  cat("Group Estimates:\n")
  grps <- c(x$info$grp_a, x$info$grp_b)

  if (has_ci) {
    ci_a <- x$credible_intervals$theta_a
    ci_b <- x$credible_intervals$theta_b
    tab <- data.frame(
      Group = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "CI y1"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 1], ci_a[2, 1]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 1], ci_b[2, 1])),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      "CI y2"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 2], ci_a[2, 2]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 2], ci_b[2, 2])),
      check.names = FALSE
    )
    names(tab)[4] <- paste0(pct_lbl, " CI y1")
    names(tab)[7] <- paste0(pct_lbl, " CI y2")
  } else {
    tab <- data.frame(
      Group   = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      check.names = FALSE
    )
  }
  print(tab, row.names = FALSE, digits = digits)
  cat(sprintf("\n  n(%s) = %d    n(%s) = %d\n\n",
              x$info$grp_a, x$sample_sizes$n_a,
              x$info$grp_b, x$sample_sizes$n_b))

  ## Treatment effect ----
  cat("Treatment Effect:\n")
  cat(sprintf("  Delta mean (y1, y2): %s\n",
              paste(round(x$delta$mean_delta, digits), collapse = ", ")))
  cat(sprintf("  Delta SE   (y1, y2): %s\n",
              paste(round(x$delta$se_delta,   digits), collapse = ", ")))

  if (has_ci) {
    ci_d <- x$credible_intervals$delta
    cat(sprintf("  %s CI delta: y1 [%.3f, %.3f]   y2 [%.3f, %.3f]\n",
                pct_lbl,
                ci_d[1, 1], ci_d[2, 1],
                ci_d[1, 2], ci_d[2, 2]))
  }

  if (x$info$rule == "Comp" && !is.null(x$delta$w_delta)) {
    cat(sprintf("  Weighted delta: %.3f\n", x$delta$w_delta))
  }

  cat(sprintf("  Posterior probability %s: %.3f\n\n", x$info$test_label,
              x$delta$pop))

  ## Regression coefficients ----
  if (!is.null(x$estimates$b)) {
    cat("Regression Coefficients (mean [SD]):\n")
    b    <- x$estimates$b
    b_sd <- x$estimates$b_sd
    Q    <- ncol(b)
    P    <- nrow(b)

    var_names <- c("Intercept",
                   x$info$grp_var,
                   x$info$x_var,
                   paste0(x$info$grp_var, ":", x$info$x_var))

    # Show all categories except the reference (last column, all zeros)
    cat_names <- c("b11", "b10", "b01")
    coef_tab <- as.data.frame(
      matrix("", nrow = P, ncol = length(cat_names),
             dimnames = list(var_names[seq_len(P)], cat_names))
    )
    for (q in seq_along(cat_names)) {
      coef_tab[, q] <- sprintf("%.3f [%.3f]", b[, q], b_sd[, q])
    }
    print(coef_tab)
    cat("\n")
  }

  ## Prior specification ----
  cat("Prior Specification (regression coefficients):\n")
  if (!is.null(x$info$b_mu0)) {
    var_names <- c("Intercept",
                   x$info$grp_var,
                   x$info$x_var,
                   paste0(x$info$grp_var, ":", x$info$x_var))
    df_mu <- setNames(as.data.frame(t(x$info$b_mu0)),
                      var_names[seq_along(x$info$b_mu0)])
    rownames(df_mu) <- NULL
    cat("  Mean:\n")
    print(df_mu)
    cat("  Variance (diagonal of inverse precision):\n")
    df_var <- setNames(
      as.data.frame(t(diag(solve(x$info$b_sigma0)))),
      var_names[seq_along(x$info$b_mu0)]
    )
    rownames(df_var) <- NULL
    print(df_var)
    cat("\n")
  }

  ## Marginalization ----
  cat("Marginalization:\n")
  cat(sprintf("  Method: %s\n", x$info$marg_method))
  cat(sprintf("  (Sub)population: [%s]\n",
              paste(x$info$sub_pop, collapse = ", ")))
  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(round(x$info$w, digits), collapse = ", ")))
  }
  cat("\n")

  ## MCMC diagnostics ----
  if (!is.null(x$mcmc_diags)) {
    cat(sprintf("MCMC Diagnostics (regression coefficients):\n"))
    cat(sprintf("  Multivariate PSRF (MPSRF): %.4f\n", x$mcmc_diags$mpsrf))
    diag_tab <- data.frame(
      Parameter = names(x$mcmc_diags$n_eff),
      ESS       = round(x$mcmc_diags$n_eff, 1),
      Rhat      = round(x$mcmc_diags$rhat, 4),
      row.names = NULL
    )
    print(diag_tab, row.names = FALSE)
    cat("\n")
  } else if (!is.null(x$effective_n)) {
    cat("Effective Sample Sizes (theta / delta):\n")
    cat(sprintf("  theta (%s): %s\n", x$info$grp_a,
                paste(round(x$effective_n$theta_a, 0), collapse = ", ")))
    cat(sprintf("  theta (%s): %s\n", x$info$grp_b,
                paste(round(x$effective_n$theta_b, 0), collapse = ", ")))
    cat(sprintf("  delta:      %s\n",
                paste(round(x$effective_n$delta, 0), collapse = ", ")))
    cat("\n")
  }

  invisible(x)
}


# -- summary.bglmm / print.summary.bglmm --------------------------------------

#' Summary Method for bglmm Objects
#'
#' Provides a comprehensive summary of a \code{\link{bglmm}} analysis,
#' including fixed and random effect tables, variance component estimates,
#' multilevel structure, and MCMC convergence diagnostics.
#'
#' @param object A \code{bglmm} object returned by \code{\link{bglmm}}.
#' @param prob Numeric. Coverage probability for credible intervals.
#'   Default is \code{0.95}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.bglmm}, a list containing:
#' \describe{
#'   \item{estimates}{Posterior means and SDs of group probabilities,
#'     fixed effects, random effects, and variance components.}
#'   \item{sample_sizes}{Group sample sizes and number of clusters.}
#'   \item{delta}{Posterior mean differences, SEs, and posterior probability.}
#'   \item{info}{Prior specification, model structure, test settings, and
#'     marginalization details.}
#'   \item{credible_intervals}{If posterior samples are available: credible
#'     interval matrices for \code{theta_a}, \code{theta_b}, and
#'     \code{delta}.}
#'   \item{effective_n}{If posterior samples are available: effective sample
#'     sizes for \code{theta_a}, \code{theta_b}, and \code{delta}.}
#'   \item{mcmc_diags}{MCMC convergence diagnostics (ESS, \eqn{\hat{R}},
#'     MPSRF) for fixed effects, random effects, and variance components.}
#' }
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package:
#' summary(bglmm_fit)
#'
#' @seealso \code{\link{bglmm}}, \code{\link{print.bglmm}}
#' @export
summary.bglmm <- function(object, prob = 0.95, ...) {

  summary_obj <- list(
    estimates    = object$estimates,
    sample_sizes = object$sample_sizes,
    delta        = object$delta,
    info         = object$info
  )

  ## Credible intervals and ESS from posterior samples (optional) ----
  if (!is.null(object$samples)) {
    alpha <- 1 - prob
    probs <- c(alpha / 2, 1 - alpha / 2)

    summary_obj$credible_intervals <- list(
      prob    = prob,
      theta_a = apply(object$samples$theta_a, 2, quantile, probs = probs),
      theta_b = apply(object$samples$theta_b, 2, quantile, probs = probs),
      delta   = apply(object$samples$delta,   2, quantile, probs = probs)
    )

    summary_obj$effective_n <- list(
      theta_a = coda::effectiveSize(object$samples$theta_a),
      theta_b = coda::effectiveSize(object$samples$theta_b),
      delta   = coda::effectiveSize(object$samples$delta)
    )
  }

  ## MCMC diagnostics from stored diagnostics ----
  mcmc_diags <- list()
  if (!is.null(object$diags$b)) {
    mcmc_diags$b <- list(
      n_eff = object$diags$b$n_eff,
      rhat  = object$diags$b$rhat,
      mpsrf = object$diags$b$convergence$mpsrf
    )
  }
  if (!is.null(object$diags$g)) {
    mcmc_diags$g <- list(
      n_eff = object$diags$g$n_eff,
      rhat  = object$diags$g$rhat,
      mpsrf = object$diags$g$convergence$mpsrf
    )
  }
  if (!is.null(object$diags$tau)) {
    mcmc_diags$tau <- list(
      n_eff = object$diags$tau$n_eff,
      rhat  = object$diags$tau$rhat,
      mpsrf = object$diags$tau$convergence$mpsrf
    )
  }
  if (length(mcmc_diags) > 0) {
    summary_obj$mcmc_diags <- mcmc_diags
  }

  class(summary_obj) <- "summary.bglmm"
  return(summary_obj)
}

#' Print Method for summary.bglmm Objects
#'
#' @param x A \code{summary.bglmm} object returned by
#'   \code{\link{summary.bglmm}}.
#' @param digits Number of digits to display. Default is \code{3}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package:
#' print(summary(bglmm_fit))
#'
#' @seealso \code{\link{summary.bglmm}}, \code{\link{print.bglmm}}
#' @export
print.summary.bglmm <- function(x, digits = 3, ...) {

  has_ci  <- !is.null(x$credible_intervals)
  pct_lbl <- if (has_ci) paste0(x$credible_intervals$prob * 100, "%") else NULL

  cat("\nBayesian Multilevel Multivariate Logistic Regression\n")
  cat("======================================================\n\n")

  ## Multilevel structure ----
  cat(sprintf("Multilevel Structure:  J = %d clusters    ",
              x$sample_sizes$J))
  cat(sprintf("n(%s) = %d    n(%s) = %d\n\n",
              x$info$grp_a, x$sample_sizes$n_a,
              x$info$grp_b, x$sample_sizes$n_b))

  ## Group estimates ----
  cat("Group Estimates:\n")
  grps <- c(x$info$grp_a, x$info$grp_b)

  if (has_ci) {
    ci_a <- x$credible_intervals$theta_a
    ci_b <- x$credible_intervals$theta_b
    tab <- data.frame(
      Group = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "CI y1"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 1], ci_a[2, 1]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 1], ci_b[2, 1])),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      "CI y2"   = c(sprintf("[%.3f, %.3f]", ci_a[1, 2], ci_a[2, 2]),
                    sprintf("[%.3f, %.3f]", ci_b[1, 2], ci_b[2, 2])),
      check.names = FALSE
    )
    names(tab)[4] <- paste0(pct_lbl, " CI y1")
    names(tab)[7] <- paste0(pct_lbl, " CI y2")
  } else {
    tab <- data.frame(
      Group   = grps,
      "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
      "sd y1"   = c(x$estimates$sd_a[1],   x$estimates$sd_b[1]),
      "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
      "sd y2"   = c(x$estimates$sd_a[2],   x$estimates$sd_b[2]),
      check.names = FALSE
    )
  }
  print(tab, row.names = FALSE, digits = digits)
  cat("\n")

  ## Treatment effect ----
  cat("Treatment Effect:\n")
  cat(sprintf("  Delta mean (y1, y2): %s\n",
              paste(round(x$delta$mean_delta, digits), collapse = ", ")))
  cat(sprintf("  Delta SE   (y1, y2): %s\n",
              paste(round(x$delta$se_delta,   digits), collapse = ", ")))

  if (has_ci) {
    ci_d <- x$credible_intervals$delta
    cat(sprintf("  %s CI delta: y1 [%.3f, %.3f]   y2 [%.3f, %.3f]\n",
                pct_lbl,
                ci_d[1, 1], ci_d[2, 1],
                ci_d[1, 2], ci_d[2, 2]))
  }

  if (x$info$rule == "Comp" && !is.null(x$delta$w_delta)) {
    cat(sprintf("  Weighted delta: %.3f\n", x$delta$w_delta))
  }

  cat(sprintf("  Posterior probability %s: %.3f\n\n",
              x$info$test_label, x$delta$pop))

  ## Fixed effects ----
  if (!is.null(x$estimates$b)) {
    cat("Fixed Effects (mean [SD]):\n")
    b    <- x$estimates$b
    b_sd <- x$estimates$b_sd
    Q    <- ncol(b)
    P    <- nrow(b)
    cat_names <- c("b11", "b10", "b01")
    coef_tab <- as.data.frame(
      matrix("", nrow = P, ncol = length(cat_names),
             dimnames = list(x$info$fixed[seq_len(P)], cat_names))
    )
    for (q in seq_along(cat_names)) {
      coef_tab[, q] <- sprintf("%.3f [%.3f]", b[, q], b_sd[, q])
    }
    print(coef_tab)
    cat("\n")
  }

  ## Random effects ----
  if (!is.null(x$estimates$g)) {
    cat("Random Effects (population mean [SD]):\n")
    g    <- x$estimates$g
    g_sd <- x$estimates$g_sd
    Q    <- ncol(g)
    P    <- nrow(g)
    cat_names <- c("g11", "g10", "g01")
    rand_tab <- as.data.frame(
      matrix("", nrow = P, ncol = length(cat_names),
             dimnames = list(x$info$random[seq_len(P)], cat_names))
    )
    for (q in seq_along(cat_names)) {
      rand_tab[, q] <- sprintf("%.3f [%.3f]", g[, q], g_sd[, q])
    }
    print(rand_tab)
    cat("\n")
  }

  ## Variance components ----
  if (!is.null(x$estimates$tau)) {
    cat("Variance Components (posterior mean):\n")
    Q_minus1 <- length(x$estimates$tau)
    cat_labels <- c("y1=1, y2=1 (b11)", "y1=1, y2=0 (b10)", "y1=0, y2=1 (b01)")
    for (q in seq_len(Q_minus1)) {
      cat(sprintf("  %s:\n", cat_labels[q]))
      tau_tab <- round(x$estimates$tau[[q]], digits)
      rownames(tau_tab) <- colnames(tau_tab) <- x$info$random
      print(tau_tab)
    }
    cat("\n")
  }

  ## Priors ----
  cat("Prior Specification:\n")
  if (!is.null(x$info$b_mu0) && length(x$info$fixed) > 0) {
    cat("  Fixed effects -- Normal prior:\n")
    cat("    Mean:     ", paste(round(x$info$b_mu0, digits), collapse = ", "), "\n")
    cat("    Variance: ", paste(round(diag(solve(x$info$b_sigma0)), digits),
                                collapse = ", "), "\n")
  }
  if (!is.null(x$info$g_mu0) && length(x$info$random) > 0) {
    cat("  Random effects -- Normal prior:\n")
    cat("    Mean:     ", paste(round(x$info$g_mu0, digits), collapse = ", "), "\n")
    cat("    Variance: ", paste(round(diag(solve(x$info$g_sigma0)), digits),
                                collapse = ", "), "\n")
    cat(sprintf("  Covariance -- Inverse-Wishart: df = %d\n", x$info$nu0))
  }

  ## Marginalization ----
  cat("\nMarginalization:\n")
  cat(sprintf("  Method: %s    (Sub)population: [%s]\n",
              x$info$marg_method,
              paste(x$info$sub_pop, collapse = ", ")))
  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(round(x$info$w, digits), collapse = ", ")))
  }
  cat("\n")

  ## MCMC diagnostics ----
  if (!is.null(x$mcmc_diags)) {
    cat("MCMC Convergence Diagnostics:\n")

    print_diag_block <- function(label, diag) {
      cat(sprintf("  %s -- MPSRF: %.4f\n", label, diag$mpsrf))
      dt <- data.frame(
        Parameter = names(diag$n_eff),
        ESS       = round(diag$n_eff, 1),
        Rhat      = round(diag$rhat,  4),
        row.names = NULL
      )
      print(dt, row.names = FALSE)
      cat("\n")
    }

    if (!is.null(x$mcmc_diags$b))   print_diag_block("Fixed effects",       x$mcmc_diags$b)
    if (!is.null(x$mcmc_diags$g))   print_diag_block("Random effects",      x$mcmc_diags$g)
    if (!is.null(x$mcmc_diags$tau)) print_diag_block("Variance components", x$mcmc_diags$tau)
  }

  invisible(x)
}
