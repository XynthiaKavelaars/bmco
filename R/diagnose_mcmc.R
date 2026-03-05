#### diagnose_mcmc.R - MCMC Diagnostics ####

#' Diagnose MCMC Chains
#'
#' Perform diagnostic checks on MCMC chains including convergence assessment
#' and optional trace plots.
#'
#' @param chains An mcmc.list object or a list of mcmc.list objects containing
#'   the MCMC chains to diagnose.
#' @param param_names Optional character vector of parameter names for plot labels.
#'   If NULL, uses column names from chains.
#'
#' @return A list with diagnostic information:
#' \describe{
#'   \item{convergence}{Multivariate potential scale reduction factor (Gelman-Rubin statistic)}
#'   \item{n_eff}{Effective sample sizes for each parameter}
#'   \item{rhat}{Univariate potential scale reduction factors for each parameter}
#'   \item{summary}{Posterior summary of MCMC chains}
#'
#' }
#' @keywords internal
#' @importFrom coda gelman.diag effectiveSize

diagnose_mcmc <- function(chains,  param_names = NULL) {

  # Input validation
  if (!inherits(chains, "mcmc.list")) {
    stop("chains must be an mcmc.list object")
  }

  # Convergence diagnostics
  gelman_result <- coda::gelman.diag(chains, multivariate = TRUE)
  convergence <- list(
    mpsrf = gelman_result$mpsrf,
    psrf = gelman_result$psrf
  )

  # Effective sample sizes
  n_eff <- coda::effectiveSize(chains)

   # Specific warnings for different convergence issues
  if (any(n_eff < 100)) {
    warning("Low effective sample size for some parameters. ",
            "Results may be unreliable.")
  }

  # Get parameter names
  if (is.null(param_names)) {
    param_names <- colnames(chains[[1]])
  }

  # Prepare output
  diagnostics <- list(
    convergence = convergence,
    n_eff = n_eff,
    rhat = convergence$psrf[, 1],
    summary = summary(chains)
  )

  class(diagnostics) <- "mcmc_diagnostics"
  return(diagnostics)
}

#' Print Method for mcmc_diagnostics Objects
#'
#' @param x A mcmc_diagnostics object.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#' @keywords internal
print.mcmc_diagnostics <- function(x, ...) {
  cat("\nMCMC Diagnostics\n")
  cat("==================\n\n")

  cat("Convergence (Gelman-Rubin):\n")
  cat(sprintf("  Multivariate PSRF: %.4f\n", x$convergence$mpsrf))
  cat("\n  Univariate PSRF:\n")
  print(round(x$convergence$psrf, 4))
  cat("\n")

  cat("Effective Sample Sizes:\n")
  print(round(x$n_eff, 1))
  cat("\n")

  # Check for convergence issues
  if (x$convergence$mpsrf > 1.1) {
    warning("Multivariate PSRF > 1.1 suggests chains have not converged. Consider running more iterations.")
  }

  if (any(x$rhat > 1.1)) {
    warning("Some parameters have Rhat > 1.1, suggesting lack of convergence.")
  }

  invisible(x)
}

#' Plot Method for bglm Objects
#'
#' @param x A bglm object returned by bglm().
#' @param type Character. Type of plot: "trace", "density", "autocorr", or "all". Default is "all".
#' @param parameters Character vector. Which parameters to plot. Default is NULL (all parameters).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns NULL (plots are displayed).
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package.
#' # Plot trace plots for the fixed-effect regression coefficients:
#' plot(bglm_fit, type = "trace")
#'
#' @export
#' @importFrom coda traceplot densplot autocorr.plot
plot.bglm <- function(x, type = "all", parameters = NULL, ...) {

  # Get chains — stored in samples$b when return_samples=TRUE or return_diagnostic_plots=TRUE
  if (!is.null(x$samples) && !is.null(x$samples$b)) {
    chains <- x$samples$b
  } else {
    stop("No MCMC chains available for plotting. Run bglm() with return_samples = TRUE or return_diagnostic_plots = TRUE.")
  }

  # Subset parameters if requested
  if (!is.null(parameters)) {
    chains <- chains[, parameters, drop = FALSE]
  }

  # Determine plot types
  if (type == "all") {
    plot_types <- c("trace", "density", "autocorr")
  } else {
    plot_types <- type
  }

  # Generate plots
  if ("trace" %in% plot_types) {
    coda::traceplot(chains, main = "Trace Plots - Fixed Effects", ...)
  }

  if ("density" %in% plot_types) {
    coda::densplot(chains, main = "Posterior Densities - Fixed Effects", ...)
  }

  if ("autocorr" %in% plot_types) {
    coda::autocorr.plot(chains, main = "Autocorrelation - Fixed Effects", ...)
  }

  invisible(NULL)
}

#' Plot Method for bglmm Objects
#'
#' @param x A bglmm object returned by bglmm().
#' @param type Character. Type of plot: "trace", "density", "autocorr", or "all". Default is "all".
#' @param which Character. Which component to plot: "fixed" (b), "random" (g), "variance" (tau), or "all". Default is "fixed".
#' @param parameters Character vector. Which parameters to plot. Default is NULL (all parameters).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns NULL (plots are displayed).
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package.
#' # Trace plots for the fixed-effect regression coefficients:
#' plot(bglmm_fit, type = "trace", which = "fixed")
#'
#' # Trace plots for the random-effect variance components:
#' plot(bglmm_fit, type = "trace", which = "variance")
#'
#' @export
#' @importFrom coda traceplot densplot autocorr.plot
plot.bglmm <- function(x, type = "all", which = "fixed", parameters = NULL, ...) {

  # Check that samples exist (chains are stored in samples$b/g/tau)
  if (is.null(x$samples)) {
    stop("No MCMC chains available for plotting. Run bglmm() with return_samples = TRUE or return_diagnostic_plots = TRUE.")
  }

  # Determine which component to plot
  components <- if (which == "all") {
    c("fixed", "random", "variance")
  } else {
    which
  }

  # Determine plot types
  if (type == "all") {
    plot_types <- c("trace", "density", "autocorr")
  } else {
    plot_types <- type
  }

  # Plot each component
  for (comp in components) {

    # Get appropriate chains from samples (stored as mcmc.list)
    chains <- switch(comp,
                     "fixed"    = x$samples$b,
                     "random"   = x$samples$g,
                     "variance" = x$samples$tau
    )

    if (is.null(chains)) {
      message("No chains available for component: ", comp)
      next
    }

    # Subset parameters if requested
    if (!is.null(parameters)) {
      chains <- chains[, parameters, drop = FALSE]
    }

    # Title prefix
    title_prefix <- switch(comp,
                           "fixed" = "Fixed Effects",
                           "random" = "Random Effects",
                           "variance" = "Variance Components"
    )

    # Generate plots
    if ("trace" %in% plot_types) {
      coda::traceplot(chains, main = paste("Trace Plots -", title_prefix), ...)
    }

    if ("density" %in% plot_types) {
      coda::densplot(chains, main = paste("Posterior Densities -", title_prefix), ...)
    }

    if ("autocorr" %in% plot_types) {
      coda::autocorr.plot(chains, main = paste("Autocorrelation -", title_prefix), ...)
    }
  }

  invisible(NULL)
}

#' Extract MCMC Chains from Model
#'
#' Extract MCMC chains from a fitted bglm or bglmm model for custom plotting or analysis.
#'
#' @param model A bglm or bglmm object.
#' @param component Character. For bglmm: "fixed" (b), "random" (g), or "variance" (tau).
#'   For bglm: only "fixed" is available. Default is "fixed".
#'
#' @return An mcmc.list object containing the MCMC chains.
#' @keywords internal
extract_chains <- function(model, component = "fixed") {

  if (is.null(model$samples)) {
    stop("No MCMC chains available. Run model with return_samples = TRUE or return_diagnostic_plots = TRUE.")
  }

  # For bglm
  if (inherits(model, "bglm")) {
    if (component != "fixed") {
      warning("bglm only has fixed effects. Returning fixed effect chains.")
    }
    if (!is.null(model$samples$b)) {
      return(model$samples$b)
    } else {
      stop("No chains available. Run bglm() with return_samples = TRUE or return_diagnostic_plots = TRUE.")
    }
  }

  # For bglmm
  if (inherits(model, "bglmm")) {
    chains <- switch(component,
                     "fixed"    = model$samples$b,
                     "random"   = model$samples$g,
                     "variance" = model$samples$tau,
                     stop("Invalid component. Choose 'fixed', 'random', or 'variance'.")
    )

    if (is.null(chains)) {
      stop("No chains available for component: ", component,
           ". Run bglmm() with return_samples = TRUE or return_diagnostic_plots = TRUE.")
    }

    return(chains)
  }

  stop("Model must be a bglm or bglmm object.")
}

#' Plot Specific MCMC Diagnostic
#'
#' Generate a specific type of diagnostic plot for MCMC chains.
#'
#' @param chains An mcmc.list object or a fitted model.
#' @param type Character. Type of plot: "trace", "density", or "autocorr".
#' @param parameters Character vector. Which parameters to plot. NULL plots all.
#' @param main Character. Main title for the plot.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @importFrom coda traceplot densplot autocorr.plot
plot_mcmc <- function(chains, type = "trace", parameters = NULL, main = NULL, ...) {

  # If chains is a model object, extract chains
  if (inherits(chains, "bglm") || inherits(chains, "bglmm")) {
    chains <- extract_chains(chains)
  }

  if (!inherits(chains, "mcmc.list")) {
    stop("chains must be an mcmc.list object or a fitted model.")
  }

  # Subset parameters if requested
  if (!is.null(parameters)) {
    chains <- chains[, parameters, drop = FALSE]
  }

  # Generate plot
  switch(type,
    "trace" = {
      if (is.null(main)) main <- "MCMC Trace Plot"
      coda::traceplot(chains, main = main, ...)
    },
    "density" = {
      if (is.null(main)) main <- "Posterior Density"
      coda::densplot(chains, main = main, ...)
    },
    "autocorr" = {
      if (is.null(main)) main <- "Autocorrelation"
      coda::autocorr.plot(chains, main = main, ...)
    },
    stop("Invalid type. Choose 'trace', 'density', or 'autocorr'.")
  )

  invisible(NULL)
}
