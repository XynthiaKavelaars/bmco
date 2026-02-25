#### bglm.R - Bayesian Generalized Linear Model Functions ####

#' Estimate Parameters using Polya-Gamma Method
#'
#' Wrapper function for MCMC procedure to estimate regression coefficients using
#' Polya-Gamma Gibbs sampling.
#'
#' @param X (n x P) design matrix.
#' @param Y (n x Q) matrix of multinomial response data.
#' @param n_burn Scalar. Number of burnin iterations.
#' @param n_it Scalar. Number of iterations.
#' @param start Vector of two starting values, each used for a different chain.
#' @param b_mu0 Scalar or vector (length(P)) of prior mean(s).
#' @param b_sigma0 Scalar or matrix (P x P) of prior variance of regression coefficients.
#' @param n_chain Scalar. Number of chains to be sampled.
#'
#' @return A list of length n_chain, being a list of length n_it with P x Q matrices with estimated regression coefficients
#' @keywords internal
estimate_parameters_pg <- function(X, Y, n_burn, n_it, start, b_mu0, b_sigma0, n_chain) {
  tryCatch(
    {
      P <- ncol(X)
      Q <- ncol(Y)
  chains_pg <- vector("list", n_chain)

      for(chain in 1:n_chain){
      chains_pg[[chain]] <- sample_beta_pg(
        X = X, Y = Y, n_burn = n_burn, n_it = n_it,
        start = start[chain], b_mu0 = b_mu0, b_sigma0 = b_sigma0
      )

#      chains_pg[[chain]] <- coda::mcmc(do.call(rbind, lapply(chain_pg[["b_draw_pg"]], function(i) {
#          as.vector(i[, -seq(Q, Q * P, Q)])
#        })))
      }

      res <- chains_pg#coda::as.mcmc.list(chains_pg)
    },
    warning = function(warn) {
      message(paste("Warning: ", warn))
      res <- list(NULL)
    },
    error = function(err) {
      message(paste("Error: ", err))
      res <- list(NULL)
    },
    finally = function(f) {
      return(res)
    }
  )
}

#' Sample Beta Coefficients using Polya-Gamma Method
#'
#' Sample regression coefficients via Polya-Gamma multinomial logistic regression
#' for a single chain.
#'
#' @param X (n x P) design matrix.
#' @param Y (n x Q) matrix of multinomial response data.
#' @param n_burn Scalar. Number of burnin iterations.
#' @param n_it Scalar. Number of iterations.
#' @param start Scalar. Starting value for this chain.
#' @param b_mu0 Scalar. Prior mean.
#' @param b_sigma0 Scalar. Prior precision of regression coefficients.
#' @param verbose Logical. If true, a progress bar is shown.
#'
#' @return b_draw_pg: List of length n_it with P x Q matrices with estimated regression coefficients
#' @keywords internal
sample_beta_pg <- function(X, Y, n_burn, n_it, start, b_mu0, b_sigma0, verbose = FALSE) {

  P <- ncol(X)
  n <- nrow(Y)
  Q <- ncol(Y)

  b_draw_pg <- list("vector", n_it)
  beta_pg <- array(start, dim = c(P, Q))
  beta_pg[, Q] <- rep(0, P)

  # kappa: Note - As opposed to the paper, kappa here is not divided by omega
  # to improve programming efficiency
  kappa <- Y - 1/2

  # Prior setup
  P0 <- b_sigma0 %*% b_mu0

  if(verbose){
  # Create progress bar
  pb <- tcltk::tkProgressBar(
    title = "MCMC Progress",
    min = 1,
    max = n_burn + n_it,
    width = 300
  )
  }
  # Start Gibbs sampler
  for (i in 2:(n_burn + n_it)) {
    for (q in 1:(Q - 1)) {
      # Compute linear predictor
      C <- log(rowSums(exp(X %*% beta_pg[, -q])))
      eta <- X %*% beta_pg[, q] - C

      # Draw Polya-Gamma variable
      omega_draw <- pgdraw::pgdraw(1, eta)

      # Draw regression coefficients
      b_sigma <- chol2inv(chol(t(X) %*% (X * omega_draw) + b_sigma0))
      b_mu <- b_sigma %*% (t(X) %*% (kappa[, q] + omega_draw * C) + P0)
      beta_pg[, q] <- b_mu + t(chol(b_sigma)) %*% rnorm(P)
    }
    if (i > n_burn) {
      b_draw_pg[[i - n_burn]] <- beta_pg
    }

    if(verbose){
  tcltk::setTkProgressBar(
    pb, i,
    label = paste("Iteration", i, "of", n_burn + n_it)
  )
  }
}

if(verbose){
close(pb)
}
  return(b_draw_pg)
}

#' Transform Regression Coefficients to Success Probabilities
#'
#' Compute theta (success probabilities) from regression coefficients via various methods.
#'
#' @param beta_draw_pg List of n_it (P x Q) matrices with posterior regression coefficients.
#' @param X (n x P) matrix with covariate data.
#' @param Y (n x Q) matrix with multinomial response data. Default is NULL.
#' @param measurement_level Character. "discrete" or "continuous".
#' @param method "Value" for vector of fixed values, "Empirical" for empirical
#'   marginalization, "Analytical" for numerical marginalization.
#' @param range If method = "Analytical" or "Empirical" and if measurement_level is
#'   continuous: range that defines the population of interest by a vector containing
#'   a lower and an upper bound.
#' @param value If method = "Value": value that defines the population of interest
#'   by a scalar value.
#'
#' @return A list with:
#' \describe{
#'   \item{m_theta_a}{List of n_it vectors of length K with multivariate probabilities for group A}
#'   \item{m_theta_b}{List of n_it vectors of length K with multivariate probabilities for group B}
#' }
#' @keywords internal
transform2theta <- function(beta_draw_pg, X, Y = NULL, grp_var, population_var, measurement_level,
                           method, range, value) {

  # Validate inputs
  if (method %in% c("Empirical", "Analytical")) {
    if (measurement_level == "discrete") {
      stop("Empirical or analytical methods not applicable to discrete covariate. ",
           "Choose method = 'Value' and specify the value parameter.")
    }

    if (is.numeric(range) && length(range) == 2 &&
        (
          (is.finite(range[1]) && is.infinite(range[2]) && range[2] > 0) ||  # (number, Inf)
          (is.infinite(range[1]) && range[1] < 0 && is.finite(range[2])) ||  # (-Inf, number)
          (is.infinite(range[1]) && range[1] < 0 && is.infinite(range[2]) && range[2] > 0) ||  # (-Inf, Inf)
          (is.finite(range[1]) && is.finite(range[2]) && range[1] < range[2])  # (lower, higher)
        )
    ) {
      TRUE  # range is valid
    } else {
      stop("Range is not specified correctly. Check whether vector with two numbers ",
           "(or -Inf/Inf) in increasing order.")
    }

  } else if (method %in% c("Value")) {
    if (!is.numeric(value)) {
      stop("Value is non-numeric.")
    }
    if (length(value) != 1) {
      stop("Length of value parameter does not equal 1.")
    }
  }

  # Apply transformation based on method
  if (method %in% c("Empirical", "Value")) {
    x_empirical <- sample_population(
      X = X, method = method, value = value, range = range,
      grp_var = grp_var, population_var = population_var
    )

    m_theta_a <- estimate_theta_empirical(
      est_pars = beta_draw_pg,
      X = x_empirical[["x_a"]]
    )
    m_theta_b <- estimate_theta_empirical(
      est_pars = beta_draw_pg,
      X = x_empirical[["x_b"]]
    )

  } else if (method == "Analytical") {
    m_theta_a <- estimate_theta_analytical(
      est_pars = beta_draw_pg,
      X = X,
      population_var = population_var,
      grp_var = grp_var,
      grp = 0,
      range_x = range
    )
    m_theta_b <- estimate_theta_analytical(
      est_pars = beta_draw_pg,
      X = X,
      population_var = population_var,
      grp_var = grp_var,
      grp = 1,
      range_x = range
    )
  }

  return(list(m_theta_a = m_theta_a, m_theta_b = m_theta_b))
}

#' Print Method for bglm Objects
#'
#' @param x A bglm object.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#' @export
print.bglm <- function(x, digits = 3, ...) {

  cat("\nBayesian Multivariate Logistic Regression\n")
  cat("============================================\n\n")

  cat("Group Estimates:\n")
  est_tab <- data.frame(
    Group   = c(x$info$grp_a, x$info$grp_b),
    "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
    "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
    "sd y1"   = c(x$estimates$sd_a[1], x$estimates$sd_b[1]),
    "sd y2"   = c(x$estimates$sd_a[2], x$estimates$sd_b[2]),
    check.names = FALSE # Allows spaces in column names
  )
  print(est_tab, row.names = FALSE)
  cat("\n")

  cat("Multinomial Regression Parameters:\n")
  par_names <- c("b11", "b10", "b01", "b00")
  var_names <- c("Intercept", x$info$grp_var, x$info$x_var, paste0(x$info$grp_var, "_", x$info$x_var))
  rc_tab <- rbind.data.frame(
    matrix(rbind(x$estimates$b, x$estimates$b_sd), nrow = nrow(x$estimates$b))
  )
  rc_colnames <- expand.grid(c("b_", "sd_"), par_names)
  colnames(rc_tab) <- c(paste0(rc_colnames$Var1, rc_colnames$Var2))
  row.names(rc_tab) <- var_names
  print(rc_tab)
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
  cat(sprintf(" Prior: normal distribution with :\n"))
  cat(sprintf(" Mean: \n"))
  df_b_mu0 <- setNames(as.data.frame(t(x$info$b_mu0)), var_names)
  rownames(df_b_mu0) <- NULL
  print(df_b_mu0)
  cat("\n")

  cat(sprintf(" Variance: \n"))
  df_b_sigma0 <- setNames(as.data.frame(cbind(var_names, solve(x$info$b_sigma0))), c("",var_names))
  rownames(df_b_sigma0) <- NULL
  print(df_b_sigma0)
  cat("\n")

  cat(sprintf("  Decision rule: %s\n", x$info$rule))
  if (x$info$rule == "Comp") {
    cat(sprintf("  Weights: %s\n", paste(x$info$w, collapse = ", ")))
  }
  cat(sprintf("  Hypothesis: %s\n", x$info$test_label))
  cat(sprintf("  Marginalization method: %s\n", x$info$marg_method))
  cat(sprintf("  Reference group: %s\n", x$info$ref_grp))
  cat(sprintf("  (Sub)population: [%s]\n\n",
              paste(x$info$sub_pop, collapse = ", ")))

  invisible(x)
}
