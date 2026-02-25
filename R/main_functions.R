#### Main Functions - bmvb, bglm, bglmm ####

#' Bayesian Multivariate Bernoulli Test
#'
#' Perform a Bayesian test for differences between two groups on multiple binary outcomes
#' using a Multivariate Bernoulli distribution.
#'
#' @param data Data frame containing the data.
#' @param grp Character string. Name of the grouping variable.
#' @param grp_a Value of grp indicating first group.
#' @param grp_b Value of grp indicating second group.
#' @param y_vars Character vector. Names of outcome variables (currently supports 2 outcomes).
#' @param test Character. Direction of test: "left_sided" for P(A>B) or "right_sided" for P(B>A).
#'   Default is "right_sided".
#' @param rule Character. Decision rule: "All" (all outcomes favor hypothesis),
#'   "Any" (at least one outcome favors hypothesis), or "Comp" (weighted combination).
#'   Default is "All".
#' @param w Numeric vector. Weights for compensatory rule. Only used if rule = "Comp".
#'   If NULL and rule = "Comp", equal weights are used. Default is NULL.
#' @param prior_a Numeric. Prior hyperparameter (Dirichlet) for group A. Default is 0.5 (Jeffreys' prior)
#' @param prior_b Numeric. Prior hyperparameter (Dirichlet) for group B. Default is 0.5 (Jeffreys' prior).
#' @param n_it Integer. Number of MCMC iterations. Default is 10000.
#' @param return_samples Logical. Should posterior samples be returned? Default is FALSE.
#'
#' @return A list of class "bmvb" containing:
#' \describe{
#'   \item{estimates}{Posterior means and standard deviations for each group}
#'   \item{sample_sizes}{Sample sizes for each group}
#'   \item{delta}{Treatment effect estimates and posterior probability}
#'   \item{info}{Test specifications}
#'   \item{samples}{Posterior samples (if return_samples = TRUE)}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Example with simulated data
#' result <- bmvb(
#'   data = my_data,
#'   grp = "treatment",
#'   grp_a = "control",
#'   grp_b = "experimental",
#'   y_vars = c("outcome1", "outcome2"),
#'   test = "right_sided",
#'   rule = "All"
#' )
#' print(result)
#' }
bmvb <- function(data, grp, grp_a, grp_b, y_vars,
                 test = c("right_sided", "left_sided"),
                 rule = c("All", "Any", "Comp"),
                 w = NULL, prior_a = 0.5, prior_b = 0.5, n_it = 1e4,
                 return_samples = FALSE) {

  # Match arguments
  test <- match.arg(test)
  rule <- match.arg(rule)

  # Input validation and data extraction
  validated <- check_input(
    data = data,
    grp = grp,
    grp_a = grp_a,
    grp_b = grp_b,
    y_vars = y_vars,
    test = test,
    rule = rule,
    w = w,
    analysis = "bmvb",
    prior_a = prior_a,
    prior_b = prior_b,
    n_burn = NULL,
    n_it = n_it
  )

  # Extract validated data
  y_a <- validated$y_a
  y_b <- validated$y_b
  w <- validated$w
  grp_a <- validated$grp_a
  grp_b <- validated$grp_b
  Q <- 2^length(y_vars)

  # Estimate theta for each group
  theta_a <- estimate_theta_mvb(y_a, prior_alpha = prior_a, n_it = n_it)
  theta_b <- estimate_theta_mvb(y_b, prior_alpha = prior_b, n_it = n_it)

  # Compute difference based on test direction
  if (test == "right_sided") {
    delta <- theta_b - theta_a
    test_label <- paste0("P(", grp_b, " > ", grp_a, ")")
  } else if (test == "left_sided") {
    delta <- theta_a - theta_b
    test_label <- paste0("P(", grp_a, " > ", grp_b, ")")
  }

  # Apply decision rule
  if (rule == "All") {
    pop <- mean(apply(delta, 1, min) > 0)
  } else if (rule == "Any") {
    pop <- mean(apply(delta, 1, max) > 0)
  } else if (rule == "Comp") {
    pop <- mean(delta %*% w > 0)
  }

  # Prepare output
  out <- list(
    estimates = list(
      mean_a = round(colMeans(theta_a), 3),
      mean_b = round(colMeans(theta_b), 3),
      sd_a = round(apply(theta_a, 2, sd), 3),
      sd_b = round(apply(theta_b, 2, sd), 3)
    ),
    sample_sizes = list(
      n_a = validated$n_a,
      n_b = validated$n_b
    ),
    delta = list(
      mean_delta = round(colMeans(delta), 3),
      se_delta = round(sqrt(apply(delta, 2, var) / nrow(delta)), 3),
      pop = round(pop, 3)
    ),
    info = list(
      rule = rule,
      test = test,
      test_label = test_label,
      grp_a = grp_a,
      grp_b = grp_b
    )
  )

  # Add weighted delta if using compensatory rule
  if (rule == "Comp") {
    out$delta$w_delta <- round(mean(delta %*% w), 3)
    out$info$w <- w
  }

  # Add posterior samples if requested
  if (return_samples) {
    out$samples <- list(
      theta_a = theta_a,
      theta_b = theta_b,
      delta = delta
    )
  }

  class(out) <- "bmvb"
  return(out)
}


#' Bayesian Generalized Linear Model
#'
#' Perform a Bayesian test for differences between two (sub)groups on multiple binary outcomes
#' using multinomial logistic regression.
#'
#' @param data Data frame containing the data.
#' @param grp Character string. Name of the grouping variable (will be treated as factor).
#' @param grp_a Value of grp indicating first group.
#' @param grp_b Value of grp indicating second group.
#' @param x_var Character string. Name of covariate variable (currently supports single continuous or binary covariate)
#' @param y_vars Character vector. Names of outcome variables (currently supports 2 outcomes).
#' @param x_method Character. Method for handling covariate: "Analytical" (numerical integration),
#'   "Empirical" (empirical marginalization), or "Value" (specific value). Default is "Empirical".
#' @param x_def Numeric vector. Defines subpopulation: length-2 vector c(lower, upper) for
#'   "Analytical"/"Empirical", or scalar for "Value". Default is c(-Inf, Inf).
#' @param test Character. Direction of test: "left_sided" for P(A>B) or "right_sided" for P(B>A).
#'   Default is "right_sided".
#' @param rule Character. Decision rule: "All" (all outcomes favor hypothesis),
#'   "Any" (at least one outcome favors hypothesis), or "Comp" (weighted combination).
#'   Default is "All".
#' @param w  Numeric vector. Weights for compensatory rule. Only used if rule = "Comp".
#'   If NULL and rule = "Comp", equal weights are used. Default is NULL.
#' @param b_mu0 Vector of prior means of fixed regression coefficients. Default is \code{rep(0, P)}, where P refers to the number of columns in the model matrix.
#' @param b_sigma0 Prior covariance matrix (PxP) of regression coefficients. Default is \code{diag(1e-2, P)}, where P refers to the number of columns in the model matrix.
#' @param n_burn Integer. Number of burn-in iterations. Default is 10000.
#' @param n_it Integer. Number of MCMC iterations. Default is 20000.
#' @param n_thin Integer. Thinning interval. Default is 1.
#' @param n_chain Integer. Number of MCMC chains to be sampled. Default is 2.
#' @param start Numeric vector. Starting values for chains. Should have length \code{n_chain}. Default is c(0.5, 1).
#' @param return_diagnostics Logical. Return MCMC diagnostics? Default is TRUE.
#' @param return_diagnostic_plots Logical. Should diagnostic plots (traceplots, autocorrelation, density) be returned? Default is \code{FALSE}. If \code{TRUE}, diagnostics are returned by default.
#' @param return_samples Logical. Should posterior samples be returned? Default is FALSE.
#'
#' @return A list of class "bglm" containing:
#' \describe{
#'   \item{estimates}{Posterior estimates for each group and regression parameters}
#'   \item{sample_sizes}{Sample sizes for each group}
#'   \item{delta}{Treatment effect estimates and posterior probability}
#'   \item{info}{Test specifications}
#'   \item{samples}{Posterior samples (if return_samples = TRUE)}
#' }
#'
#' @export
bglm <- function(data, grp, grp_a, grp_b, x_var, y_vars,
                 x_method = c("Empirical", "Analytical", "Value"),
                 x_def = c(-Inf, Inf),
                 test = c("right_sided", "left_sided"),
                 rule = c("All", "Any", "Comp"),
                 w = NULL, b_mu0 = NULL, b_sigma0 = NULL,
                 n_burn = 1e4, n_it = 2e4, n_thin = 1, n_chain = 2,
                 start = c(0.5, 1),
                 return_diagnostics = TRUE, return_diagnostic_plots = FALSE, return_samples = FALSE) {

  test <- match.arg(test)
  rule <- match.arg(rule)
  x_method <- match.arg(x_method)

  validated <- check_input(
    data = data,
    grp = grp,
    grp_a = grp_a,
    grp_b = grp_b,
    y_vars = y_vars,
    test = test,
    rule = rule,
    w = w,
    b_mu0 = b_mu0,
    b_sigma0 = b_sigma0,
    analysis = "bglm",
    x_var = x_var,
    x_method = x_method,
    x_def = x_def,
    n_burn = n_burn,
    n_it = n_it
  )

  y_a <- validated$y_a
  y_b <- validated$y_b
  w <- validated$w
  grp_a <- validated$grp_a
  grp_b <- validated$grp_b
  x_data <- validated$x_data
  x_a <- x_data[x_data[,grp] == 0,, drop = FALSE]
  x_b <- x_data[x_data[,grp] == 1,, drop = FALSE]
  x_names <- validated$x_names
  Q <- 2^length(y_vars)
  ml <- validated$ml
  b_mu0 <- validated$b_mu0
  b_sigma0 <- validated$b_sigma0

  # Transform to multinomial
  y_mult_a <- multivariate2multinomial(y_a)
  y_mult_b <- multivariate2multinomial(y_b)

    # Draw posterior samples
  pars <- estimate_parameters_pg(
    X = rbind(x_a, x_b),
    Y = rbind(y_mult_a, y_mult_b),
    n_burn = n_burn,
    n_it = n_it,
    start = start,
    b_mu0 = b_mu0,
    b_sigma0 = b_sigma0,
    n_chain = n_chain
  )

    chains_b <- coda::mcmc.list(lapply(1:n_chain, function(chain) {
      y <- coda::mcmc(matrix(
        unlist(pars[[chain]]),
        nrow = n_it / n_thin,
        ncol = ncol(x_data) * Q,
        byrow = TRUE
      )[, -((length(x_names) * (Q - 1)) + 1:(length(x_names)))])
      colnames(y) <- unlist(lapply(c("b11", "b10", "b01"), function(s) {
        paste0(s, "_", x_names, "[", seq_along(x_names), "]")
      }), use.names = FALSE)
      return(y)
    }))
    diags_b <- diagnose_mcmc(chains_b, plots = return_diagnostic_plots)

  diags <- list(
     diags_b = diags_b
  )

  # Check convergence
  if (any(sapply(diags, function(x) x$convergence$mpsrf > 1.1))) {
    warning("MCMC chains may not have converged (MPSRF > 1.1). ",
            "Consider increasing n_burn or n_it.")
  }

  # Transform to theta
  theta <- transform2theta(
    beta_draw_pg = do.call(c, pars),
    X = x_data,
    Y = rbind(y_mult_a, y_mult_b),
    grp_var = grp,
    population_var = x_var,
    measurement_level = ml,
    method = x_method,
    range = validated$range,
    value = validated$value
  )

  theta_a <- do.call(rbind, theta$m_theta_a)
  theta_b <- do.call(rbind, theta$m_theta_b)

  # Compute difference
  if (test == "right_sided") {
    delta <- theta_b - theta_a
    test_label <- paste0("P(", grp_b, " > ", grp_a, ")")
  } else if (test == "left_sided") {
    delta <- theta_a - theta_b
    test_label <- paste0("P(", grp_a, " > ", grp_b, ")")
  }

  # Apply decision rule
  if (rule == "All") {
    pop <- mean(apply(delta, 1, min) > 0)
  } else if (rule == "Any") {
    pop <- mean(apply(delta, 1, max) > 0)
  } else if (rule == "Comp") {
    pop <- mean(delta %*% w > 0)
  }

  # Prepare output
  out <- list(
    estimates = list(
      mean_a = round(colMeans(theta_a), 3),
      mean_b = round(colMeans(theta_b), 3),
      sd_a = round(apply(theta_a, 2, sd), 3),
      sd_b = round(apply(theta_b, 2, sd), 3),
      b = Reduce("+", do.call(c, pars)) / length(do.call(c, pars)),
      b_sd = apply(simplify2array(do.call(c, pars)), 1:2, sd)
    ),

    sample_sizes = list(
      n_a = validated$n_a,
      n_b = validated$n_b
    ),
    delta = list(
      mean_delta = round(colMeans(delta), 3),
      se_delta = round(sqrt(apply(delta, 2, var) / nrow(delta)), 3),
      pop = round(pop, 3)
    ),
    info = list(
      b_mu0 = b_mu0,
      b_sigma0 = b_sigma0,
      rule = rule,
      test = test,
      test_label = test_label,
      grp_a = as.character(grp_a),
      grp_b = as.character(grp_b),
      x_var = x_var,
      grp_var = grp,
      marg_method = x_method,
      sub_pop = x_def
    )
  )

  # Add weighted delta if using compensatory rule
  if (rule == "Comp") {
    out$delta$w_delta <- round(mean(delta %*% w), 3)
    out$info$w <- w
  }

  # Add diagnostics if requested
  if (return_diagnostics | return_diagnostic_plots) {
    out$diags$b <- diags_b
  }


  # Add posterior samples if requested
  if (return_samples) {
    out$samples <- list(
      theta_a = theta_a,
      theta_b = theta_b,
      delta = delta,
      b = pars
    )
  }

  class(out) <- "bglm"
  return(out)
}

#' Bayesian Generalized Linear Mixed Model
#'
#' Perform a Bayesian test for differences between two (sub)groups on multiple binary outcomes
#' using multilevel multinomial logistic regression.
#'
#' @param data Data frame containing the data.
#' @param grp Character string. Name of the grouping variable.
#' @param grp_a Value of \code{grp} indicating first group (will be determined from factor levels if NULL).
#' @param grp_b Value of \code{grp} indicating second group (will be determined from factor levels if NULL).
#' @param id_var Character string. Name of cluster/ID variable.
#' @param x_var Character string. Name of covariate variable.
#' @param y_vars Character vector. Names of outcome variables (currently supports 2 outcomes).
#' @param x_method Character. Method for handling covariate. Default is "Empirical".
#' @param x_def Numeric. Defines subpopulation. Default is \code{c(-Inf, Inf)}.
#' @param test Character. Direction of test: "left_sided" for P(A>B) or "right_sided" for P(B>A).
#'   Default is "right_sided".
#' @param rule Character. Decision rule: "All" (all outcomes favor hypothesis),
#'   "Any" (at least one outcome favors hypothesis), or "Comp" (weighted combination).
#'   Default is "All".
#' @param w  Numeric vector. Weights for compensatory rule. Only used if rule = "Comp".
#'   If NULL and rule = "Comp", equal weights are used. Default is NULL.
#' @param n_burn Integer. Number of burn-in iterations. Default is 10000.
#' @param n_it Integer. Number of MCMC iterations. Default is 50000 (takes long running time!).
#' @param start Numeric vector. Starting values for chains. Default is \code{c(0.5, 1)}.
#' @param fixed Character vector. Names of fixed effect variables. Default is c(x_var, grp_x_var).
#' @param random Character vector. Names of random effect variables. Default is c("Intercept", grp).
#' @param b_mu0 Numeric vector. Prior means for fixed effects. Default is \code{rep(0, length(fixed))}.
#' @param b_sigma0 Matrix. Prior covariance for fixed effects. Default is \code{diag(0.1, length(fixed))}.
#' @param g_mu0 Numeric vector. Prior means for random effects. Default is \code{rep(0, length(random))}.
#' @param g_sigma0 Matrix. Prior covariance for random effects. Default is \code{diag(0.1, length(random))}.
#' @param nu0 Numeric. Prior degrees of freedom for inverse-Wishart. Default is \code{length(random)}.
#' @param tau0 Matrix. Prior scale matrix of dimension \code{length(random) x length(random)} for inverse-Wishart. Default is \code{diag(1e-1, length(random))}.
#' @param n_chain Integer. Number of MCMC chains. Default is 2.
#' @param return_thinned Logical. Return thinned chains? Default is TRUE.
#' @param n_thin Integer. Thinning interval. Default is 10.
#' @param return_diagnostics Logical. Return MCMC diagnostics? Default is TRUE.
#' @param return_diagnostic_plots Logical. Should diagnostic plots (traceplots, autocorrelation, density) be returned? Default is \code{FALSE}. If \code{TRUE}, diagnostics are returned by default.
#' @param return_samples Logical. Return posterior samples? Default is FALSE.
#'
#' @return A list of class "bglmm".
#' @export
bglmm <- function(data, grp, grp_a = NULL, grp_b = NULL, id_var, x_var, y_vars,
                  x_method = c("Empirical", "Analytical", "Value"),
                  x_def = c(-Inf, Inf),
                  test = c("right_sided", "left_sided"),
                  rule = c("All", "Any", "Comp"),
                  w = NULL,
                  n_burn = 1e4, n_it = 5e4, start = c(0.5, 1),
                  fixed = NULL,
                  random = NULL,
                  b_mu0 = NULL, b_sigma0 = NULL,
                  g_mu0 = NULL, g_sigma0 = NULL,
                  nu0 = NULL, tau0 = NULL,
                  n_chain = 2, return_thinned = TRUE, n_thin = 1e1,
                  return_diagnostics = TRUE, return_diagnostic_plots = FALSE, return_samples = FALSE) {

  x_method <- match.arg(x_method)
  test <- match.arg(test)
  rule <- match.arg(rule)

  # Input validation and data extraction
  validated <- check_input(
    data = data,
    grp = grp,
    grp_a = grp_a,
    grp_b = grp_b,
    y_vars = y_vars,
    test = test,
    rule = rule,
    w = w,
    analysis = "bglmm",
    b_mu0 = b_mu0,
    b_sigma0 = b_sigma0,
    g_mu0 = g_mu0,
    g_sigma0 = g_sigma0,
    nu0 = nu0,
    tau0 = tau0,
    fixed = fixed,
    random = random,
    x_var = x_var,
    x_method = x_method,
    x_def = x_def,
    id_var = id_var,
    n_burn = n_burn,
    n_it = n_it,
    n_thin = n_thin
  )


  w <- validated$w
  grp_a <- validated$grp_a
  grp_b <- validated$grp_b
  ml <- validated$ml
  Q <- 2^length(y_vars)
  J <- validated$J
  grp_lvl <- levels(factor(data[[grp]]))

  # Fixed / random effect structure and design matrices (built in check_input)
  fixed      <- validated$fixed
  random     <- validated$random
  x_names    <- validated$x_names
  data_list  <- validated$data_list
  x_data     <- validated$x_data
  y_data     <- validated$y_data

  # Prior
  b_mu0 <- validated$b_mu0
  b_sigma0 <- validated$b_sigma0
  g_mu0 <- validated$g_mu0
  g_sigma0 <- validated$g_sigma0
  nu0 <- validated$nu0
  tau0 <- validated$tau0

  # Sample parameters
  pars <- estimate_parameters_ml(
    X = x_data, Y = y_data, fixed = fixed, random = random,
    n_burn = n_burn, n_it = n_it, start = start,
    b_mu0 = b_mu0, b_sigma0 = b_sigma0,
    g_mu0 = g_mu0, g_sigma0 = g_sigma0,
    nu0 = nu0, tau0 = tau0, n_chain = n_chain,
    return_thinned = return_thinned, n_thin = n_thin
  )

  # Diagnose MCMC samples
     chains_g <- coda::mcmc.list(lapply(1:n_chain, function(chain) {
      y <- coda::mcmc(matrix(
        unlist(pars[["Pars"]][[chain]][["g_draw_pg"]]),
        nrow = n_it / n_thin,
        ncol = length(random) * Q,
        byrow = TRUE
      )[, -((length(random) * (Q - 1)) + 1:(length(random)))])
      colnames(y) <- unlist(lapply(c("g11", "g10", "g01"), function(s) {
        paste0(s, "_", random, "[", seq_along(random), "]")
      }), use.names = FALSE)
      return(y)
    }))

     chains_tau <- coda::mcmc.list(lapply(1:n_chain, function(chain) {
      y <- coda::mcmc(matrix(
        unlist(lapply(1:(n_it / n_thin), function(i) {
          lapply(1:(Q - 1), function(q) {
            x <- pars[["Pars"]][[chain]][["tau_draw_pg"]][[i]][[q]]
            x[lower.tri(x, diag = TRUE)]
          })
        })),
        nrow = n_it / n_thin,
        byrow = TRUE
      ))
      colnames(y) <- unlist(lapply(c("g11", "g10", "g01"), function(s) {
        name_mat <- outer(random, random, paste0)
        random_names <- name_mat[lower.tri(name_mat, diag = TRUE)]
        paste0(s, "_", random_names)
      }), use.names = FALSE)
      return(y)
    }))


    diags_g <- diagnose_mcmc(chains_g, plots = return_diagnostic_plots)
    diags_tau <- diagnose_mcmc(chains_tau, plots = return_diagnostic_plots)

    if (length(fixed) > 0) {
      chains_b <- coda::mcmc.list(lapply(1:n_chain, function(chain) {
        y <- coda::mcmc(matrix(
          unlist(pars[["Pars"]][[chain]][["b_draw_pg"]]),
          nrow = n_it / n_thin,
          ncol = length(fixed) * Q,
          byrow = TRUE
        )[, -((length(fixed) * (Q - 1)) + 1:(length(fixed)))])
        colnames(y) <- unlist(lapply(c("b11", "b10", "b01"), function(s) {
          paste0(s, "_", fixed, "[", seq_along(fixed), "]")
        }), use.names = FALSE)
        return(y)
      }))
      diags_b <- diagnose_mcmc(chains_b, plots = return_diagnostic_plots)
    }

    diags <- list(
      diags_g = diags_g,
      diags_tau = diags_tau,
      diags_b = if (length(fixed > 0)) {diags_b} else {NULL}
    )
    diags <- Filter(Negate(is.null), diags)

  # Check convergence
     if (any(sapply(diags, function(x) x$convergence$mpsrf > 1.1))) {
      warning("MCMC chains may not have converged (MPSRF > 1.1). ",
              "Consider increasing n_burn or n_it.")
     }

  # Extract draws of fixed and random parameters
  draws_rc <- lapply(1:(n_it / n_thin), function(i) {
    lapply(1:J, function(j) {
      rbind(
        if (length(fixed) > 0) {pars[["Pars"]][[1]][["b_draw_pg"]][[i]]},
        if (length(random) > 0) {pars[["Pars"]][[1]][["gj_draw_pg"]][[i]][, , j]}
      )
    })
  })

  # Transform to success probabilities
  theta <- transform2theta_lr_ml(
    est_pars = draws_rc, X = x_data,
    measurement_levels = ml, population_var = x_var, grp_var = grp, grp_lvl = grp_lvl,
    method = x_method, range = validated$range, value = validated$value, fixed = fixed, random = random
  )

  theta_a <- do.call(rbind, theta$m_theta_a)
  theta_b <- do.call(rbind, theta$m_theta_b)

  # Compute difference
  if (test == "right_sided") {
    delta <- theta_b - theta_a
    test_label <- paste0("P(", grp_b, " > ", grp_a, ")")
  } else if (test == "left_sided") {
    delta <- theta_a - theta_b
    test_label <- paste0("P(", grp_a, " > ", grp_b, ")")
  }

  # Apply decision rule
  if (rule == "All") {
    pop <- mean(apply(delta, 1, min) > 0)
  } else if (rule == "Any") {
    pop <- mean(apply(delta, 1, max) > 0)
  } else if (rule == "Comp") {
    pop <- mean(delta %*% w > 0)
  }

  if (length(fixed) > 0) {
    draws_b <- do.call(c, lapply(pars[["Pars"]], function(chain) {
        chain[["b_draw_pg"]][seq_len(n_it / n_thin)]}))
    b_draw_mean <- Reduce(`+`, draws_b) / length(draws_b)
    b_draw_sd <- sqrt(Reduce(`+`, lapply(draws_b, function(x) (x - b_draw_mean)^2)) / (length(draws_b) - 1))

    if(!return_samples){
      rm(draws_b)
    }
  }

  if (length(random) > 0) {
  draws_g <- do.call(c, lapply(pars[["Pars"]], function(chain) {
    chain[["g_draw_pg"]][seq_len(n_it / n_thin)]}))
  g_draw_mean <- Reduce(`+`, draws_g) / length(draws_g)
  g_draw_sd <- sqrt(Reduce(`+`, lapply(draws_g, function(x) (x - g_draw_mean)^2)) / (length(draws_g) - 1))

  draws_tau <- do.call(c, lapply(pars[["Pars"]], function(chain) {
    chain[["tau_draw_pg"]][seq_len(n_it / n_thin)]}))
  tau_draw_mean <- lapply(seq_len(length(draws_tau[[1]])), function(k) {
    mat <- lapply(draws_tau, `[[`, k)
    Reduce(`+`, mat) / length(mat)
  })

  tau_draw_sd <- lapply(seq_len(length(draws_tau[[1]])), function(k) {
    mat <- lapply(draws_tau, `[[`, k)
   sqrt(Reduce(`+`, lapply(mat, function(x) (x - tau_draw_mean[[k]])^2)) / (length(mat) - 1))
  })

 if(!return_samples){
  rm(draws_g, draws_tau)
 }
  }

  # Prepare output
  out <- list(
    estimates = list(
      mean_a = round(colMeans(theta_a), 3),
      mean_b = round(colMeans(theta_b), 3),
      sd_a = round(apply(theta_a, 2, sd), 3),
      sd_b = round(apply(theta_b, 2, sd), 3)
    ),

    sample_sizes = list(
      n_a = validated$n_a,
      n_b = validated$n_b,
      J = J
    ),
    delta = list(
      mean_delta = round(colMeans(delta), 3),
      se_delta = round(sqrt(apply(delta, 2, var) / nrow(delta)), 3),
      pop = round(pop, 3)
    ),
    info = list(
      b_mu0 = b_mu0,
      b_sigma0 = b_sigma0,
      g_mu0 = g_mu0,
      g_sigma0 = g_sigma0,
      nu0 = nu0,
      tau0 = tau0,
      rule = rule,
      test = test,
      test_label = test_label,
      grp_a = as.character(grp_a),
      grp_b = as.character(grp_b),
      x_var = x_var,
      grp_var = grp,
      marg_method = x_method,
      sub_pop = x_def,
      fixed = fixed,
      random = random
    ),
    diags = list()
  )

  # Add weighted delta if using compensatory rule
  if (rule == "Comp") {
    out$delta$w_delta <- round(mean(delta %*% w), 3)
    out$info$w <- w
  }

  # Add diagnostics if requested
  if (return_diagnostics | return_diagnostic_plots) {
    if (length(fixed) > 0) {
      out$diags$b <- diags_b
      out$estimates$b <- b_draw_mean
      out$estimates$b_sd <- b_draw_sd
        }
    if (length(random) > 0) {
      out$diags$g <- diags_g
      out$diags$tau <- diags_tau
      out$estimates$g <- g_draw_mean
      out$estimates$g_sd <- g_draw_sd
      out$estimates$tau <- tau_draw_mean
      out$estimates$tau_sd <- tau_draw_sd

    }
  }

  # Add posterior samples if requested
  if (return_samples) {
     samples <- list(
      theta_a = theta_a,
      theta_b = theta_b,
      delta = delta,
      draws_fixed =
        if (length(fixed) > 0) {
          draws_b
          } else {
            NULL
            },
      draws_random = draws_g,
      draws_random_var = draws_tau
    )
    out$samples <- Filter(Negate(is.null), samples)
  }

  class(out) <- "bglmm"
  return(out)
}
