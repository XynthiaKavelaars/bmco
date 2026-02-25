#### utils.R - Shared helper functions for bayesmulticat package ####

#' Compute Pairwise Correlation
#'
#' Compute pairwise correlation from joint response probabilities.
#' Based on Dai, Ding, & Wahba (2013). Multivariate Bernoulli distribution.
#'
#' @param phi Numeric vector of length 4 with joint response probabilities.
#'
#' @return Scalar. Correlation between outcome variables.
#' @keywords internal
compute_rho <- function(phi) {
  rho <- (phi[1] * phi[4] - phi[2] * phi[3]) /
    sqrt(sum(phi[c(1, 2)]) * sum(phi[c(1, 3)]) * sum(phi[c(3, 4)]) * sum(phi[c(2, 4)]))
  return(rho)
}

#' Transform Theta to Phi
#'
#' Transform bivariate success probabilities to joint response probabilities.
#' Currently supported for K=2/Q=4 only.
#'
#' @param theta Numeric vector of bivariate success probabilities.
#' @param rho Scalar pairwise correlation between success probabilities in theta.
#'
#' @return Numeric vector of 4 joint response probabilities, summing to 1 and
#'   ordered as \code{11}, \code{10}, \code{01}, \code{00}.
#' @keywords internal
theta2phi <- function(theta, rho) {
  phi11 <- rho * sqrt(prod(theta) * prod(1 - theta)) + prod(theta)
  phi <- c(phi11, theta[1] - phi11, theta[2] - phi11, 1 - theta[1] - theta[2] + phi11)

  if (!all(phi >= 0)) {
    warning("One or more negative probabilities. Rho adjusted.")
    while (!all(phi >= 0)) {
      if (rho < 0) rho <- rho + 0.01
      if (rho > 0) rho <- rho - 0.01
      phi11 <- rho * sqrt(prod(theta) * prod(1 - theta)) + prod(theta)
      phi <- c(phi11, theta[1] - phi11, theta[2] - phi11, 1 - theta[1] - theta[2] + phi11)
    }
  }
  return(phi)
}

#' Transform Phi to Theta
#'
#' Transform joint response probabilities to bivariate success probabilities.
#' Currently supported for K=2/Q=4 only.
#'
#' @param phi Numeric vector of 4 joint response probabilities, summing to 1 and
#'   ordered as \code{11}, \code{10}, \code{01}, \code{00}.
#'
#' @return Numeric vector of bivariate success probabilities.
#' @keywords internal
phi2theta <- function(phi) {
  c(sum(phi[c(1, 2)]), sum(phi[c(1, 3)]))
}

#' Transform Theta to Weighted Treatment Difference
#'
#' Transform bivariate success probabilities to weighted treatment difference.
#' Currently supported for K=2 only.
#'
#' @param theta_a Numeric vector of bivariate success probabilities for group A.
#' @param theta_b Numeric vector of bivariate success probabilities for group B.
#' @param weights Numeric vector of length K with weights for linear combination
#'   of treatment differences (Compensatory rule).
#'
#' @return Scalar. Weighted treatment difference.
#' @keywords internal
theta2delta_w <- function(theta_a, theta_b, weights) {
  delta_w <- (theta_a - theta_b) %*% weights
  return(delta_w)
}

#' Transform Phi to Weighted Treatment Difference
#'
#' Transform joint response probabilities to weighted treatment difference.
#' Currently supported for Q=4/K=2 only.
#'
#' @param phi_a Numeric vector of Q joint response probabilities for group A,
#'   ordered as \code{11}, \code{10}, \code{01}, \code{00}.
#' @param phi_b Numeric vector of Q joint response probabilities for group B,
#'   ordered as \code{11}, \code{10}, \code{01}, \code{00}.
#' @param weights Numeric vector of length K with weights for linear combination
#'   of treatment differences (Compensatory rule).
#'
#' @return Scalar. Weighted treatment difference.
#' @keywords internal
phi2delta_w <- function(phi_a, phi_b, weights) {
  delta <- phi_a - phi_b
  delta_theta <- cbind(rowSums(delta[, c(1, 2)]), rowSums(delta[, c(1, 3)]))
  delta_w <- delta_theta %*% weights
  return(delta_w)
}

#' Transform Bivariate Binomial to Multinomial Response
#'
#' Transform bivariate binomial response data to multinomial response format.
#'
#' @param y_bv n x 2 matrix with bivariate binomial responses.
#'
#' @return n x 4 matrix with multinomial responses, ordered as \code{11}, \code{10}, \code{01}, \code{00}.
#' @keywords internal
multivariate2multinomial <- function(y_bv) {
  answers <- rev(expand.grid(rev(list(c(0, 1), c(0, 1)))))
  answers <- answers[nrow(answers):1, ]
  y_mult_vec <- apply(y_bv, 1, function(x) which(apply(answers, 1, function(y) all(y == x))))
  y_mult <- matrix(0, length(y_mult_vec), 2^ncol(y_bv))
  for (i in 1:length(y_mult_vec)) {
    y_mult[i, y_mult_vec[i]] <- 1
  }
  return(y_mult)
}

#' Round to Chosen Interval
#'
#' Function to round number up or down to a chosen interval.
#' Adapted from: https://stackoverflow.com/a/32508105
#'
#' @param x Scalar. Number to be rounded.
#' @param round_to Scalar. Interval to be rounded to. E.g. 5, to round to the next 5th number.
#' @param dir Integer. "1" for rounding up; "0" for rounding down. Defaults to 1.
#'
#' @return Scalar. Rounded number.
#' @keywords internal
round_choose <- function(x, round_to, dir = 1) {
  if (dir == 1) {
    # ROUND UP
    rounded_x <- x + (round_to - x %% round_to)
  } else {
    if (dir == 0) {
      # ROUND DOWN
      rounded_x <- x - (x %% round_to)
    }
  }
  return(rounded_x)
}

#' Sample Population Subsets
#'
#' Extract subsets of covariate data for treatment groups based on specified method.
#'
#' @param X Design matrix with covariate data. Can be (n x P) matrix or list of J
#'   (n_j x P) matrices for multilevel data.
#' @param population_var Optional character string with the variable name that defines
#'   the subpopulation. Required when range is not c(-Inf, Inf).
#' @param grp_var Character string with name of variable that defines groups.
#' @param method Character. "Empirical" for empirical marginalization, "Value" for
#'   vector of fixed values.
#' @param value Scalar. Value of x representing subpopulation. Required when method = "Value".
#' @param range Numeric vector of lower and upper bound of covariate that represents
#'   subpopulation. Required when method = "Empirical".
#' @param fixed Character vector with names of fixed variables in covariate vector.
#'   Default is NULL.
#' @param random Character vector with names of random variables in covariate vector.
#'   Default is NULL.
#'
#' @return A list with:
#' \describe{
#'   \item{x_a}{Design matrix with covariate data for group A}
#'   \item{x_b}{Design matrix with covariate data for group B}
#' }
#' @keywords internal
sample_population <- function(X, population_var = NULL, grp_var, method,
                             value = NULL, range = NULL, fixed = NULL, random = NULL) {
  if (method == "Value") {
    # Create design matrix for specific value
    Intercept <- 1
    assign(grp_var, 1)
    assign(population_var, value)
    assign(paste0(grp_var, "_", population_var), get(grp_var) * value)
    if(!is.null(fixed) & !is.null(random)){
    x_vals_b <- c(fixed, random)
    } else {
      x_vals_b <- c("Intercept", grp_var, population_var, paste0(grp_var, "_", population_var))
    }
    x_b <- t(do.call(rbind, mget(x_vals_b)))

    assign(grp_var, 0)
    assign(population_var, value)
    assign(paste0(grp_var, "_", population_var), get(grp_var) * value)
    if(!is.null(fixed) & !is.null(random)){
      x_vals_a <- c(fixed, random)
    } else {
      x_vals_a <- c("Intercept", grp_var, population_var, paste0(grp_var, "_", population_var))
    }
    x_a <- t(do.call(rbind, mget(x_vals_a)))

  } else if (method == "Empirical") {
    # Extract empirical subset
    if (ncol(X) > 2) {
      x_b <- X[X[, grp_var] == 1 & X[, population_var] >= min(range) & X[, population_var] <= max(range), ]
      x_a <- X[X[, grp_var] == 0 & X[, population_var] >= min(range) & X[, population_var] <= max(range), ]
    } else {
      x_b <- X[X[, grp_var] == 1, ]
      x_a <- X[X[, grp_var] == 0, ]
    }
  }

  return(list(x_a = x_a, x_b = x_b))
}

#' Estimate Theta via Empirical Marginalization
#'
#' Estimate success probabilities via empirical marginalization over covariate values.
#'
#' @param est_pars List of nIt (P x Q) matrices/arrays of posterior draws of
#'   regression coefficients.
#' @param X Design matrix with covariate data.
#'
#' @return List of nIt vectors of length K containing bivariate probabilities.
#'   Currently supported for K=2 only.
#' @keywords internal
estimate_theta_empirical <- function(est_pars, X) {
  n_it <- length(est_pars)

  # Compute predicted probabilities
  m_phi <- lapply(1:n_it, function(i) {
    exp(X %*% est_pars[[i]]) / rowSums(exp(X %*% est_pars[[i]]))
  })

  # Transform to theta (marginal probabilities)
  if (length(m_phi[[1]]) == 4 || ncol(m_phi[[1]]) == 4) {
    # Bivariate case (Q=4)
    m_theta <- lapply(1:n_it, function(i) {
      colMeans(cbind(
        rowSums(m_phi[[i]][, c(1, 2), drop = FALSE]),
        rowSums(m_phi[[i]][, c(1, 3), drop = FALSE])
      ))
    })
  } else if (length(m_phi[[1]]) == 2 || ncol(m_phi[[1]]) == 2) {
    # Univariate case (Q=2)
    m_theta <- lapply(1:n_it, function(i) colMeans(m_phi[[i]]))
  }

  return(m_theta)
}

#' Integrand for Analytical Marginalization
#'
#' Integrand function for integration over a range of a covariate in numerical
#' marginalization. Works for both GLM and GLMM.
#'
#' @param x Scalar. Value of covariate x.
#' @param est_rc (P x Q) matrix of regression coefficients.
#' @param mu_x Scalar. Mean of distribution of covariate.
#' @param sigma_x Scalar. Standard deviation of distribution of covariate.
#' @param range_x Numeric vector of lower and upper bound of range to integrate over.
#' @param grp Scalar. Value of group indicator (0 or 1).
#' @param q Scalar. Response category in 1 to Q.
#' @param fixed Optional character vector with names of fixed variables. Default is NULL.
#' @param random Optional character vector with names of random variables. Default is NULL.
#'
#' @return Scalar. Joint response probability for response category q.
#' @keywords internal
integrand <- function(x, est_rc, mu_x, sigma_x, range_x, grp, grp_var, population_var, q,
                     fixed = NULL, random = NULL) {
  # Create design matrix for integration point
  if (!is.null(fixed) && !is.null(random)) {
    # Multilevel case
    Intercept <- 1
    assign(grp_var, grp)
    assign(population_var, x)
    assign(paste0(grp_var, "_", population_var), grp * x)
    x_vals <- c(fixed, random)
    x_int <- t(do.call(rbind, mget(x_vals)))
  } else {
    # Standard GLM case
    x_int <- cbind(1, grp, x, grp * x)
  }

  # Compute predicted probability
  psi <- x_int %*% est_rc
  p_x <- msm::dtnorm(x, mean = mu_x, sd = sigma_x,
                     lower = range_x[1], upper = range_x[2], log = TRUE)
  phi <- exp(psi - log(rowSums(exp(psi))) + p_x)

  return(phi[, q])
}

#' Estimate Theta via Analytical Integration
#'
#' Estimate success probabilities by integrating over a range of covariate values.
#'
#' @param est_pars List of nIt (P x Q) matrices of posterior regression coefficients.
#' @param X Design matrix of covariate data, where the covariate of interest is in
#'   the third column.
#' @param grp Scalar. Value of group indicator (0 or 1).
#' @param range_x Numeric vector of lower and upper bound of range to integrate over.
#' @param fixed Optional character vector with names of fixed variables. Default is NULL.
#' @param random Optional character vector with names of random variables. Default is NULL.
#'
#' @return List of nIt vectors of length K with multivariate probabilities.
#'   Currently supported for K=2 only.
#' @keywords internal
estimate_theta_analytical <- function(est_pars, X, grp, range_x, grp_var, population_var,
                                      fixed = NULL, random = NULL) {
  Q <- ncol(est_pars[[1]])
  n_it <- length(est_pars)

  m_theta <- vector("list", n_it)
  m_phi <- array(NA, dim = c(1, Q))

  for (i in 1:n_it) {
    for (q in 1:(Q - 1)) {
      m_phi[, q] <- integrate(
        integrand,
        lower = range_x[1],
        upper = range_x[2],
        est_rc = est_pars[[i]],
        mu_x = mean(X[, population_var]),
        sigma_x = sd(X[, population_var]),
        range_x = range_x,
        grp_var = grp_var,
        population_var = population_var,
        grp = grp,
        q = q,
        fixed = fixed,
        random = random
      )$value
    }

    if (Q == 4) {
      m_theta[[i]] <- cbind(
        rowSums(m_phi[, c(1, 2), drop = FALSE]),
        rowSums(m_phi[, c(1, 3), drop = FALSE])
      )
    } else if (Q == 2) {
      m_theta[[i]] <- m_phi
    }
  }

  return(m_theta)
}

#### check_input.R - Input validation for all main functions ####

#' Check and Validate Input Parameters
#'
#' Comprehensive input validation for bmvb, bglm, and bglmm functions.
#' Performs all necessary checks and returns cleaned/validated parameters.
#'
#' @param data Data frame containing the data.
#' @param grp Character string. Name of the grouping variable.
#' @param grp_a Value of grp indicating first group.
#' @param grp_b Value of grp indicating second group.
#' @param y_vars Character vector. Names of outcome variables.
#' @param test Character. Direction of test ("right_sided" or "left_sided").
#' @param rule Character. Decision rule ("All", "Any", or "Comp").
#' @param w Numeric vector. Weights for compensatory rule (can be NULL).
#' @param analysis Character. Type of analysis: "bmvb", "bglm", or "bglmm".
#' @param prior_a Numeric. Prior for group A (bmvb only). Default is NULL.
#' @param prior_b Numeric. Prior for group B (bmvb only). Default is NULL.
#' @param b_mu0 Vector (length = no. of fixed covariates) of prior means of fixed regression coefficients (bglm/bglmm only). Default is NULL.
#' @param b_sigma0 Prior precision matrix (P_fixed x P_fixed) of fixed regression coefficients (bglm/bglmm only). Default is NULL.
#' @param g_mu0 Vector (length = no. of random covariates) of prior means of random regression coefficients (bglmm only) . Default is NULL.
#' @param g_sigma0 Prior precision matrix (P_random x P_random) of random regression coefficients (bglmm only). Default is NULL.
#' @param nu0 Scalar. Prior df for random covariance matrix (bglmm only). Default is NULL.
#' @param tau0 Numeric. Prior scale matrix (length(random) x length(random)) for random covariance matrix (bglmm only). Default is NULL.
#' @param fixed Character vector. Names of fixed effect variables. Default is c("x", "grp_x").
#' @param random Character vector. Names of random effect variables. Default is c("Intercept", grp).
#' @param x_var Character string. Name of covariate (bglm/bglmm only). Default is NULL.
#' @param x_method Character. Method for handling covariate (bglm/bglmm only). Default is NULL.
#' @param x_def Numeric. Defines subpopulation (bglm/bglmm only). Default is NULL.
#' @param id_var Character string. Name of cluster ID (bglmm only). Default is NULL.
#' @param n_burn Integer. Number of burnin iterations.
#' @param n_it Integer. Number of MCMC iterations.
#' @param n_thin Integer. Thinning interval (bglmm only). Default is 1.
#' @importFrom stats contr.treatment "contrasts<-" relevel setNames

#' @return List with validated parameters and cleaned data subsets:
#' \describe{
#'   \item{y_a}{Matrix of outcomes for group A (with missing data removed)}
#'   \item{y_b}{Matrix of outcomes for group B (with missing data removed)}
#'   \item{w}{Weights (generated if NULL and rule = "Comp")}
#'   \item{grp_a}{Validated group A value}
#'   \item{grp_b}{Validated group B value}
#'   \item{J}{Number of clusters (bglmm only)}
#' }
#' @keywords internal
check_input <- function(data, grp, grp_a, grp_b, y_vars, test, rule, w,
                        analysis = c("bmvb", "bglm", "bglmm"),
                        prior_a = NULL, prior_b = NULL,
                        b_mu0 = NULL, b_sigma0 = NULL,
                        g_mu0 = NULL, g_sigma0 = NULL,
                        nu0 = NULL, tau0 = NULL,
                        fixed = NULL, random = NULL,
                        x_var = NULL, x_method = NULL, x_def = NULL,
                        id_var = NULL, n_burn = NULL, n_it, n_thin = 1) {

  # Match analysis type
  analysis <- match.arg(analysis)

  # ============================================================================
  # 1. BASIC EXISTENCE CHECKS
  # ============================================================================

  # Check that data is a data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Check group variable exists
  if (!grp %in% names(data)) {
    stop("Group variable '", grp, "' not found in data.")
  }

  # Check outcome variables exist
  if (!all(y_vars %in% names(data))) {
    missing_vars <- y_vars[!y_vars %in% names(data)]
    stop("Outcome variable(s) not found in data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Check group values exist
  if (!grp_a %in% data[[grp]]) {
    stop("Group value '", grp_a, "' not found in grouping variable '", grp, "'.")
  }

  if (!grp_b %in% data[[grp]]) {
    stop("Group value '", grp_b, "' not found in grouping variable '", grp, "'.")
  }

  # Check if group is factor and assign reference level grp_a
    data[,grp] <- factor(data[,grp])
    data[,grp] <- relevel(data[,grp], ref = grp_a)

  # ============================================================================
  # 2. OUTCOME VARIABLE CHECKS
  # ============================================================================

  # Check number of outcomes
  if (length(y_vars) > 2) {
    stop("Currently only supported for 2 outcomes. ",
         "Please provide exactly 2 outcome variables.")
  }

  # Check for binary outcomes (0/1)
  for (var in y_vars) {
    unique_vals <- unique(data[[var]][!is.na(data[[var]])])
    if (!all(unique_vals %in% c(0, 1))) {
      stop("Outcome variable '", var, "' must be binary (0/1). ",
           "Found values: ", paste(unique_vals, collapse = ", "))
    }
  }

  # Check for outcome variation (warn if all 0 or all 1)
  for (var in y_vars) {
    non_missing <- data[[var]][!is.na(data[[var]])]
    if (all(non_missing == 0)) {
      warning("Outcome variable '", var, "' has no variation (all 0). ",
              "This will cause poor MCMC convergence and unreliable results. ",
              "Consider removing this outcome or collecting more data.")
    } else if (all(non_missing == 1)) {
      warning("Outcome variable '", var, "' has no variation (all 1). ",
              "This will cause poor MCMC convergence and unreliable results. ",
              "Consider removing this outcome or collecting more data.")
    }
  }


  # ============================================================================
  # 3. EXTRACT AND CLEAN GROUP DATA
  # ============================================================================

  # Extract data for each group
  y_a <- as.matrix(data[data[[grp]] == grp_a, y_vars, drop = FALSE])
  y_b <- as.matrix(data[data[[grp]] == grp_b, y_vars, drop = FALSE])

  # Check for and remove missing data
  if (any(is.na(y_a))) {
    n_missing <- sum(!complete.cases(y_a))
    warning("Missing data detected in group ", grp_a, " (", n_missing,
            " rows). Rows with missing data will be removed.")
    y_a <- y_a[complete.cases(y_a), , drop = FALSE]
  }

  if (any(is.na(y_b))) {
    n_missing <- sum(!complete.cases(y_b))
    warning("Missing data detected in group ", grp_b, " (", n_missing,
            " rows). Rows with missing data will be removed.")
    y_b <- y_b[complete.cases(y_b), , drop = FALSE]
  }

  # Get sample sizes
  n_a <- nrow(y_a)
  n_b <- nrow(y_b)

  # Check for outcome variation within each group
  for (i in seq_along(y_vars)) {
    var <- y_vars[i]
    if (n_a > 0 && all(y_a[, i] == y_a[1, i])) {
      warning("Outcome variable '", var, "' has no variation in group ", grp_a,
              " (all ", y_a[1, i], "). This may cause convergence issues.")
    }
    if (n_b > 0 && all(y_b[, i] == y_b[1, i])) {
      warning("Outcome variable '", var, "' has no variation in group ", grp_b,
              " (all ", y_b[1, i], "). This may cause convergence issues.")
    }
  }

  # Check for sufficient sample size
  if (n_a < 10) {
    warning("Small sample size in group ", grp_a, " (n = ", n_a,
            "). Results may be unreliable.")
  }
  if (n_b < 10) {
    warning("Small sample size in group ", grp_b, " (n = ", n_b,
            "). Results may be unreliable.")
  }

  # Check that groups are not empty
  if (n_a == 0) {
    stop("No observations in group ", grp_a, " after removing missing data.")
  }
  if (n_b == 0) {
    stop("No observations in group ", grp_b, " after removing missing data.")
  }

  # ============================================================================
  # 4. DECISION RULE AND WEIGHTS CHECKS
  # ============================================================================

  # Handle weights for compensatory rule
  if (rule == "Comp") {
    if (is.null(w)) {
      w <- rep(1 / length(y_vars), length(y_vars))
      message("No weights specified for compensatory rule. Using equal weights: ",
              paste(round(w, 3), collapse = ", "))
    } else {
      # Check weight vector length
      if (length(w) != length(y_vars)) {
        stop("Length of weight vector (", length(w),
             ") not equal to number of outcomes (", length(y_vars), ")")
      }

      # Check weights are non-negative
      if (any(w < 0)) {
        stop("Weights must be non-negative")
      }

      # Warn if weights don't sum to 1
      if (abs(sum(w) - 1) > 0.01) {
        warning("Weights do not sum to 1 (sum = ", round(sum(w), 3), "). ",
                "Consider normalizing: w / sum(w)")
      }
    }
  }

  # ============================================================================
  # 5. MCMC CHECKS
  # ============================================================================

  if (n_it == 0) {
    stop("Number of iterations is zero.")
  }

  if(!is.null(n_burn) && n_burn == 0) {
   stop("Number of burnin iterations is zero.")
  }
  # ============================================================================
  # 6. BMVB-SPECIFIC CHECKS
  # ============================================================================

  if (analysis == "bmvb") {
    # Check prior values are positive
    if (!is.null(prior_a) && prior_a <= 0) {
      stop("prior_a must be positive (found: ", prior_a, ")")
    }
    if (!is.null(prior_b) && prior_b <= 0) {
      stop("prior_b must be positive (found: ", prior_b, ")")
    }
  }

  # ============================================================================
  # 7. BGLM AND BGLMM-SPECIFIC CHECKS
  # ============================================================================

  if (analysis %in% c("bglm", "bglmm")) {

    # Check covariate exists
    if (is.null(x_var)) {
      stop("x_var must be specified for ", analysis, " analysis")
    }

    if (!x_var %in% names(data)) {
      stop("Covariate '", x_var, "' not found in data.")
    }

    # Determine and validate measurement level
    if (is.factor(data[[x_var]]) || is.character(data[[x_var]])) {
      ml <- "discrete"
    } else if (is.numeric(data[[x_var]])) {
      ml <- "continuous"
    } else {
      stop("Covariate '", x_var, "' has unknown measurement level. ",
           "Provide as factor (discrete) or numeric (continuous).")
    }

    # For discrete covariates, check maximum 2 categories
    if (ml == "discrete") {
      n_categories <- length(unique(data[[x_var]][!is.na(data[[x_var]])]))
      if (n_categories > 2) {
        stop("Discrete covariate '", x_var, "' has ", n_categories,
             " categories. Currently only binary covariates (2 categories) are supported.")
      }
      if (n_categories < 2) {
        stop("Discrete covariate '", x_var, "' has only ", n_categories,
             " category. At least 2 categories are required.")
      }
    }

    # Validate x_def based on x_method
    if (!is.null(x_method)) {
      if (length(x_def) == 1) {
        if (x_method != "Value") {
          stop("When x_method = 'Value', x_def should have length 1. ",
               "For x_method = 'Empirical' or 'Analytical', x_def should be length 2.")
        }
      } else if (length(x_def) == 2) {
        if (!(x_method %in% c("Empirical", "Analytical"))) {
          stop("When x_method = 'Empirical' or 'Analytical', x_def should be a vector of length 2. ",
               "For x_method = 'Value', x_def should be length 1.")
        }
        # Check that range is in correct order
        if (x_def[1] >= x_def[2] && is.finite(x_def[1]) && is.finite(x_def[2])) {
          stop("x_def range must be in increasing order: [lower, upper]")
        }
      }
    }

     # Check for zero variation in covariate
    if (is.numeric(data[[x_var]])) {
      covar_sd <- sd(data[[x_var]], na.rm = TRUE)

      if (covar_sd == 0 || is.na(covar_sd)) {
        warning("Covariate '", x_var, "' has no variation (all values are identical). ",
                "This will cause problems in regression.")
      } else {
        # Check for collinearity (convert to numeric for correlation)
        grp_numeric <- as.numeric(as.factor(data[[grp]]))
        cor_val <- cor(grp_numeric, data[[x_var]], use = "complete.obs")
        if (!is.na(cor_val) && abs(cor_val) > 0.9) {
          warning("High correlation detected between group and covariate (r = ",
                  round(cor_val, 2), "). This may cause convergence issues.")
        }

        # Only check unique values if there IS variation
        if (ml == "continuous" && !is.null(x_method) && x_method == "Empirical") {
          n_unique <- length(unique(data[[x_var]][!is.na(data[[x_var]])]))
          if (n_unique < 5) {
            warning("Few unique values in continuous covariate (", n_unique,
                    " values). Consider treating as discrete.")
          }
        }
      }
    }

     # Validate range is within data (only for continuous)
    if (ml == "continuous" && !is.null(x_method) &&
        x_method %in% c("Empirical", "Analytical")) {
      data_range <- range(data[[x_var]], na.rm = TRUE)

      # Check for extreme values or wide range first (applies to ALL methods)
      covar_range_width <- data_range[2] - data_range[1]
      if (covar_range_width > 1e4) {
        stop("Covariate '", x_var, "' has very wide range (",
             round(covar_range_width, 0), "). ",
             "This causes numerical instability in MCMC. ",
             "Consider: (1) standardizing the covariate, ",
             "(2) using a subset with narrower range, or (3) using bmvb().")
      }

      if (abs(data_range[1]) > 1e4 || abs(data_range[2]) > 1e4) {
        stop("Covariate '", x_var, "' has extreme values (range: ",
             round(data_range[1], 0), " to ", round(data_range[2], 0), "). ",
             "This causes numerical instability in MCMC. ",
             "Consider: (1) standardizing the covariate, ",
             "(2) centering around the mean, or (3) using bmvb().")
      }

      # Check x_def bounds
      if (is.finite(x_def[1]) && x_def[1] < data_range[1]) {
        warning("Lower bound of x_def (", x_def[1],
                ") is below minimum observed value (", round(data_range[1], 2), ")")
      }
      if (is.finite(x_def[2]) && x_def[2] > data_range[2]) {
        warning("Upper bound of x_def (", x_def[2],
                ") is above maximum observed value (", round(data_range[2], 2), ")")
      }

      # For Empirical method, x_def must overlap with data
      if (x_method == "Empirical") {
        if (is.finite(x_def[1]) && x_def[1] > data_range[2]) {
          stop("Lower bound of x_def (", x_def[1],
               ") is above maximum observed value (", round(data_range[2], 2),
               "). No data in specified range. Use x_method = 'Analytical'.")
        }
        if (is.finite(x_def[2]) && x_def[2] < data_range[1]) {
          stop("Upper bound of x_def (", x_def[2],
               ") is below minimum observed value (", round(data_range[1], 2),
               "). No data in specified range. Use x_method = 'Analytical'.")
        }
      }
    }

    # Also check extreme values for Value method (continuous only)
    if (ml == "continuous" && !is.null(x_method) && x_method == "Value") {
      data_range <- range(data[[x_var]], na.rm = TRUE)

      # Same extreme value checks
      covar_range_width <- data_range[2] - data_range[1]
      if (covar_range_width > 1e4) {
        stop("Covariate '", x_var, "' has very wide range (",
             round(covar_range_width, 0), "). ",
             "This causes numerical instability in MCMC. ",
             "Consider standardizing the covariate or using bmvb().")
      }

      if (abs(data_range[1]) > 1e4 || abs(data_range[2]) > 1e4) {
        stop("Covariate '", x_var, "' has extreme values (range: ",
             round(data_range[1], 0), " to ", round(data_range[2], 0), "). ",
             "This causes numerical instability in MCMC. ",
             "Consider standardizing the covariate or using bmvb().")
      }
    }


    # Suggest thinning for large iterations
    if (n_it > 1e6 && n_thin == 1) {
      message("Consider using thinning (n_thin > 1) with very large n_it ",
              "to reduce time and memory usage")
    }

  # Derive value and range from x_def based on x_method.
  # These are passed to transform2theta() / transform2theta_lr_ml().
  if (!is.null(x_method) && x_method == "Value") {
    pop_value <- x_def          # scalar to evaluate at
    pop_range <- c(-Inf, Inf)   # not used for Value, kept as safe default
  } else {
    pop_value <- NULL           # not used for Empirical / Analytical
    pop_range <- x_def          # integration / subsetting range
  }

  # Create prior parameters if needed
  if(analysis == "bglm"){
    fixed <- c("Intercept", "group", "x", "group_x")
  } else if(analysis == "bglmm") {
    if(is.null(fixed)) {
      fixed  <- c(x_var, paste0(grp, "_", x_var))
    }
  }

  if(is.null(b_mu0)){
    b_mu0 <- rep(0, length(fixed))
  } else if (length(b_mu0) == 1){
    b_mu0 <- rep(b_mu0, length(fixed))
  } else if (!(length(b_mu0) %in% c(0, length(fixed)))){
    stop("b_mu0 must have length(fixed)")
  }

   if(is.null(b_sigma0)){
    b_sigma0 <- as.matrix(diag(1e-1, length(fixed)))
  } else if (is.vector(b_sigma0) && length(b_sigma0) == 1){
    b_sigma0 <- as.matrix(diag(b_sigma0, length(fixed)))
  } else if (is.vector(b_sigma0) && length(b_sigma0) == length(fixed)){
    b_sigma0 <- as.matrix(diag(b_sigma0, length(fixed)))
  } else if (is.vector(b_sigma0) && length(b_sigma0) > length(fixed)){
    stop("b_sigma0 must be a vector of length(fixed) or a matrix with dimensions length(fixed) x length(fixed).")
  } else if(is.matrix(b_sigma0) && !identical(dim(b_sigma0), c(length(fixed), length(fixed)))){
    stop("b_sigma0 must be a vector of length(fixed) or a matrix with dimensions length(fixed) x length(fixed).")
  }

  if(any(b_sigma0 < 0)){
    stop("b_sigma0 must be non-negative (found: ", b_sigma0, ")")
  }
}
  # ============================================================================
  # 8. BGLMM-SPECIFIC CHECKS
  # ============================================================================

  if (analysis == "bglmm") {

    # Check cluster ID exists
    if (is.null(id_var)) {
      stop("id_var must be specified for bglmm analysis")
    }

    if (!id_var %in% names(data)) {
      stop("Cluster ID variable '", id_var, "' not found in data.")
    }

    # Get number of clusters
    J <- length(unique(data[[id_var]]))

    # Check for sufficient clusters
    if (J < 5) {
      warning("Very few clusters (J = ", J, "). Random effects may not be reliable. ",
              "Consider using bglm() instead.")
    }

    # Check cluster sizes
    cluster_sizes <- table(data[[id_var]])
    min_cluster_size <- min(cluster_sizes)
    if (min_cluster_size < 3) {
      warning("Some clusters have very few observations (minimum = ",
              min_cluster_size, "). This may affect random effect estimates.")
    }

    # Check that each cluster has observations from both groups
    cluster_group_counts <- table(data[[id_var]], data[[grp]])
    clusters_with_one_group <- rowSums(cluster_group_counts > 0) == 1
    if (any(clusters_with_one_group)) {
      n_unbalanced <- sum(clusters_with_one_group)
      warning(n_unbalanced, " cluster(s) have observations from only one group. ",
              "This may affect estimation.")
    }

    # Define random variables
    if(is.null(random)) {
      random <- c("Intercept", grp)
    }

    # Create prior parameters if needed
    if(is.null(g_mu0)){
      g_mu0 <- rep(0, length(random))
    } else if (length(g_mu0) == 1){
      g_mu0 <- rep(g_mu0, length(random))
    } else if (length(g_mu0) > length(random)){
      stop("g_mu0 must have length(random)")
    }

    if(is.null(g_sigma0)){
      g_sigma0 <- diag(1e-1, length(random))
    } else if (is.vector(g_sigma0) && length(g_sigma0) == 1){
      g_sigma0 <- diag(g_sigma0, length(random))
    } else if (is.vector(g_sigma0) && length(g_sigma0) == length(random)){
      g_sigma0 <- diag(g_sigma0, length(random))
    } else if (is.vector(g_sigma0) && length(g_sigma0) > length(random)){
      stop("g_sigma0 must be a vector of length(random) or a matrix with dimensions length(random) x length(random).")
    } else if(is.matrix(g_sigma0) && !identical(dim(g_sigma0), c(length(random), length(random)))){
      stop("g_sigma0 must be a vector of length(random) or a matrix with dimensions length(random) x length(random).")
    }

    if(any(g_sigma0 < 0)){
      stop("g_sigma0 must be non-negative (found: ", g_sigma0, ")")
    }

    if(is.null(nu0)){
      nu0 <- length(random)
    } else if (length(nu0) > 1){
      stop("nu0 must have length 1")
    }

    if(any(nu0) < 0){
      stop("nu0 must be non-negative (found: ", nu0, ")")
    }

    if(is.null(tau0)){
      tau0 <- diag(1e-1, length(random))
    } else if (is.vector(tau0) && length(tau0) == length(random)){
      tau0 <- diag(tau0, length(random))
    } else if (is.vector(tau0) && length(tau0) == 1){
      tau0 <- diag(tau0, length(random))
    } else if (is.vector(tau0) && length(tau0) > length(random)){
      stop("tau0 must be a vector of length(random) or a matrix with dimensions length(random) x length(random).")
    } else if(is.matrix(tau0) && !identical(dim(tau0), c(length(random), length(random)))){
      stop("tau0 must be a vector of length(random) or a matrix with dimensions length(random) x length(random).")
    }

    if(any(tau0 < 0)){
      stop("tau0 must be non-negative (found: ", tau0, ")")
    }

  } else {
    J <- NULL
  }


  # ============================================================================
  # 9. DESIGN MATRIX CREATION (bglm and bglmm)
  # ============================================================================

  if (analysis == "bglm") {

    # Identify complete cases across both x and y
    complete_rows <- complete.cases(data[, c(x_var, y_vars)])
    n_missing_x <- sum(!complete.cases(data[, x_var, drop = FALSE]) &
                         complete.cases(data[, y_vars, drop = FALSE]))
    if (n_missing_x > 0) {
      warning("Missing data detected in covariate '", x_var, "' (", n_missing_x,
              " additional rows removed due to missing covariate values).")
    }
    data_clean <- data[complete_rows, ]

    # Order factor levels explicitly, so that grp_a == 0 and grp_b == 1 in design matrix
    data_clean[[grp]] <- factor(data_clean[[grp]], levels = c(grp_a, grp_b))
    contrasts(data_clean[[grp]]) <- contr.treatment(2, base = 1)

    # Design matrix
    x_data <- model.matrix(
        ~ 1 + data_clean[, grp] + data_clean[, x_var] + data_clean[, grp]:data_clean[, x_var]
      )
    colnames(x_data) <- x_names <- c("Intercept", grp, x_var, paste0(grp, "_", x_var))

    # Recreate y_a and y_b from the jointly clean dataset
    y_a <- data_clean[data_clean[[grp]] == grp_a, y_vars, drop = FALSE]
    y_b <- data_clean[data_clean[[grp]] == grp_b, y_vars, drop = FALSE]
    n_a <- nrow(y_a)
    n_b <- nrow(y_b)
  }

  if (analysis == "bglmm") {

    # Ensure id_var is a factor so levels() is meaningful
    if (!is.factor(data[[id_var]])) {
      data[[id_var]] <- factor(data[[id_var]])
    }

    x_names <- c(fixed, random)

    # Warn about missing covariate values (silently dropped per cluster below)
    n_missing_x <- sum(is.na(data[[x_var]]))
    if (n_missing_x > 0) {
      warning("Missing data detected in covariate '", x_var, "' (", n_missing_x,
              " rows). Rows with missing covariate values will be removed per cluster.")
    }

    # Per-cluster datasets: keep only complete cases on x AND y
    data_list <- lapply(seq_len(J), function(j) {
      data_j <- data[data[[id_var]] == levels(data[[id_var]])[j], ]
      data_j[complete.cases(data_j[, c(x_var, y_vars)]), ]
    })

    # All group levels from the full dataset (needed so per-cluster model.matrix
    # does not crash when a cluster only contains one of the two groups)
    all_grp_levels <- levels(factor(data[[grp]]))

    # Per-cluster design matrices (fixed columns first, then random)
    x_data <- lapply(seq_len(J), function(j) {
      data_j <- data_list[[j]]
      # Use explicit factor with all levels so contrasts work even in
      # single-group clusters
      grp_j <- factor(data_j[[grp]], levels = all_grp_levels)
      x_j <- model.matrix(
        ~ 1 + grp_j + data_j[, x_var] + grp_j:data_j[, x_var]
      )
      colnames(x_j) <- c("Intercept", grp, x_var, paste0(grp, "_", x_var))
      x_j[, c(fixed, random), drop = FALSE]
    })

    # Per-cluster multinomial response matrices
    y_data <- lapply(seq_len(J), function(j) {
      data_j <- data_list[[j]]
      multivariate2multinomial(as.matrix(data_j[, y_vars, drop = FALSE]))
    })

    # Update n_a and n_b to reflect complete-case filtering
    n_a <- sum(sapply(data_list, function(d) sum(d[[grp]] == grp_a)))
    n_b <- sum(sapply(data_list, function(d) sum(d[[grp]] == grp_b)))

    # Update y_a, y_b for consistent reporting
    y_a <- do.call(rbind, lapply(data_list, function(d) {
      as.matrix(d[d[[grp]] == grp_a, y_vars, drop = FALSE])
    }))
    y_b <- do.call(rbind, lapply(data_list, function(d) {
      as.matrix(d[d[[grp]] == grp_b, y_vars, drop = FALSE])
    }))
  }

  # ============================================================================
  # 10. RETURN VALIDATED PARAMETERS
  # ============================================================================

  result <- list(
    y_a   = y_a,
    y_b   = y_b,
    n_a   = n_a,
    n_b   = n_b,
    w     = w,
    grp_a = as.character(grp_a),
    grp_b = as.character(grp_b)
  )

  if (analysis == "bglmm") {
    result$J         <- J
    result$fixed     <- fixed
    result$random    <- random
    result$data_list <- data_list
    result$y_data    <- y_data
    result$g_mu0    <- g_mu0
    result$g_sigma0 <- g_sigma0
    result$nu0      <- nu0
    result$tau0   <- tau0
  }


  if (analysis %in% c("bglm", "bglmm")) {
    result$x_data   <- x_data
    result$x_names  <- x_names
    result$ml <- ml
    result$value  <- pop_value
    result$range  <- pop_range
    result$b_mu0  <- b_mu0
    result$b_sigma0 <- b_sigma0
  }


  return(result)
}



#' Compute Credible Intervals
#' @keywords internal
credible_interval <- function(samples, prob = 0.95) {
  alpha <- 1 - prob
  quantile(samples, probs = c(alpha/2, 1 - alpha/2))
}
