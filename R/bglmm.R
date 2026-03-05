#### bglmm.R - Bayesian Generalized Linear Mixed Model Functions ####

#' Sample Beta Coefficients for Multilevel Model
#'
#' Sample regression coefficients via Polya-Gamma multinomial logistic regression
#' for multilevel data (single chain).
#'
#' @param X List of J (n_j x P) design matrices.
#' @param Y List of J (n_j x Q) response matrices.
#' @param fixed Character vector with names of fixed variables in covariate vector.
#' @param random Character vector with names of random variables in covariate vector.
#' @param n_burn Scalar. Number of burnin iterations.
#' @param n_it Scalar. Number of iterations.
#' @param start Vector of two starting values, each used for a different chain.
#' @param b_mu0 Vector of prior means of fixed regression coefficients (length = no. of fixed covariates).
#' @param b_sigma0 Prior precision matrix of fixed regression coefficients (length(fixed) x length(fixed)).
#' @param g_mu0 Vector of prior means of random regression coefficients (length = no. of random covariates).
#' @param g_sigma0 Prior precision matrix of random regression coefficients (length(random) x length(random)).
#' @param nu0 Scalar. Degrees of freedom of prior inverse-Wishart distribution.
#' @param tau0 Prior matrix of inverse-Wishart distribution.
#' @param return_thinned Logical. Should thinned chains be returned? Default is FALSE.
#' @param n_thin Thinning rate. Default is 1. Adjustable to integers > 1 when return_thinned = TRUE.
#'
#' @return A named list with:
#' \describe{
#'   \item{b_draw_pg}{List of length n_it/n_thin with fixed regression coefficients (if fixed effects present)}
#'   \item{g_draw_pg}{List of length n_it/n_thin with random regression coefficients (if random effects present)}
#'   \item{gj_draw_pg}{List of length n_it/n_thin with cluster-specific random effects (if return_thinned = TRUE)}
#'   \item{tau_draw_pg}{List of length n_it/n_thin with covariance matrices of random effects}
#' }
#' @keywords internal
sample_beta_ml <- function(X, Y, fixed, random, n_burn, n_it, start,
                          b_mu0 = NULL, b_sigma0 = NULL, g_mu0 = NULL,
                          g_sigma0 = NULL, nu0 = NULL, tau0 = NULL,
                          return_thinned = FALSE, n_thin = 1) {

  out <- tryCatch(
    {
      J <- length(X)
      z <- lapply(1:J, function(j) Y[[j]] - 1/2)

      # Fixed effects setup
      p_f <- length(fixed)
      if (p_f > 0) {
        ind_f <- 1:p_f
        X_f <- lapply(1:J, function(j) {
          matrix(X[[j]], ncol = length(c(fixed, random)))[, ind_f, drop = FALSE]
        })
        p_b0 <- b_sigma0 %*% b_mu0
      }

      # Random effects setup
      p_r <- length(random)
      if (p_r > 0) {
        ind_r <- 1:p_r + p_f
        X_r <- lapply(1:J, function(j) {
          matrix(X[[j]], ncol = length(c(fixed, random)))[, ind_r, drop = FALSE]
        })
        p_g0 <- g_sigma0 %*% g_mu0
      }

      n <- sapply(Y, nrow)
      Q <- ncol(Y[[1]])

      # Initialize storage
      if (p_f > 0) {
        b_draw_pg <- vector("list", n_it)
      }
      if (p_r > 0) {
        g_draw_pg <- gj_draw_pg <- tau_draw_pg <- omega_draw_pg <- w_g_draw_pg <-
          list("vector", n_it)
      }

      # Initialize parameters
      omega_draw <- lapply(1:(Q - 1), function(q) {
        lapply(1:J, function(j) rep(start, n[j]))
      })

      if (p_f > 0) {
        b_draw <- cbind(matrix(start, nrow = p_f, ncol = Q - 1), 0)
      }

      if (p_r > 0) {
        g_draw <- cbind(matrix(start, nrow = p_r, ncol = Q - 1), 0)
        gj_draw <- w_g_draw <- abind::abind(
          array(start, dim = c(p_r, Q - 1, J)),
          matrix(0, p_r, J),
          along = 2
        )
        var_mat <- lapply(1:(Q - 1), function(q) diag(start, p_r))
        tau_mat <- lapply(var_mat, function(x) chol2inv(chol(x)))
      }

      # Create progress bar
      pb <- tcltk::tkProgressBar(
        title = "MCMC Progress",
        min = 1,
        max = n_burn + n_it,
        width = 300
      )

      # Start Gibbs sampler
      for (i in 2:(n_burn + n_it)) {
        for (q in 1:(Q - 1)) {

          # Draw auxiliary variable (variance of kappa)
          RC <- lapply(1:J, function(j) {
            rbind(
              if (p_f > 0) {b_draw},
              if (p_r > 0) {matrix(gj_draw[, , j], nrow = p_r, ncol = Q)}
            )
          })
          C <- lapply(1:J, function(j) {
            log(rowSums(exp(X[[j]] %*% RC[[j]][, -q, drop = FALSE])))
          })
          omega_draw[[q]] <- lapply(1:J, function(j) {
            pgdraw::pgdraw(1, X[[j]] %*% RC[[j]][, q, drop = FALSE] - C[[j]])
          })

          # Draw fixed effects
          if (p_f > 0) {
            if (p_r > 0) {
              C_b <- lapply(1:J, function(j) {
                -X_r[[j]] %*% matrix(gj_draw[, q, j], nrow = p_r) +
                  log(rowSums(exp(X[[j]] %*% RC[[j]][, -q])))
              })
            } else {
              C_b <- C
            }

            b_sigma <- chol2inv(chol(
              Reduce('+', lapply(1:J, function(j) {
                t(X_f[[j]]) %*% (X_f[[j]] * omega_draw[[q]][[j]])
              })) + b_sigma0
            ))
            b_mu <- b_sigma %*% (
              Reduce('+', lapply(1:J, function(j) {
                t(X_f[[j]]) %*% (z[[j]][, q, drop = FALSE] + omega_draw[[q]][[j]] * C_b[[j]])
              })) + p_b0
            )
            b_draw[, q] <- b_mu + t(chol(b_sigma)) %*% rnorm(p_f)
          }

          # Draw random effects
          if (p_r > 0) {
            if (p_f > 0) {
              C_g <- lapply(1:J, function(j) {
                -X_f[[j]] %*% b_draw[, q, drop = FALSE] +
                  log(rowSums(exp(X[[j]] %*% rbind(
                    b_draw[, -q, drop = FALSE],
                    matrix(gj_draw[, -q, j], nrow = p_r, ncol = Q - 1)
                  ))))
              })
            } else {
              C_g <- C
            }

            # Cluster-specific random effects
            for (j in 1:J) {
              gj_sigma <- chol2inv(chol(
                t(X_r[[j]]) %*% (X_r[[j]] * omega_draw[[q]][[j]]) + tau_mat[[q]]
              ))
              gj_mu <- gj_sigma %*% (
                t(X_r[[j]]) %*% (z[[j]][, q, drop = FALSE] + omega_draw[[q]][[j]] * C_g[[j]]) +
                  tau_mat[[q]] %*% g_draw[, q, drop = FALSE]
              )
              gj_draw[, q, j] <- c(gj_mu + t(chol(gj_sigma)) %*% rnorm(p_r))
            }

            # Population-level random effects
            g_sigma <- chol2inv(chol(J * tau_mat[[q]] + g_sigma0))
            g_mu <- g_sigma %*% (
              tau_mat[[q]] %*% as.matrix(rowSums(matrix(gj_draw[, q, ], nrow = p_r, ncol = J))) + p_g0
            )
            g_draw[, q] <- g_mu + t(chol(g_sigma)) %*% rnorm(p_r)

            # Variance components
            tau_shape <- nu0 + J
            e <- matrix(NA, p_r, J)
            for (j in 1:J) {
              e[, j] <- gj_draw[, q, j] - g_draw[, q]
            }
            tau_scale <- tau0 + e %*% t(e)
            var_mat[[q]] <- MCMCpack::riwish(tau_shape, tau_scale)
            tau_mat[[q]] <- chol2inv(chol(var_mat[[q]]))  # Precision
          }
        }

        # Store draws after burn-in
        if (i > n_burn) {
          if (p_f > 0) {
            b_draw_pg[[i - n_burn]] <- b_draw
          }
          if (p_r > 0) {
            g_draw_pg[[i - n_burn]] <- g_draw
            gj_draw_pg[[i - n_burn]] <- gj_draw
            tau_draw_pg[[i - n_burn]] <- var_mat
          }
        }

        tcltk::setTkProgressBar(
          pb, i,
          label = paste("Iteration", i, "of", n_burn + n_it)
        )
      }
      close(pb)

      # Apply thinning if requested
      if (return_thinned) {
        if (p_f > 0) {
          thinned_b_draw_pg <- lapply(seq(1, n_it, n_thin), function(i) {
            b_draw_pg[[i]]
          })
          out_fixed <- list(b_draw_pg = thinned_b_draw_pg)
        }
        if (p_r > 0) {
          thinned_g_draw_pg <- lapply(seq(1, n_it, n_thin), function(i) {
            g_draw_pg[[i]]
          })
          thinned_gj_draw_pg <- lapply(seq(1, n_it, n_thin), function(i) {
            gj_draw_pg[[i]]
          })
          thinned_tau_draw_pg <- lapply(seq(1, n_it, n_thin), function(i) {
            tau_draw_pg[[i]]
          })
          out_random <- list(
            g_draw_pg = thinned_g_draw_pg,
            gj_draw_pg = thinned_gj_draw_pg,
            tau_draw_pg = thinned_tau_draw_pg
          )
        }
      } else {
        if (p_f > 0) {
          out_fixed <- list(b_draw_pg = b_draw_pg)
        }
        if (p_r > 0) {
          out_random <- list(g_draw_pg = g_draw_pg, tau_draw_pg = tau_draw_pg)
        }
      }

      return(c(if (p_f > 0) {out_fixed}, if (p_r > 0) {out_random}))
    },
    error = function(e) {
      message("Error in MCMC sampling:")
      message(e)
      return(NA)
    }
  )

  return(out)
}

#' Estimate Parameters for Multilevel Model
#'
#' Estimate regression coefficients using multiple chains of Polya-Gamma Gibbs
#' sampling for multilevel data.
#'
#' @param X List of J (n_j x P) design matrices.
#' @param Y List of J (n_j x Q) response matrices.
#' @param fixed Character vector with names of fixed variables.
#' @param random Character vector with names of random variables.
#' @param n_burn Scalar. Number of burnin iterations.
#' @param n_it Scalar. Number of iterations.
#' @param start Vector of starting values for each chain.
#' @param b_mu0 Vector of prior means of fixed regression coefficients.
#' @param b_sigma0 Prior covariance matrix of fixed regression coefficients.
#' @param g_mu0 Vector of prior means of random regression coefficients.
#' @param g_sigma0 Prior covariance matrix of random regression coefficients.
#' @param nu0 Scalar. Degrees of freedom of prior inverse-Wishart distribution.
#' @param tau0 Prior matrix of inverse-Wishart distribution.
#' @param n_chain Scalar. Number of chains. Default is 2.
#' @param return_thinned Logical. Should thinned chains be returned? Default is FALSE.
#' @param n_thin Thinning rate. Default is 1.
#'
#' @return A named list with:
#' \describe{
#'   \item{Pars}{List of parameter draws from each chain}
#' }
#' @keywords internal
estimate_parameters_ml <- function(X, Y, fixed, random, n_burn, n_it, start,
                                  b_mu0 = NULL, b_sigma0 = NULL, g_mu0 = NULL,
                                  g_sigma0 = NULL, nu0 = NULL, tau0 = NULL,
                                  n_chain = 2, return_thinned = FALSE, n_thin = 1) {

  Q <- ncol(Y[[1]])
  chain <- vector("list", n_chain)

  for (ch in 1:n_chain) {
    chain[[ch]] <- sample_beta_ml(
      X = X, Y = Y, fixed = fixed, random = random,
      n_burn = n_burn, n_it = n_it, start = start[ch],
      b_mu0 = b_mu0, b_sigma0 = b_sigma0,
      g_mu0 = g_mu0, g_sigma0 = g_sigma0,
      nu0 = nu0, tau0 = tau0,
      return_thinned = return_thinned, n_thin = n_thin
    )
  }

  return(list(Pars = chain))
}

#' Transform Regression Coefficients to Theta for Multilevel Model
#'
#' Compute theta (success probabilities) from regression coefficients for
#' multilevel data via various methods.
#'
#' @param est_pars List of n_it (P x Q) arrays with posterior regression coefficients.
#' @param X List of J (n_j x P) design matrices.
#' @param measurement_levels Vector of measurement levels per predictor ("discrete" or "continuous").
#' @param population_var Optional character string with variable name that defines
#'   the subpopulation. Required when range is not c(-Inf, Inf).
#' @param grp_var Name of variable that defines groups.
#' @param grp_lvl Vector of group level names.
#' @param method "Value" for fixed values, "Empirical" for empirical marginalization,
#'   "Analytical" for numerical integration. Default is "Empirical".
#' @param range Optional range that defines the population of interest. Default is c(-Inf, Inf).
#' @param value Optional scalar that defines the population of interest.
#'   Required when method = "Value".
#' @param fixed Character vector with names of fixed variables.
#' @param random Character vector with names of random variables.
#'
#' @return A list with:
#' \describe{
#'   \item{m_theta_a}{List of n_it vectors with multivariate probabilities for group A}
#'   \item{m_theta_b}{List of n_it vectors with multivariate probabilities for group B}
#' }
#' @keywords internal
transform2theta_lr_ml <- function(est_pars, X, measurement_levels, population_var,
                                  grp_var, grp_lvl, method = c("Empirical", "Analytical", "Value"),
                                  range = NULL, value = NULL, fixed, random) {

  method <- match.arg(method)
  n_it <- length(est_pars)
  J <- length(X)

  theta_a <- theta_b <- vector("list", J)
  nj_a <- nj_b <- rep(NA, J)

  indices <- vector("list", J)

  if (method %in% c("Empirical")) {
    for (j in 1:J) {
      if (ncol(X[[j]]) <= 2) {
        indices[[j]] <- 1:nrow(X[[j]])
      } else if (measurement_levels == "discrete" && ncol(X[[j]]) > 2) {
        indices[[j]] <- which(X[[j]][, population_var] == range)
      } else if (measurement_levels == "continuous" && ncol(X[[j]]) > 2) {
        indices[[j]] <- which(
          X[[j]][, population_var] > min(range) &
            X[[j]][, population_var] < max(range)
        )
      }

      nj_a[j] <- length(intersect(indices[[j]], which(X[[j]][, grp_var] == 0)))
      nj_b[j] <- length(intersect(indices[[j]], which(X[[j]][, grp_var] == 1)))

      x_empirical <- sample_population(
        X = X[[j]], population_var = population_var, grp_var = grp_var,
        method = method, value = value, range = range,
        fixed = fixed, random = random
      )

      theta_a[[j]] <- estimate_theta_empirical(
        est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
        X = x_empirical[["x_a"]]
      )
      theta_b[[j]] <- estimate_theta_empirical(
        est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
        X = x_empirical[["x_b"]]
      )
    }

    # Weighted average across clusters
    m_theta_a <- lapply(1:n_it, function(k) {
      colSums(do.call(rbind, lapply(which(nj_a > 0), function(j) {
        theta_a[[j]][[k]] * nj_a[j] / sum(nj_a)
      })))
    })
    m_theta_b <- lapply(1:n_it, function(k) {
      colSums(do.call(rbind, lapply(which(nj_b > 0), function(j) {
        theta_b[[j]][[k]] * nj_b[j] / sum(nj_b)
      })))
    })

  } else if (method == "Value") {
    for (j in 1:J) {
      x_empirical <- sample_population(
        X = X[[j]], population_var = population_var, grp_var = grp_var,
        method = method, value = value, range = range,
        fixed = fixed, random = random
      )
      theta_a[[j]] <- estimate_theta_empirical(
        est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
        X = x_empirical[["x_a"]]
      )
      theta_b[[j]] <- estimate_theta_empirical(
        est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
        X = x_empirical[["x_b"]]
      )
    }
    m_theta_a <- lapply(1:n_it, function(k) {
      colMeans(do.call(rbind, lapply(1:J, function(j) theta_a[[j]][[k]])))
    })
    m_theta_b <- lapply(1:n_it, function(k) {
      colMeans(do.call(rbind, lapply(1:J, function(j) theta_b[[j]][[k]])))
    })

  } else if (method == "Analytical") {
    for (j in 1:J) {
      if (ncol(X[[1]]) > 2) {
        theta_a[[j]] <- estimate_theta_analytical(
          est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
          X = do.call(rbind, X),
          grp_var = grp_var,
          population_var = population_var,
          grp = 0,
          range_x = range,
          fixed = fixed,
          random = random
        )
        theta_b[[j]] <- estimate_theta_analytical(
          est_pars = lapply(1:n_it, function(i) est_pars[[i]][[j]]),
          X = do.call(rbind, X),
          grp_var = grp_var,
          population_var = population_var,
          grp = 1,
          range_x = range,
          fixed = fixed,
          random = random
        )
      } else {
        theta_a[[j]] <- theta_b[[j]] <- NULL
      }
    }
    m_theta_a <- lapply(1:n_it, function(k) {
      colMeans(do.call(rbind, lapply(1:J, function(j) theta_a[[j]][[k]])))
    })
    m_theta_b <- lapply(1:n_it, function(k) {
      colMeans(do.call(rbind, lapply(1:J, function(j) theta_b[[j]][[k]])))
    })
  }

  return(list(m_theta_a = m_theta_a, m_theta_b = m_theta_b))
}

#' Print Method for bglmm Objects
#'
#' @param x A bglmm object.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' # Uses the pre-computed example object shipped with the package:
#' print(bglmm_fit)
#'
#' @seealso \code{\link{summary.bglmm}}
#' @export
print.bglmm <- function(x, digits = 3, ...) {

  cat("\nBayesian Multilevel Multivariate Logistic Regression\n")
  cat("======================================================\n\n")

  est_tab <- data.frame(
    Group   = c(x$info$grp_a, x$info$grp_b),
    "mean y1" = c(x$estimates$mean_a[1], x$estimates$mean_b[1]),
    "mean y2" = c(x$estimates$mean_a[2], x$estimates$mean_b[2]),
    check.names = FALSE
  )
  print(est_tab, row.names = FALSE, digits = digits)
  cat(sprintf("  J = %d clusters    n(%s) = %d    n(%s) = %d\n\n",
              x$sample_sizes$J,
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

  cat(sprintf("Marginalization: %s over [%s]\n",
              x$info$marg_method,
              paste(x$info$sub_pop, collapse = ", ")))

  # MCMC convergence at a glance
  if (!is.null(x$diags)) {
    mpsrf_vals <- vapply(
      x$diags[c("b", "g", "tau")],
      function(d) if (!is.null(d)) round(d$convergence$mpsrf, 4) else NA_real_,
      numeric(1)
    )
    mpsrf_vals <- mpsrf_vals[!is.na(mpsrf_vals)]
    nms <- c(b = "fixed", g = "random", tau = "variance")[names(mpsrf_vals)]
    cat(sprintf("MPSRF: %s\n",
                paste(sprintf("%s = %.4f", nms, mpsrf_vals), collapse = "    ")))
  }

  cat("\nUse summary() for full coefficient tables, priors and MCMC diagnostics.\n\n")

  invisible(x)
}
