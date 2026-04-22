# Sample Beta Coefficients for Multilevel Model

Sample regression coefficients via Polya-Gamma multinomial logistic
regression for multilevel data (single chain).

## Usage

``` r
sample_beta_ml(
  X,
  Y,
  fixed,
  random,
  n_burn,
  n_it,
  start,
  b_mu0 = NULL,
  b_sigma0 = NULL,
  g_mu0 = NULL,
  g_sigma0 = NULL,
  nu0 = NULL,
  tau0 = NULL,
  return_thinned = FALSE,
  n_thin = 1
)
```

## Arguments

- X:

  List of J (n_j x P) design matrices.

- Y:

  List of J (n_j x Q) response matrices.

- fixed:

  Character vector with names of fixed variables in covariate vector.

- random:

  Character vector with names of random variables in covariate vector.

- n_burn:

  Scalar. Number of burnin iterations.

- n_it:

  Scalar. Number of iterations.

- start:

  Vector of two starting values, each used for a different chain.

- b_mu0:

  Vector of prior means of fixed regression coefficients (length = no.
  of fixed covariates).

- b_sigma0:

  Prior precision matrix of fixed regression coefficients (length(fixed)
  x length(fixed)).

- g_mu0:

  Vector of prior means of random regression coefficients (length = no.
  of random covariates).

- g_sigma0:

  Prior precision matrix of random regression coefficients
  (length(random) x length(random)).

- nu0:

  Scalar. Degrees of freedom of prior inverse-Wishart distribution.

- tau0:

  Prior matrix of inverse-Wishart distribution.

- return_thinned:

  Logical. Should thinned chains be returned? Default is FALSE.

- n_thin:

  Thinning rate. Default is 1. Adjustable to integers \> 1 when
  return_thinned = TRUE.

## Value

A named list with:

- b_draw_pg:

  List of length n_it/n_thin with fixed regression coefficients (if
  fixed effects present)

- g_draw_pg:

  List of length n_it/n_thin with random regression coefficients (if
  random effects present)

- gj_draw_pg:

  List of length n_it/n_thin with cluster-specific random effects (if
  return_thinned = TRUE)

- tau_draw_pg:

  List of length n_it/n_thin with covariance matrices of random effects
