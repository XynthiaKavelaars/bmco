# Estimate Parameters for Multilevel Model

Estimate regression coefficients using multiple chains of Polya-Gamma
Gibbs sampling for multilevel data.

## Usage

``` r
estimate_parameters_ml(
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
  n_chain = 2,
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

  Character vector with names of fixed variables.

- random:

  Character vector with names of random variables.

- n_burn:

  Scalar. Number of burnin iterations.

- n_it:

  Scalar. Number of iterations.

- start:

  Vector of starting values for each chain.

- b_mu0:

  Vector of prior means of fixed regression coefficients.

- b_sigma0:

  Prior covariance matrix of fixed regression coefficients.

- g_mu0:

  Vector of prior means of random regression coefficients.

- g_sigma0:

  Prior covariance matrix of random regression coefficients.

- nu0:

  Scalar. Degrees of freedom of prior inverse-Wishart distribution.

- tau0:

  Prior matrix of inverse-Wishart distribution.

- n_chain:

  Scalar. Number of chains. Default is 2.

- return_thinned:

  Logical. Should thinned chains be returned? Default is FALSE.

- n_thin:

  Thinning rate. Default is 1.

## Value

A named list with:

- Pars:

  List of parameter draws from each chain
