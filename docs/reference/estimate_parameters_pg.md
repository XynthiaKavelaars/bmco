# Estimate Parameters using Polya-Gamma Method

Wrapper function for MCMC procedure to estimate regression coefficients
using Polya-Gamma Gibbs sampling.

## Usage

``` r
estimate_parameters_pg(X, Y, n_burn, n_it, start, b_mu0, b_sigma0, n_chain)
```

## Arguments

- X:

  (n x P) design matrix.

- Y:

  (n x Q) matrix of multinomial response data.

- n_burn:

  Scalar. Number of burnin iterations.

- n_it:

  Scalar. Number of iterations.

- start:

  Vector of two starting values, each used for a different chain.

- b_mu0:

  Scalar or vector (length(P)) of prior mean(s).

- b_sigma0:

  Scalar or matrix (P x P) of prior variance of regression coefficients.

- n_chain:

  Scalar. Number of chains to be sampled.

## Value

A list of length n_chain, being a list of length n_it with P x Q
matrices with estimated regression coefficients
