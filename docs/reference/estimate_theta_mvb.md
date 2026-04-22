# Estimate Theta with Multivariate Bernoulli Distribution

Estimate success probabilities with a Multivariate Bernoulli
distribution (reference approach).

## Usage

``` r
estimate_theta_mvb(Y, prior_alpha = 0.5, n_it = 10000)
```

## Arguments

- Y:

  n x K matrix with bivariate binary responses.

- prior_alpha:

  Numeric vector of length Q with prior hyperparameters, or scalar if
  all are equal. Default is 0.5.

- n_it:

  Scalar. Number of draws from posterior distributions. Default is
  10000.

## Value

nIt x K matrix with bivariate Bernoulli probabilities. Currently
supported for K=2 only.
