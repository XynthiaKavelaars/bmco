# Sample Beta Coefficients using Polya-Gamma Method

Sample regression coefficients via Polya-Gamma multinomial logistic
regression for a single chain.

## Usage

``` r
sample_beta_pg(X, Y, n_burn, n_it, start, b_mu0, b_sigma0, verbose = FALSE)
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

  Scalar. Starting value for this chain.

- b_mu0:

  Scalar. Prior mean.

- b_sigma0:

  Scalar. Prior precision of regression coefficients.

- verbose:

  Logical. If true, a progress bar is shown.

## Value

b_draw_pg: List of length n_it with P x Q matrices with estimated
regression coefficients
