# Estimate Theta via Empirical Marginalization

Estimate success probabilities via empirical marginalization over
covariate values.

## Usage

``` r
estimate_theta_empirical(est_pars, X)
```

## Arguments

- est_pars:

  List of nIt (P x Q) matrices/arrays of posterior draws of regression
  coefficients.

- X:

  Design matrix with covariate data.

## Value

List of nIt vectors of length K containing bivariate probabilities.
Currently supported for K=2 only.
