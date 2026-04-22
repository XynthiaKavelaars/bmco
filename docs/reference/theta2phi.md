# Transform Theta to Phi

Transform bivariate success probabilities to joint response
probabilities. Currently supported for K=2/Q=4 only.

## Usage

``` r
theta2phi(theta, rho)
```

## Arguments

- theta:

  Numeric vector of bivariate success probabilities.

- rho:

  Scalar pairwise correlation between success probabilities in theta.

## Value

Numeric vector of 4 joint response probabilities, summing to 1 and
ordered as `11`, `10`, `01`, `00`.
