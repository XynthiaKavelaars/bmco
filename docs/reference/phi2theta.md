# Transform Phi to Theta

Transform joint response probabilities to bivariate success
probabilities. Currently supported for K=2/Q=4 only.

## Usage

``` r
phi2theta(phi)
```

## Arguments

- phi:

  Numeric vector of 4 joint response probabilities, summing to 1 and
  ordered as `11`, `10`, `01`, `00`.

## Value

Numeric vector of bivariate success probabilities.
