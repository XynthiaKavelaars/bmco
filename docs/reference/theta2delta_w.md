# Transform Theta to Weighted Treatment Difference

Transform bivariate success probabilities to weighted treatment
difference. Currently supported for K=2 only.

## Usage

``` r
theta2delta_w(theta_a, theta_b, weights)
```

## Arguments

- theta_a:

  Numeric vector of bivariate success probabilities for group A.

- theta_b:

  Numeric vector of bivariate success probabilities for group B.

- weights:

  Numeric vector of length K with weights for linear combination of
  treatment differences (Compensatory rule).

## Value

Scalar. Weighted treatment difference.
