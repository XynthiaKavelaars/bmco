# Transform Phi to Weighted Treatment Difference

Transform joint response probabilities to weighted treatment difference.
Currently supported for Q=4/K=2 only.

## Usage

``` r
phi2delta_w(phi_a, phi_b, weights)
```

## Arguments

- phi_a:

  Numeric vector of Q joint response probabilities for group A, ordered
  as `11`, `10`, `01`, `00`.

- phi_b:

  Numeric vector of Q joint response probabilities for group B, ordered
  as `11`, `10`, `01`, `00`.

- weights:

  Numeric vector of length K with weights for linear combination of
  treatment differences (Compensatory rule).

## Value

Scalar. Weighted treatment difference.
