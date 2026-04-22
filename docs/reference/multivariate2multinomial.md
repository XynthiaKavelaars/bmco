# Transform Bivariate Binomial to Multinomial Response

Transform bivariate binomial response data to multinomial response
format.

## Usage

``` r
multivariate2multinomial(y_bv)
```

## Arguments

- y_bv:

  n x 2 matrix with bivariate binomial responses.

## Value

n x 4 matrix with multinomial responses, ordered as `11`, `10`, `01`,
`00`.
