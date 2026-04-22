# Sample Population Subsets

Extract subsets of covariate data for treatment groups based on
specified method.

## Usage

``` r
sample_population(
  X,
  population_var = NULL,
  grp_var,
  method,
  value = NULL,
  range = NULL,
  fixed = NULL,
  random = NULL
)
```

## Arguments

- X:

  Design matrix with covariate data. Can be (n x P) matrix or list of J
  (n_j x P) matrices for multilevel data.

- population_var:

  Optional character string with the variable name that defines the
  subpopulation. Required when range is not c(-Inf, Inf).

- grp_var:

  Character string with name of variable that defines groups.

- method:

  Character. "Empirical" for empirical marginalization, "Value" for
  vector of fixed values.

- value:

  Scalar. Value of x representing subpopulation. Required when method =
  "Value".

- range:

  Numeric vector of lower and upper bound of covariate that represents
  subpopulation. Required when method = "Empirical".

- fixed:

  Character vector with names of fixed variables in covariate vector.
  Default is NULL.

- random:

  Character vector with names of random variables in covariate vector.
  Default is NULL.

## Value

A list with:

- x_a:

  Design matrix with covariate data for group A

- x_b:

  Design matrix with covariate data for group B
