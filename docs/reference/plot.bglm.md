# Plot Method for bglm Objects

Plot Method for bglm Objects

## Usage

``` r
# S3 method for class 'bglm'
plot(x, type = "all", parameters = NULL, ...)
```

## Arguments

- x:

  A bglm object returned by bglm().

- type:

  Character. Type of plot: "trace", "density", "autocorr", or "all".
  Default is "all".

- parameters:

  Character vector. Which parameters to plot. Default is NULL (all
  parameters).

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns NULL (plots are displayed).

## Examples

``` r
# Uses the pre-computed example object shipped with the package.
# Plot trace plots for the fixed-effect regression coefficients:
plot(bglm_fit, type = "trace")












```
