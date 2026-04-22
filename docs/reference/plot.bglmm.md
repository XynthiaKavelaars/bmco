# Plot Method for bglmm Objects

Plot Method for bglmm Objects

## Usage

``` r
# S3 method for class 'bglmm'
plot(x, type = "all", which = "fixed", parameters = NULL, ...)
```

## Arguments

- x:

  A bglmm object returned by bglmm().

- type:

  Character. Type of plot: "trace", "density", "autocorr", or "all".
  Default is "all".

- which:

  Character. Which component to plot: "fixed" (b), "random" (g),
  "variance" (tau), or "all". Default is "fixed".

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
# Trace plots for the fixed-effect regression coefficients:
plot(bglmm_fit, type = "trace", which = "fixed")










# Trace plots for the random-effect variance components:
plot(bglmm_fit, type = "trace", which = "variance")



```
