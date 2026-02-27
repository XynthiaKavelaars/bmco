
# bmco

Bayesian methods for comparing groups on multiple binary outcomes.

## Installation

``` r
# From CRAN (when accepted)
install.packages("bmco")

# Development version
# devtools::install_github("XynthiaKavelaars/bmco")
```

## Quick Example

``` r
library(bmco)

# Generate data
set.seed(123)
data <- data.frame(
  treatment = rep(c("control", "drug"), each = 50),
  outcome1 = rbinom(100, 1, 0.5),
  outcome2 = rbinom(100, 1, 0.5)
)

# Analyze
result <- bmvb(
  data = data,
  grp = "treatment",
  grp_a = "control",
  grp_b = "drug",
  y_vars = c("outcome1", "outcome2"),
  n_it = 10000
)

print(result)
```

## Functions

- `bmvb()`: Basic group comparison
- `bglm()`: Subgroup analysis  
- `bglmm()`: Multilevel data

## Getting Help

See `vignette("introduction")` for detailed examples.
