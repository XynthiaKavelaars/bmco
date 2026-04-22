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

- [`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md):
  Basic group comparison
- [`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md):
  Subgroup analysis  
- [`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md):
  Multilevel data

## Getting Help

See
[`vignette("introduction")`](https://xynthiakavelaars.github.io/bmco/articles/introduction.md)
for detailed examples.

## Acknowledgements

The statistical underpinnings of this package were developed with
financial support of a NWO (Dutch Research Council) research talent
grant (no. 406.18.505) and the theoretical insights of Maurits Kaptein
(Eindhoven University of Technology, The Netherlands) and Joris Mulder
(Tilburg University, The Netherlands).
