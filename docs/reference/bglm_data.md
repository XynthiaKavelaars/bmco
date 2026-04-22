# Simulated Single-Level Clinical Trial Data

A simulated dataset representing a two-arm clinical trial with 200
subjects, one continuous covariate, and two binary outcomes. It serves
as the underlying data for
[`bglm_fit`](https://xynthiakavelaars.github.io/bmco/reference/bglm_fit.md)
and can be used to illustrate
[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md).

## Usage

``` r
bglm_data
```

## Format

A data frame with 200 rows and 4 columns:

- group:

  Character. Treatment arm: `"placebo"` (n = 100) or `"drug"` (n = 100).

- age:

  Numeric. Continuous covariate drawn from \\N(50,\\10^2)\\.

- y1:

  Integer (0/1). First binary outcome.

- y2:

  Integer (0/1). Second binary outcome.

## Details

Data were generated with `set.seed(2024)` using logistic models: \$\$
P(y_1 = 1) = \text{logit}^{-1}(-0.50 + 0.75\\\text{drug} +
0.10\\\text{age}/10) \$\$ \$\$ P(y_2 = 1) = \text{logit}^{-1}(-0.50 +
0.80\\\text{drug} + 0.05\\\text{age}/10) \$\$ where `drug` is 1 for the
drug arm and 0 for placebo. See `data-raw/generate_examples.R` for the
full script.

## See also

[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md),
[`bglm_fit`](https://xynthiakavelaars.github.io/bmco/reference/bglm_fit.md),
[`bglmm_data`](https://xynthiakavelaars.github.io/bmco/reference/bglmm_data.md)

## Examples

``` r
head(bglm_data)
#>     group      age y2 y1
#> 1 placebo 59.81969  1  1
#> 2 placebo 54.68715  1  0
#> 3 placebo 48.92029  0  0
#> 4 placebo 47.87122  1  1
#> 5 placebo 61.58098  1  0
#> 6 placebo 62.92355  0  1
table(bglm_data$group, bglm_data$y1)
#>          
#>            0  1
#>   drug    34 66
#>   placebo 48 52
```
