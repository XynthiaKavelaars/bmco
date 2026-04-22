# Simulated Multilevel Clinical Trial Data

A simulated dataset representing a two-arm clinical trial with 300
subjects nested within 20 clusters (e.g., hospitals), one continuous
covariate, and two binary outcomes. It serves as the underlying data for
[`bglmm_fit`](https://xynthiakavelaars.github.io/bmco/reference/bglmm_fit.md)
and can be used to illustrate
[`bglmm`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md).

## Usage

``` r
bglmm_data
```

## Format

A data frame with 300 rows and 5 columns:

- id:

  Factor with 20 levels (`1`–`20`). Cluster identifier (e.g., hospital).
  Each cluster contains 15 subjects.

- group:

  Character. Treatment arm: `"placebo"` (clusters 1–10) or `"drug"`
  (clusters 11–20).

- age:

  Numeric. Continuous covariate drawn from \\N(50,\\10^2)\\.

- y1:

  Integer (0/1). First binary outcome.

- y2:

  Integer (0/1). Second binary outcome.

## Details

Data were generated with `set.seed(2024)` using logistic models with
cluster-specific random intercepts \\u\_{j1},\\u\_{j2} \sim
N(0,\\0.25)\\: \$\$ P(y_1 = 1) = \text{logit}^{-1}(-0.50 +
0.75\\\text{drug} + 0.10\\\text{age}/10 + u\_{j1}) \$\$ \$\$ P(y_2 = 1)
= \text{logit}^{-1}(-0.50 + 0.80\\\text{drug} + 0.05\\\text{age}/10 +
u\_{j2}) \$\$ where `drug` is 1 for the drug arm and 0 for placebo. See
`data-raw/generate_examples.R` for the full script.

## See also

[`bglmm`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md),
[`bglmm_fit`](https://xynthiakavelaars.github.io/bmco/reference/bglmm_fit.md),
[`bglm_data`](https://xynthiakavelaars.github.io/bmco/reference/bglm_data.md)

## Examples

``` r
head(bglmm_data)
#>   id   group      age y2 y1
#> 1  1 placebo 41.27538  0  1
#> 2  1 placebo 49.95309  1  1
#> 3  1 placebo 54.95375  1  1
#> 4  1 placebo 69.72819  1  1
#> 5  1 placebo 38.03138  1  0
#> 6  1 placebo 51.30438  1  0
table(bglmm_data$group, bglmm_data$id)
#>          
#>            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#>   drug     0  0  0  0  0  0  0  0  0  0 15 15 15 15 15 15 15 15 15 15
#>   placebo 15 15 15 15 15 15 15 15 15 15  0  0  0  0  0  0  0  0  0  0
```
