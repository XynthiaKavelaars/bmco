# Extract MCMC Chains from Model

Extract MCMC chains from a fitted bglm or bglmm model for custom
plotting or analysis.

## Usage

``` r
extract_chains(model, component = "fixed")
```

## Arguments

- model:

  A bglm or bglmm object.

- component:

  Character. For bglmm: "fixed" (b), "random" (g), or "variance" (tau).
  For bglm: only "fixed" is available. Default is "fixed".

## Value

An mcmc.list object containing the MCMC chains.
