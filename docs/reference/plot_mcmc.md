# Plot Specific MCMC Diagnostic

Generate a specific type of diagnostic plot for MCMC chains.

## Usage

``` r
plot_mcmc(chains, type = "trace", parameters = NULL, main = NULL, ...)
```

## Arguments

- chains:

  An mcmc.list object or a fitted model.

- type:

  Character. Type of plot: "trace", "density", or "autocorr".

- parameters:

  Character vector. Which parameters to plot. NULL plots all.

- main:

  Character. Main title for the plot.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns NULL.
