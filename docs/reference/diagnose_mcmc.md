# Diagnose MCMC Chains

Perform diagnostic checks on MCMC chains including convergence
assessment and optional trace plots.

## Usage

``` r
diagnose_mcmc(chains, param_names = NULL)
```

## Arguments

- chains:

  An mcmc.list object or a list of mcmc.list objects containing the MCMC
  chains to diagnose.

- param_names:

  Optional character vector of parameter names for plot labels. If NULL,
  uses column names from chains.

## Value

A list with diagnostic information:

- convergence:

  Multivariate potential scale reduction factor (Gelman-Rubin statistic)

- n_eff:

  Effective sample sizes for each parameter

- rhat:

  Univariate potential scale reduction factors for each parameter

- summary:

  Posterior summary of MCMC chains
