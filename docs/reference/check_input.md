# Check and Validate Input Parameters

Comprehensive input validation for bmvb, bglm, and bglmm functions.
Performs all necessary checks and returns cleaned/validated parameters.

## Usage

``` r
check_input(
  data,
  grp,
  grp_a,
  grp_b,
  y_vars,
  test,
  rule,
  w,
  analysis = c("bmvb", "bglm", "bglmm"),
  prior_a = NULL,
  prior_b = NULL,
  b_mu0 = NULL,
  b_sigma0 = NULL,
  g_mu0 = NULL,
  g_sigma0 = NULL,
  nu0 = NULL,
  tau0 = NULL,
  fixed = NULL,
  random = NULL,
  x_var = NULL,
  x_method = NULL,
  x_def = NULL,
  id_var = NULL,
  n_burn = NULL,
  n_it,
  n_thin = 1
)
```

## Arguments

- data:

  Data frame containing the data.

- grp:

  Character string. Name of the grouping variable.

- grp_a:

  Value of grp indicating first group.

- grp_b:

  Value of grp indicating second group.

- y_vars:

  Character vector. Names of outcome variables.

- test:

  Character. Direction of test ("right_sided" or "left_sided").

- rule:

  Character. Decision rule ("All", "Any", or "Comp").

- w:

  Numeric vector. Weights for compensatory rule (can be NULL).

- analysis:

  Character. Type of analysis: "bmvb", "bglm", or "bglmm".

- prior_a:

  Numeric. Prior for group A (bmvb only). Default is NULL.

- prior_b:

  Numeric. Prior for group B (bmvb only). Default is NULL.

- b_mu0:

  Vector (length = no. of fixed covariates) of prior means of fixed
  regression coefficients (bglm/bglmm only). Default is NULL.

- b_sigma0:

  Prior precision matrix (P_fixed x P_fixed) of fixed regression
  coefficients (bglm/bglmm only). Default is NULL.

- g_mu0:

  Vector (length = no. of random covariates) of prior means of random
  regression coefficients (bglmm only) . Default is NULL.

- g_sigma0:

  Prior precision matrix (P_random x P_random) of random regression
  coefficients (bglmm only). Default is NULL.

- nu0:

  Scalar. Prior df for random covariance matrix (bglmm only). Default is
  NULL.

- tau0:

  Numeric. Prior scale matrix (length(random) x length(random)) for
  random covariance matrix (bglmm only). Default is NULL.

- fixed:

  Character vector. Names of fixed effect variables. Default is c("x",
  "grp_x").

- random:

  Character vector. Names of random effect variables. Default is
  c("Intercept", grp).

- x_var:

  Character string. Name of covariate (bglm/bglmm only). Default is
  NULL.

- x_method:

  Character. Method for handling covariate (bglm/bglmm only). Default is
  NULL.

- x_def:

  Numeric. Defines subpopulation (bglm/bglmm only). Default is NULL.

- id_var:

  Character string. Name of cluster ID (bglmm only). Default is NULL.

- n_burn:

  Integer. Number of burnin iterations.

- n_it:

  Integer. Number of MCMC iterations.

- n_thin:

  Integer. Thinning interval (bglmm only). Default is 1.

## Value

List with validated parameters and cleaned data subsets:

- y_a:

  Matrix of outcomes for group A (with missing data removed)

- y_b:

  Matrix of outcomes for group B (with missing data removed)

- w:

  Weights (generated if NULL and rule = "Comp")

- grp_a:

  Validated group A value

- grp_b:

  Validated group B value

- J:

  Number of clusters (bglmm only)
