# Subgroup Analysis with Multivariate Binary Outcomes in Multilevel Data

## Introduction

When data are clustered (e.g., students within schools, patients within
hospitals), observations are not independent. Ignoring this clustering
can lead to:

- Underestimated standard errors
- Inflated Type I error rates
- Invalid inference
- Under- or overpowered decisions

The
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
function extends
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md) to
handle clustered data by incorporating **random effects**.

## When to Use Multilevel Models

Use
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
when:

- Observations are grouped into clusters (schools, hospitals, families,
  etc.);
- Cluster variation affects outcomes
- Clustering is theoretically meaningful

Use
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
when:

- Data are not clustered; AND
- Clustering is not theoretically meaningful

## Example: Educational Intervention Study

We’ll analyze a study comparing two teaching methods across 20 schools
on reading and math ability, differentiating on student baseline
ability.

### Generate Example Data

``` r
library(bmco)
set.seed(2020)

n_schools <- 20 
n_per_school <- 15
n_total <- n_schools * n_per_school

# Generate random intercepts
school_intercepts_math <- rnorm(n_schools, mean = 0, sd = 1.50)
school_intercepts_reading <- rnorm(n_schools, mean = 0, sd = 1.50)

# Cluster randomization: whole schools receive the same method.
school_method <- rep(c("traditional", "new_method"), each = n_schools/2)

study_data <- data.frame(
  school_id = factor(rep(1:n_schools, each = n_per_school)),
  method = rep(school_method, each = n_per_school),
  baseline_ability = rnorm(n_total, mean = 50, sd = 10),
  stringsAsFactors = FALSE
)

# Standardize covariate
study_data$baseline_ability_std <- scale(study_data$baseline_ability)[, 1]

study_data$math_pass <- NA
study_data$reading_pass <- NA

p_math <- p_reading <- rep(NA, nrow(study_data))

for (i in 1:nrow(study_data)) {
  school <- as.numeric(study_data$school_id[i])
  is_new <- ifelse(study_data$method[i] == "new_method", 1, 0)

  p_math[i]    <- plogis(-0.50 + 0.75 * is_new + 0.10 * study_data$baseline_ability_std[i] + 
                           0.20 * is_new * study_data$baseline_ability_std[i] + 
                           school_intercepts_math[school])
  p_reading[i] <- plogis(-0.50 + 0.80 * is_new + 0.05 * study_data$baseline_ability_std[i] + 
                           0.15 * is_new * study_data$baseline_ability_std[i] + 
                           school_intercepts_reading[school])
  
  study_data$math_pass[i]    <- rbinom(1, 1, p_math[i])
  study_data$reading_pass[i] <- rbinom(1, 1, p_reading[i])
}
```

### Fit Multilevel Model

``` r
set.seed(2025)
result <- bglmm(
  data = study_data,
  grp = "method",
  grp_a = "traditional",
  grp_b = "new_method",
  id_var = "school_id",
  x_var = "baseline_ability",
  y_vars = c("math_pass", "reading_pass"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf),
  fixed = c("method", "baseline_ability", "method_baseline_ability"),
  random = c("Intercept"),
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500,
  n_thin = 1
)
```

``` r
print(result)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#>        Group mean y1 mean y2
#>  traditional   0.416   0.546
#>   new_method   0.537   0.634
#>   J = 20 clusters    n(traditional) = 150    n(new_method) = 150
#> 
#> Posterior probability P(new_method > traditional) [All rule]: 0.980
#> Marginalization: Empirical over [-Inf, Inf]
#> MPSRF: fixed = 10.9836    random = 1.3948    variance = 5.4158
#> 
#> Use summary() for full coefficient tables, priors and MCMC diagnostics.
summary(result)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#> Multilevel Structure:  J = 20 clusters    n(traditional) = 150    n(new_method) = 150
#> 
#> Group Estimates:
#>        Group mean y1 sd y1 mean y2 sd y2
#>  traditional   0.416 0.030   0.546 0.033
#>   new_method   0.537 0.031   0.634 0.032
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.121, 0.088
#>   Delta SE   (y1, y2): 0.002, 0.002
#>   Posterior probability P(new_method > traditional): 0.980
#> 
#> Fixed Effects (mean [SD]):
#>                                   b11           b10            b01
#> method                  3.693 [5.041] 2.346 [3.905]  3.596 [4.233]
#> baseline_ability        0.021 [0.032] 0.045 [0.024] -0.005 [0.027]
#> method_baseline_ability 0.233 [0.205] 0.155 [0.194]  0.216 [0.191]
#> 
#> Random Effects (population mean [SD]):
#>                      g11            g10           g01
#> Intercept -0.241 [2.605] -1.000 [2.675] 1.010 [2.475]
#> 
#> Variance Components (posterior mean):
#>   y1=1, y2=1 (b11):
#>           Intercept
#> Intercept   1029.96
#>   y1=1, y2=0 (b10):
#>           Intercept
#> Intercept   544.973
#>   y1=0, y2=1 (b01):
#>           Intercept
#> Intercept   712.445
#> 
#> Prior Specification:
#>   Fixed effects -- Normal prior:
#>     Mean:      0, 0, 0 
#>     Variance:  10, 10, 10 
#>   Random effects -- Normal prior:
#>     Mean:      0 
#>     Variance:  10 
#>   Covariance -- Inverse-Wishart: df = 1
#> 
#> Marginalization:
#>   Method: Empirical    (Sub)population: [-Inf, Inf]
#>   Decision rule: All
#> 
#> MCMC Convergence Diagnostics:
#>   Fixed effects -- MPSRF: 10.9836
#>                       Parameter  ESS   Rhat
#>                   b11_method[1] 36.9 7.0948
#>         b11_baseline_ability[2] 10.7 1.2712
#>  b11_method_baseline_ability[3] 12.9 5.4743
#>                   b10_method[1] 41.2 3.5012
#>         b10_baseline_ability[2] 17.5 1.1125
#>  b10_method_baseline_ability[3] 19.8 5.4723
#>                   b01_method[1] 19.1 2.4033
#>         b01_baseline_ability[2] 11.4 1.6787
#>  b01_method_baseline_ability[3] 10.6 5.5842
#> 
#>   Random effects -- MPSRF: 1.3948
#>         Parameter   ESS   Rhat
#>  g11_Intercept[1] 565.3 1.1999
#>  g10_Intercept[1] 515.1 1.5451
#>  g01_Intercept[1] 512.2 1.2905
#> 
#>   Variance components -- MPSRF: 5.4158
#>               Parameter   ESS   Rhat
#>  g11_InterceptIntercept  89.5 4.4046
#>  g10_InterceptIntercept 106.7 4.9981
#>  g01_InterceptIntercept 305.3 4.5160
```

### Interpretation

The posterior probability indicates the probability that the new
teaching method improves **both** math and reading outcomes (rule =
“All”), accounting for:

1.  **Fixed effects**: Overall method effect, baseline ability effect,
    and their interaction
2.  **Random effects**: School-specific deviations from the overall
    pattern
3.  **Clustering**: Within-school correlation of student outcomes

## Specifying Population of Interest

Like
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md),
you can define subpopulations:

### Specific Ability Level

``` r
set.seed(2027)
result_value <- bglmm(
  # ... other arguments ...
  x_method = "Value",
  x_def = 50,  # Average ability
)
```

``` r
cat("Effect for students with baseline ability = 50:\n")
cat("  Treatment difference:", result_value$delta$mean_delta, "\n")
cat("  Posterior probability:", result_value$delta$pop, "\n")
```

### Ability Range (Empirical)

``` r
set.seed(2028)
result_empirical <- bglmm(
  # ... other arguments ...
  x_method = "Empirical",
  x_def = c(40,60) # Ability range 40-60
)
```

### Ability Range (Analytical)

``` r
set.seed(2029)
result_analytical <- bglmm(
  # ... other arguments ...
    x_method = "Analytical",
    x_def = c(40,60)
 )
```

## Decision Rules with Multilevel Data

### All Rule (Conjunctive)

Both outcomes must favor the new method:

``` r
result_all <- bglmm(
  data = study_data,
  # ... other arguments ...
  rule = "All"  # P(new_method better on BOTH outcomes)
)
```

### Any Rule (Disjunctive)

At least one outcome must favor the new method:

``` r
result_any <- bglmm(
  data = study_data,
  # ... other arguments ...
  rule = "Any"  # P(new_method better on AT LEAST ONE outcome)
)
```

### Compensatory Rule

Weighted combination (e.g., math weighted 60%, reading 40%):

``` r
result_comp <- bglmm(
  data = study_data,
  # ... other arguments ...
  rule = "Comp",
  w = c(0.6, 0.4)  # Weights for math and reading
)
```

## Specifying Prior Distributions

[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
uses **weakly informative priors** by default. You can customize these
priors for:

- Fixed effects regression coefficients (`b_mu0`, `b_sigma0`)
- Random effects means (`g_mu0`, `g_sigma0`)
- Random effects covariance (`nu0`, `tau0`)

### Fixed Effects Priors (bglm and bglmm)

Fixed effects have a **multivariate normal prior**:

``` math
\beta \sim N(\mu_0, \Sigma_0)
```

#### Default Priors

``` r
# Default
p_fixed <- 2
b_mu0 <- rep(0, p_fixed)                 # Zero prior mean
b_sigma0 <- solve(diag(1e1, p_fixed))    # Small prior precision (large prior variance on the log-odds scale). 
                                         # Weakly informative. 
                                         # solve() generates precision matrix needed as input
```

Where `p_fixed` is the number of fixed effects. In our example, we
use: - Covariate effect(s) - Interaction term(s)

#### Custom Fixed Effects Priors

``` r
# Assume treatment improves outcomes (positive coefficient on group)

custom_b_mu0 <- c(0, 0.5)                 # Expect positive group effect 
custom_b_sigma0 <- solve(diag(c(1e1, 5))) # More certainty on group effect
```

### Random Effects Priors

#### Population-Level Random Effects

Random effects means have a **multivariate normal prior**:

``` math
\gamma \sim N(\mu_g, \Sigma_g)
```

**Default**:

``` r
p_random <- 2
g_mu0 <- rep(0, p_random)               # Prior mean
g_sigma0 <- solve(diag(1e1, p_random))  # Weakly informative prior variance; 
                                        # solve() transforms variance matrix to precision matrix needed as input
```

Where `p_random` is the number of fixed effects. In our example, we
use: - Intercept - Group effect

#### Random Effects Covariance

The covariance matrix has an **inverse-Wishart prior**:

``` math
\Sigma \sim \text{InvWishart}(\nu_0, \mathbf{T}_0)
```

**Default**:

``` r
nu0 <- p_random                # Degrees of freedom (dimension)
tau0 <- diag(1e-1, p_random)   # Scale matrix
```

**Important**: - `nu0` must be ≥ `p_random` - `tau0` must be **positive
definite** (never use `diag(0, p_random)`) - Smaller `nu0` = less
informative - Larger `tau0` diagonal values = expect more
between-cluster variation

### Prior Sensitivity Analysis

Compare results under different priors:

``` r
# -----------------------------------------------------------------
# Three priors that pull in clearly different directions.
# -----------------------------------------------------------------

# 1. Weakly informative (neutral starting point)
#    Diffuse on fixed effects, diffuse IW for random effects.
result_weak <- bglmm(
  # ... other arguments ...
  b_mu0    = rep(0, p_fixed),
  b_sigma0 = solve(diag(10, p_fixed)),
  g_mu0    = rep(0, p_random),
  g_sigma0 = solve(diag(10, p_random)),
  nu0  = p_random,
  tau0 = diag(0.1, p_random),
)

# 2. Skeptical prior
#    Tight variance so the prior resists positive or negative evidence. 
#    Also a small tau0, implying homogeneity between schools.
result_skeptical <- bglmm(
  # ... other arguments ...
  b_mu0    = rep(0, p_fixed),               
  b_sigma0 = solve(diag(c(1, 1), p_fixed)),  # Tight — hard(er) to overcome
  g_mu0    = rep(0, p_random),
  g_sigma0 = solve(diag(c(1, 1), p_random)), # Tight — hard(er) to overcome
  nu0  = p_random + 4,                 # More informative IW
  tau0 = diag(0.05, p_random),         # Expects small school variation
)
```

**Interpretation**:

- If results are **similar** → data dominate, prior has little influence
- If results **differ substantially** → prior is influential

### Guidelines for Choosing Priors

1.  **Default priors**: Use for exploratory analysis or when no prior
    information exists

2.  **Informative priors**: Use when:

    - You have strong theoretical expectations
    - Previous studies provide evidence
    - You want to regularize extreme estimates
    - Sample size is very small

3.  **Skeptical priors**: Use for:

    - Conservative hypothesis testing
    - Replication studies (favor null)
    - High-stakes decisions requiring strong evidence

4.  **Prior sensitivity**: Check:

``` r
   # Run with default priors
   # Run with informative priors
   # Compare posterior probabilities and effect sizes
   # If similar → robust; if different → prior-dependent
```

### Common Mistakes to Avoid

- **Don’t**: Use `tau0 = diag(0, p)` (singular matrix)

- **Do**: Use `tau0 = diag(1e-1, p)` or `diag(0.5, p)`

- **Don’t**: Set `nu0 < p` (too few degrees of freedom)

- **Do**: Use \`nu0 =$`\geq p`$

- **Don’t**: Use extremely small prior variances without justification

- **Do**: Use weakly informative priors unless you have strong evidence

- **Don’t**: Ignore prior sensitivity when sample size is small

- **Do**: Report results under multiple priors for transparency

### Further Reading

For Bayesian logistic regression priors: - Gelman, A. Jakulin, A.,
Pittau, M. G. & Su, Y-S (2008). A weakly informative default prior
distribution for logistic and other regression models. *Annals of
Applied Statistics, 2*(4), 1360-1383. .

For inverse-Wishart priors: - Gelman, A. (2006). Prior distributions for
variance parameters in hierarchical models. *Bayesian Analysis, 1*(3),
515-534.

## MCMC Diagnostics

The function internally checks the multivariate potential scale
reduction factor (MPSRF) and warns if \> 1.1. Additional MCMC
diagnostics will be returned When `return_diagnostics = TRUE`. When
`return_diagnostic_plots = TRUE`,
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
returns diagnostic plots as well.

``` r
set.seed(2030)
result_diag <- bglmm(
  data = study_data,
  grp = "method",
  grp_a = "traditional",
  grp_b = "new_method",
  id_var = "school_id",
  x_var = "baseline_ability",
  y_vars = c("math_pass", "reading_pass"),
  x_method = "Value",
  x_def = 50,
  n_burn = 200,
  n_it = 500,
  n_thin = 1,
  start = c(0.5, 1),
  return_diagnostics = TRUE
)
```

Key diagnostics:

- **MPSRF** (Multivariate PSRF): Ideally \< 1.1
- **Effective sample size**: Should be \> 100 per parameter
- **Rhat**: Ideally \< 1.1 for all parameters

``` r
# Check convergence
cat("Random effects convergence:\n")
#> Random effects convergence:
print(result_diag$diags$g$convergence)
#> $mpsrf
#> [1] 6.27427
#> 
#> $psrf
#>                  Point est. Upper C.I.
#> g11_Intercept[1]   1.368093   2.748530
#> g11_method[2]      1.146185   1.151254
#> g10_Intercept[1]   1.135531   1.218255
#> g10_method[2]      2.271145   7.825439
#> g01_Intercept[1]   1.683002   3.983673
#> g01_method[2]      3.646306   8.278269

cat("\nRandom effects covariance convergence:\n")
#> 
#> Random effects covariance convergence:
print(result_diag$diags$tau$convergence)
#> $mpsrf
#> [1] 5.167838
#> 
#> $psrf
#>                        Point est. Upper C.I.
#> g11_InterceptIntercept   4.983014  29.190062
#> g11_methodIntercept      5.274581  21.553906
#> g11_methodmethod         1.201419   1.731316
#> g10_InterceptIntercept   2.539835  13.514819
#> g10_methodIntercept      2.659497  14.379082
#> g10_methodmethod         2.672124  14.473397
#> g01_InterceptIntercept   4.175515  20.831837
#> g01_methodIntercept      4.298845  20.442173
#> g01_methodmethod         4.247503  19.116352
```

Available plots:

- **Autocorrelation**: Is ideally high at lag 0, then drops off quickly
  toward 0 within a few lags.
- **Traceplot**: Chains should fluctuate randomly around a stable mean,
  with no clear trends, drifts, or long flat stretches. Good mixing.
- **Density**: Posterior density of individual parameters (regression
  coefficients and variance parameters) if preferably smooth, unimodal
  (unless theory suggests otherwise), and without irregular spikes.

``` r
plot(result_diag)                                    # All
plot(result_diag, type = "trace")                    # Trace only
plot(result_diag, type = "density")                  # Density only

plot(result_diag, which = "fixed")                   # Fixed effects
plot(result_diag, which = "random")                  # Random effects
plot(result_diag, which = "variance")                # Variance components
plot(result_diag, which = "all")                     # All

# Combine type + which
plot(result_diag, type = "trace", which = "all")          # Trace for all parameters
plot(result_diag, type = "density", which = "variance")   # Variance densities
```

If convergence issues occur:

1.  Increase `n_burn` (more warmup)
2.  Increase `n_it` (more samples)
3.  Check for very few clusters
4.  Check for very small cluster sizes
5.  Check range of covariate and consider standardizing
6.  Check for problems and warnings in data (perfect separation, low
    variation, etc.)

If autocorrelation is high, usually the effective sample size is much
lower than `n_it`:

1.  Consider thinning
2.  Consider increasing `n_it`

For more information on MCMC diagnostics, consult

- Brooks, S. P. & Gelman, A. (1998). General methods for monitoring
  convergence of iterative simulations. *Journal of Computational and
  Graphical Statistics, 7*(4), 434–455.

- Gelman, A. & Rubin, D. B. (1992). Inference from iterative simulation
  using multiple sequences. *Statistical Science, 7*(4), 457–472.

## Extracting Posterior Samples

``` r
set.seed(2031)
result_samples <- bglmm(
# ... other arguments ...
  return_samples = TRUE
)
```

``` r
# Treatment effect samples
str(result_samples$samples$delta)

# Compute custom summaries
cat("Treatment effect quantiles:\n")
print(apply(result_samples$samples$delta, 2, quantile, probs = c(0.025, 0.5, 0.975)))

# Probability that effect on math > 0.1
cat("\nP(delta_math > 0.1):", mean(result_samples$samples$delta[, 1] > 0.1), "\n")
```

## Data Requirements

For reliable multilevel modeling:

1.  **Number of clusters**: Ideally there are sufficient clusters
    - More clusters → Better random effect estimates
2.  **Cluster sizes**:
    - Very small clusters provide little information about random
      effects
    - Unbalanced cluster sizes are acceptable
3.  **Group balance within clusters**: Ideally, each cluster has both
    treatment groups
    - If some clusters have only one group, the model can still fit but
      with a warning
4.  **Missing data**: Handled automatically
    - Rows with missing outcomes or covariates are removed per cluster
    - Complete-case analysis within each cluster

## Common Issues and Solutions

### Warning: “Very few clusters (J \< 5)”

**Problem**: Too few clusters might hinder reliable random effect
estimation

**Solutions**: - Collect more clusters if possible - Use
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
instead (ignores clustering), but only if appropriate

### Warning: “MCMC chains may not have converged”

**Problem**: MPSRF \> 1.1

**Solutions**: 1. Increase burn-in: `n_burn = 2000` or higher 2.
Increase iterations: `n_it = 5000` or higher 3. Reduce thinning:
`n_thin = 1` (keeps more samples) 4. Check data quality (outliers,
sufficient variation)

### Slow computation

**Problem**: Large number of clusters or long chains

**Solutions**: - Start with small `n_it` to verify model runs - Increase
`n_thin` to reduce stored samples: `n_thin = 10` - Run overnight for
production analyses - Typical times: 5-10 minutes for 10 clusters, 10k
iterations

## Summary

[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
extends Bayesian multivariate analysis to clustered data by:

1.  **Accounting for clustering** via random effects
2.  **Allowing for subgroup effects** via logistic regression
3.  **Flexible population definitions** (Value, Empirical, Analytical)
4.  **Multivariate decision rules** (All, Any, Compensatory)

Choose your analysis:

- **[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md)**:
  Full sample, no clustering
- **[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)**:
  Subgroup analysis, no clustering
- **[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)**:
  Subgroup analysis + clustering

All three functions support multivariate binary outcomes and provide
posterior probabilities for treatment effects.

## References

For more details on statistical methodology or when using the
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)-function,
please cite:

Kavelaars, X., J. Mulder, and M. Kaptein (2023). “Bayesian multilevel,
multivariate logistic regression for superiority decision-making under,
observable treatment heterogeneity”. In: *BMC Medical Research,
Methodology* 23.1. DOI:,
[10.1186/s12874-023-02034-z](https://doi.org/10.1186%2Fs12874-023-02034-z).
