# Introduction to bmco

## Overview

The **bmco** package provides Bayesian methods for analyzing multiple
binary outcomes simultaneously. It addresses a common problem in
clinical trials, educational research, and other fields: how to compare
groups when success is defined by multiple endpoints.

### The Problem

Suppose you’re evaluating a new drug and want to know if it improves
**both** symptom relief and quality of life.

`bmco` provides a principled Bayesian framework that:

- Analyzes **multiple binary outcomes jointly**
- Accounts for **correlations between outcomes**
- Offers **flexible decision rules** (All, Any, or weighted
  combinations)
- Provides **posterior probabilities**
- Estimates **subgroup effects** using data from the entire sample
- Handles **clustered data**

## Quick Start Example

``` r
library(bmco)
set.seed(2024)

# Simulate a clinical trial with two binary outcomes
trial_data <- data.frame(
  treatment = rep(c("placebo", "drug"), each = 50),
  symptom_relief = rbinom(100, 1, prob = rep(c(0.4, 0.6), each = 50)),
  quality_of_life = rbinom(100, 1, prob = rep(c(0.5, 0.7), each = 50))
)

# Test if drug improves BOTH outcomes
result <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  test = "right_sided",
  rule = "All",
  n_it = 1000
)

print(result)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.520   0.481
#>     drug   0.616   0.639
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.795
#> 
#> Use summary() for credible intervals and ESS.
summary(result)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#> Group Estimates:
#>    Group mean y1 sd y1 mean y2 sd y2
#>  placebo   0.520 0.067   0.481 0.072
#>     drug   0.616 0.068   0.639 0.067
#> 
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.096, 0.157
#>   Delta SE   (y1, y2): 0.003, 0.003
#>   Posterior probability P(drug > placebo): 0.795
#> 
#> Test Information:
#>   Decision rule: All
#>   Hypothesis: P(drug > placebo)
#> 
#>   (Run bmvb() with return_samples = TRUE for credible intervals and ESS.)
```

The posterior probability (`pop`) of this example tells us: **“What is
the probability that the drug is better than placebo on BOTH
outcomes?”**

## Three Analysis Functions

The package provides three main functions for different scenarios:

### 1. `bmvb()`: Basic Comparison

**Use when**: Comparing two groups on multiple binary outcomes.
Currently supported for two outcomes only.

``` r
# Compare two teaching methods on two outcomes
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  rule = "All",
  n_it = 1000
)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.521   0.482
#>     drug   0.616   0.631
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.791
#> 
#> Use summary() for credible intervals and ESS.
```

For more details on statistical methodology or when using this function,
please cite:

\[1\] X. Kavelaars, J. Mulder, and M. Kaptein. “Decision-making with,
multiple correlated binary outcomes in clinical trials”. In:,
*Statistical Methods in Medical Research* 29.11 (2020), pp. 3265-3277.,
DOI: 10.1177/0962280220922256.

### 2. `bglm()`: Subgroup analysis

**Use when**: Comparing two subgroups on multiple binary outcomes (e.g.,
age range, particular baseline severity). Currently supported for two
outcomes only.

``` r
# Add age covariate
trial_data$age <- rnorm(100, mean = 50, sd = 10)

result_bglm <- bglm(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  x_var = "age",
  y_vars = c("symptom_relief", "quality_of_life"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf),
  rule = "All",
  n_burn = 200,
  n_it = 500
)

print(result_bglm)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.518   0.479
#>     drug   0.637   0.653
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.859
#> Marginalization: Empirical over [-Inf, Inf]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
summary(result_bglm) 
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#> Group Estimates:
#>    Group mean y1 sd y1 mean y2 sd y2
#>  placebo   0.518 0.069   0.479 0.068
#>     drug   0.637 0.077   0.653 0.073
#> 
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.119, 0.173
#>   Delta SE   (y1, y2): 0.003, 0.003
#>   Posterior probability P(drug > placebo): 0.859
#> 
#> Regression Coefficients (mean [SD]):
#>                          b11            b10            b01
#> Intercept     -0.731 [1.667] -0.093 [1.727] -0.502 [1.771]
#> treatment     -0.717 [2.296] -1.080 [2.280] -0.290 [2.350]
#> age            0.014 [0.033] -0.005 [0.034] -0.001 [0.035]
#> treatment:age  0.119 [0.226]  0.124 [0.224]  0.114 [0.226]
#> 
#> Prior Specification (regression coefficients):
#>   Mean:
#>   Intercept treatment age treatment:age
#> 1         0         0   0             0
#>   Variance (diagonal of inverse precision):
#>   Intercept treatment age treatment:age
#> 1        10        10  10            10
#> 
#> Marginalization:
#>   Method: Empirical
#>   (Sub)population: [-Inf, Inf]
#>   Decision rule: All
#> 
#> MCMC Diagnostics (regression coefficients):
#>   Multivariate PSRF (MPSRF): 1.0352
#>             Parameter   ESS   Rhat
#>      b11_Intercept[1] 460.8 0.9985
#>      b11_treatment[2] 393.8 1.0084
#>            b11_age[3] 445.6 0.9997
#>  b11_treatment_age[4] 199.7 1.0040
#>      b10_Intercept[1] 502.1 1.0002
#>      b10_treatment[2] 363.6 1.0004
#>            b10_age[3] 486.9 0.9992
#>  b10_treatment_age[4] 163.7 0.9997
#>      b01_Intercept[1] 394.3 1.0015
#>      b01_treatment[2] 326.3 1.0176
#>            b01_age[3] 440.9 1.0021
#>  b01_treatment_age[4] 260.4 1.0140
```

For more details on statistical methodology or when using this function,
please cite:

\[1\] X. Kavelaars, J. Mulder, and M. Kaptein. “Bayesian Multivariate,
Logistic Regression for Superiority and Inferiority Decision-Making,
under Observable Treatment Heterogeneity”. In: *Multivariate Behavioral,
Research* 59.4 (2024), pp. 859-882. DOI: 10.1080/00273171.2024.2337340.

Also see the “Subgroup analysis” vignette for details.

### 3. `bglmm()`: With Clustering

**Use when**: Comparing two (sub)groups on multiple binary outcomes,
when data are clustered (patients within hospitals, students within
schools).

``` r
# Add cluster variable
trial_data$hospital <- rep(1:10, each = 10)

result_bglmm <- suppressWarnings( # Warnings due to small number of iterations suppressed
  bglmm(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  id_var = "hospital",
  x_var = "age",
  y_vars = c("symptom_relief", "quality_of_life"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf),
  rule = "All",
  n_burn = 200,
  n_it = 500,
  n_thin = 3
)
)

print(result_bglmm)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.539   0.487
#>     drug   0.623   0.634
#>   J = 10 clusters    n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.765
#> Marginalization: Empirical over [-Inf, Inf]
#> MPSRF: fixed = 4.7581    random = 7.3981    variance = 1.2611
#> 
#> Use summary() for full coefficient tables, priors and MCMC diagnostics.
summary(result_bglmm)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#> Multilevel Structure:  J = 10 clusters    n(placebo) = 50    n(drug) = 50
#> 
#> Group Estimates:
#>    Group mean y1 sd y1 mean y2 sd y2
#>  placebo   0.539 0.069   0.487 0.070
#>     drug   0.623 0.062   0.634 0.059
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.084, 0.147
#>   Delta SE   (y1, y2): 0.007, 0.007
#>   Posterior probability P(drug > placebo): 0.765
#> 
#> Fixed Effects (mean [SD]):
#>                         b11           b10           b01
#> age           0.008 [0.022] 0.014 [0.030] 0.034 [0.023]
#> treatment_age 0.118 [0.202] 0.077 [0.206] 0.070 [0.225]
#> 
#> Random Effects (population mean [SD]):
#>                      g11            g10            g01
#> Intercept -0.365 [1.074] -1.025 [1.497] -2.321 [1.154]
#> treatment -0.854 [1.385]  1.226 [1.021]  1.910 [1.576]
#> 
#> Variance Components (posterior mean):
#>   y1=1, y2=1 (b11):
#>           Intercept treatment
#> Intercept     0.133    -0.008
#> treatment    -0.008     0.225
#>   y1=1, y2=0 (b10):
#>           Intercept treatment
#> Intercept     0.182    -0.075
#> treatment    -0.075     0.244
#>   y1=0, y2=1 (b01):
#>           Intercept treatment
#> Intercept     0.208    -0.083
#> treatment    -0.083     0.211
#> 
#> Prior Specification:
#>   Fixed effects -- Normal prior:
#>     Mean:      0, 0 
#>     Variance:  10, 10 
#>   Random effects -- Normal prior:
#>     Mean:      0, 0 
#>     Variance:  10, 10 
#>   Covariance -- Inverse-Wishart: df = 2
#> 
#> Marginalization:
#>   Method: Empirical    (Sub)population: [-Inf, Inf]
#>   Decision rule: All
#> 
#> MCMC Convergence Diagnostics:
#>   Fixed effects -- MPSRF: 4.7581
#>             Parameter  ESS   Rhat
#>            b11_age[1] 21.2 4.6628
#>  b11_treatment_age[2]  6.2 1.7603
#>            b10_age[1] 23.0 5.2521
#>  b10_treatment_age[2] 17.1 1.1535
#>            b01_age[1] 12.2 1.8520
#>  b01_treatment_age[2]  6.5 1.0558
#> 
#>   Random effects -- MPSRF: 7.3981
#>         Parameter  ESS   Rhat
#>  g11_Intercept[1] 14.3 5.2441
#>  g11_treatment[2]  4.5 1.9938
#>  g10_Intercept[1]  9.2 5.8404
#>  g10_treatment[2] 11.2 1.2483
#>  g01_Intercept[1]  6.9 1.9428
#>  g01_treatment[2]  5.5 1.0851
#> 
#>   Variance components -- MPSRF: 1.2611
#>               Parameter   ESS   Rhat
#>  g11_InterceptIntercept 164.4 1.2415
#>  g11_treatmentIntercept 106.1 1.4007
#>  g11_treatmenttreatment 113.9 1.0590
#>  g10_InterceptIntercept 136.3 1.1900
#>  g10_treatmentIntercept 124.5 0.9952
#>  g10_treatmenttreatment 108.0 1.1740
#>  g01_InterceptIntercept  74.9 1.0785
#>  g01_treatmentIntercept  49.6 1.1052
#>  g01_treatmenttreatment  85.3 1.1856
```

For more details on statistical methodology or when using this function,
please cite:

\[1\] X. Kavelaars, J. Mulder, and M. Kaptein. “Bayesian multilevel,
multivariate logistic regression for superiority decision-making under,
observable treatment heterogeneity”. In: *BMC Medical Research,
Methodology* 23.1 (2023). DOI: 10.1186/s12874-023-02034-z.

Also see the “Multilevel Models” vignette for details.

## Decision Rules

Choose how to define “treatment success”:

### All Rule (Conjunctive)

Treatment must improve **ALL** outcomes:

``` r
result_all <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  rule = "All",  # BOTH outcomes must improve
  n_it = 5000
)

cat("P(drug better on BOTH outcomes):", result_all$delta$pop, "\n")
#> P(drug better on BOTH outcomes): 0.793
```

### Any Rule (Disjunctive)

Treatment must improve **AT LEAST ONE** outcome:

``` r
result_any <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  rule = "Any",  # At least one outcome must improve
  n_it = 5000
)

cat("P(drug better on AT LEAST ONE outcome):", result_any$delta$pop, "\n")
#> P(drug better on AT LEAST ONE outcome): 0.992
```

### Compensatory Rule (Weighted)

Weight outcomes by importance:

``` r
result_comp <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  rule = "Comp",
  w = c(0.7, 0.3),  # Symptom relief weighted 70%, QoL 30%
  n_it = 5000
)

cat("P(weighted improvement):", result_comp$delta$pop, "\n")
#> P(weighted improvement): 0.936
cat("Weighted effect:", result_comp$delta$w_delta, "\n")
#> Weighted effect: 0.114
```

## Test Directions

### Right-sided Test

Test if group B is better than group A:

``` r
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",      
  grp_b = "drug", # Put the group you expect to be better here
  y_vars = c("symptom_relief", "quality_of_life"),
  test = "right_sided",  # P(drug > placebo)
  n_it = 5000
)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.519   0.481
#>     drug   0.617   0.635
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.798
#> 
#> Use summary() for credible intervals and ESS.
```

### Left-sided Test

Test if group A is better than group B:

``` r
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  test = "left_sided",  # P(placebo > drug)
  n_it = 5000
)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.521   0.481
#>     drug   0.615   0.634
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(placebo > drug) [All rule]: 0.012
#> 
#> Use summary() for credible intervals and ESS.
```

**Note**: Left-sided with (A, B) is equivalent to right-sided with (B,
A).

## Understanding the Output

``` r
result <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  n_it = 5000
)

# Key components
names(result)
#> [1] "estimates"    "sample_sizes" "delta"        "info"

# Posterior estimates for each group
result$estimates
#> $mean_a
#> [1] 0.520 0.481
#> 
#> $mean_b
#> [1] 0.615 0.635
#> 
#> $sd_a
#> [1] 0.069 0.069
#> 
#> $sd_b
#> [1] 0.066 0.066

# Sample sizes
result$sample_sizes
#> $n_a
#> [1] 50
#> 
#> $n_b
#> [1] 50

# Treatment effect
result$delta
#> $mean_delta
#> [1] 0.095 0.155
#> 
#> $se_delta
#> [1] 0.001 0.001
#> 
#> $pop
#> [1] 0.795

# Test information
result$info
#> $rule
#> [1] "All"
#> 
#> $test
#> [1] "right_sided"
#> 
#> $test_label
#> [1] "P(drug > placebo)"
#> 
#> $grp_a
#> [1] "placebo"
#> 
#> $grp_b
#> [1] "drug"
```

### Key Output Elements

- `estimates$mean_a` / `mean_b`: Posterior mean success probabilities
  for each group
- `estimates$sd_a` / `sd_b`: Posterior standard deviations
- `delta$mean_delta`: Mean treatment effect (difference in
  probabilities)
- `delta$se_delta`: Standard error of treatment effect
- `delta$pop`: **Posterior probability** (the main result!)
- `info$test_label`: Human-readable description of the test

## Posterior Samples

For custom analyses, extract posterior samples:

``` r
result_samples <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  n_it = 5000,
  return_samples = TRUE
)

# Access posterior samples
str(result_samples$samples)
#> List of 3
#>  $ theta_a: num [1:5000, 1:2] 0.379 0.541 0.596 0.592 0.468 ...
#>  $ theta_b: num [1:5000, 1:2] 0.743 0.514 0.722 0.634 0.696 ...
#>  $ delta  : num [1:5000, 1:2] 0.3641 -0.0268 0.1256 0.0416 0.2285 ...

# Custom probability calculations
samples <- result_samples$samples$delta

# P(effect on symptom relief > 0.1)
cat("P(delta[symptom] > 0.1):", mean(samples[, 1] > 0.1), "\n")
#> P(delta[symptom] > 0.1): 0.4832

# P(effect on QoL > 0.15)
cat("P(delta[QoL] > 0.15):", mean(samples[, 2] > 0.15), "\n")
#> P(delta[QoL] > 0.15): 0.5162

# Credible intervals
cat("\n95% Credible Intervals:\n")
#> 
#> 95% Credible Intervals:
print(apply(samples, 2, quantile, probs = c(0.025, 0.975)))
#>              [,1]        [,2]
#> 2.5%  -0.09942862 -0.03655631
#> 97.5%  0.28142777  0.34018042
```

## Practical Considerations

### Sample Size

For basic comparison: see @Kavelaars2020. For subgroup analysis: see
@Kavelaars2024.

### MCMC Settings

Default settings usually work, but you can adjust:

``` r
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  n_it = 10000  # More iterations for smoother estimates
)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.520   0.481
#>     drug   0.616   0.634
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.791
#> 
#> Use summary() for credible intervals and ESS.
```

For
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
and
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md):

``` r
bglm(
  # ... arguments ...
  n_burn = 2000,  # Burn-in / warm-up iterations
  n_it = 5000,    # Sampling iterations
  n_thin = 1      # Keep every nth sample (reduces memory)
)
```

### Missing Data

Missing values are automatically handled via complete-case analysis:

``` r
# Introduce missing data
trial_data_missing <- trial_data
trial_data_missing$symptom_relief[1:5] <- NA

result_missing <- bmvb(
  data = trial_data_missing,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  n_it = 1000
)
#> Warning in check_input(data = data, grp = grp, grp_a = grp_a, grp_b = grp_b, :
#> Missing data detected in group placebo (5 rows). Rows with missing data will be
#> removed.

# Sample sizes are reduced
result_missing$sample_sizes
#> $n_a
#> [1] 45
#> 
#> $n_b
#> [1] 50
```

## Comparison with Frequentist Approaches

### Traditional approach

``` r
# Separate tests for each outcome
chisq_symptom <- chisq.test(table(trial_data$treatment, trial_data$symptom_relief))
chisq_qol <- chisq.test(table(trial_data$treatment, trial_data$quality_of_life))

cat("Frequentist p-values:\n")
#> Frequentist p-values:
cat("  Symptom relief:", chisq_symptom$p.value, "\n")
#>   Symptom relief: 0.4191152
cat("  Quality of life:", chisq_qol$p.value, "\n")
#>   Quality of life: 0.1584835
cat("  Bonferroni-adjusted alpha: 0.025 (needed for Any-rule only; not needed for All-rule)\n\n")
#>   Bonferroni-adjusted alpha: 0.025 (needed for Any-rule only; not needed for All-rule)
```

### Bayesian approach

``` r
result_bayes <- bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  rule = "All",
  n_it = 1000
)

cat("Bayesian result:\n")
#> Bayesian result:
cat("  P(drug better on BOTH):", result_bayes$delta$pop, "\n")
#>   P(drug better on BOTH): 0.787
```

### Advantages of Bayesian approach

1.  **Direct probability statements**: “95% probability drug is better”
    vs. “reject null at alpha=0.05”
2.  **Flexible decision rules**: All/Any/Weighted combinations
3.  **Coherent uncertainty**: Full sample of posterior distribution
    available

## Next Steps

Explore the other vignettes:

- **Subgroup analysis**: Using
  [`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
  to perform subgroup analysis
- **Multilevel Models**: Using
  [`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
  for clustered data

## References

**Statistical methodology**: More details about statistical methodology
can be found here in the papers below.

When using the
[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md)
function, please cite: Kavelaars, X., J. Mulder, and M. Kaptein (2020).
“Decision-making with, multiple correlated binary outcomes in clinical
trials”. In:, *Statistical Methods in Medical Research* 29.11,
pp. 3265-3277. DOI:,
[10.1177/0962280220922256](https://doi.org/10.1177%2F0962280220922256).

When using the
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
function, please cite: Kavelaars, X., J. Mulder, and M. Kaptein (2024).
“Bayesian Multivariate, Logistic Regression for Superiority and
Inferiority Decision-Making, under Observable Treatment Heterogeneity”.
In: *Multivariate Behavioral, Research* 59.4, pp. 859-882. DOI:,
[10.1080/00273171.2024.2337340](https://doi.org/10.1080%2F00273171.2024.2337340).

When using the
[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)
function, please cite: Kavelaars, X., J. Mulder, and M. Kaptein (2023).
“Bayesian multilevel, multivariate logistic regression for superiority
decision-making under, observable treatment heterogeneity”. In: *BMC
Medical Research, Methodology* 23.1. DOI:,
[10.1186/s12874-023-02034-z](https://doi.org/10.1186%2Fs12874-023-02034-z).
