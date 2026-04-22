# Subgroup Analysis with Multivariate Binary Outcomes

## Introduction

In many studies, treatment effects may depend on patient characteristics
(e.g., age, disease severity). Ignoring this patient information can,
among others, lead to:

- **Masked treatment effects**: Specific subgroup effects cancel each
  other out
- **Reduced precision**: Unexplained variability
- **Limited generalizability**: Results only apply to the specific
  sample

The
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
function extends
[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md) by
allowing for subgroup analysis while using full sample information
through multinomial logistic regression.

## When to Use Subgroup Analysis

Use
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
when:

- You want to estimate effects for specific subpopulations (e.g.,
  elderly patients, severe cases)
- You want to improve precision by explaining outcome variability
- Treatment effects may vary by covariate level (interaction effects)

Use
[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md)
when:

- No covariates are available or theoretically relevant

## Example: Clinical Trial with Age Effects

We’ll simulate a trial where treatment effectiveness varies by patient
age.

### Generate Data

``` r
library(bmco)
set.seed(2024)

n <- 150
clinical_data <- data.frame(
  treatment = rep(c("standard", "new"), each = n/2),
  age = rnorm(n, mean = 55, sd = 12)
)

# Treatment is more effective in younger patients
clinical_data$pain_relief <- NA
clinical_data$mobility_improved <- NA

for (i in 1:nrow(clinical_data)) {
  is_new <- ifelse(clinical_data$treatment[i] == "new", 1, 0)
  age_centered <- (clinical_data$age[i] - 55) / 10
  
  # Younger patients benefit more from new treatment (interaction effect)
  p_pain <- plogis(-0.5 + 0.8 * is_new - 0.3 * age_centered - 0.4 * is_new * age_centered)
  p_mobility <- plogis(-0.3 + 0.6 * is_new - 0.2 * age_centered - 0.3 * is_new * age_centered)
  
  clinical_data$pain_relief[i] <- rbinom(1, 1, p_pain)
  clinical_data$mobility_improved[i] <- rbinom(1, 1, p_mobility)
}

# View data
head(clinical_data)
#>   treatment      age pain_relief mobility_improved
#> 1  standard 66.78363           0                 1
#> 2  standard 60.62458           0                 1
#> 3  standard 53.70434           1                 0
#> 4  standard 52.44546           1                 1
#> 5  standard 68.89718           0                 1
#> 6  standard 70.50826           0                 0
summary(clinical_data$age)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   15.71   46.48   54.80   54.43   63.16   81.50
table(clinical_data$treatment)
#> 
#>      new standard 
#>       75       75
```

### Full sample analysis

First, analyze treatment effect on full sample:

``` r
result_full_sample <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf), # Full sample
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)

print(result_full_sample)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>     Group mean y1 mean y2
#>  standard    0.41   0.496
#>       new    0.57   0.492
#>   n(standard) = 75    n(new) = 75
#> 
#> Posterior probability P(new > standard) [All rule]: 0.467
#> Marginalization: Empirical over [-Inf, Inf]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
summary(result_full_sample)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#> Group Estimates:
#>     Group mean y1 sd y1 mean y2 sd y2
#>  standard    0.41 0.056   0.496 0.055
#>       new    0.57 0.064   0.492 0.061
#> 
#>   n(standard) = 75    n(new) = 75
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.16, -0.004
#>   Delta SE   (y1, y2): 0.003, 0.003
#>   Posterior probability P(new > standard): 0.467
#> 
#> Regression Coefficients (mean [SD]):
#>                          b11            b10            b01
#> Intercept      1.096 [1.175] -0.417 [1.315] -0.683 [1.226]
#> treatment      2.521 [1.729]  3.145 [1.785]  1.339 [1.848]
#> age           -0.028 [0.022] -0.004 [0.024]  0.009 [0.022]
#> treatment:age -0.030 [0.037] -0.033 [0.038] -0.013 [0.035]
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
#>   Multivariate PSRF (MPSRF): 1.0825
#>             Parameter   ESS   Rhat
#>      b11_Intercept[1] 475.2 1.0041
#>      b11_treatment[2] 314.4 1.0058
#>            b11_age[3] 430.2 1.0058
#>  b11_treatment_age[4] 272.2 1.0001
#>      b10_Intercept[1] 490.9 1.0022
#>      b10_treatment[2] 508.6 1.0086
#>            b10_age[3] 483.1 1.0050
#>  b10_treatment_age[4] 295.9 1.0063
#>      b01_Intercept[1] 475.4 1.0034
#>      b01_treatment[2] 417.3 1.0177
#>            b01_age[3] 583.1 1.0090
#>  b01_treatment_age[4] 364.0 1.0089
```

This gives an average effect across all ages in the sample.

### Subgroup Analysis

Now estimate the treatment difference for people aged 60+:

``` r
result_subgroup <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(60, Inf),
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)
```

``` r
print(result_subgroup)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>     Group mean y1 mean y2
#>  standard   0.350   0.482
#>       new   0.434   0.467
#>   n(standard) = 75    n(new) = 75
#> 
#> Posterior probability P(new > standard) [All rule]: 0.338
#> Marginalization: Empirical over [60, Inf]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
```

This analysis shows the treatment difference for a part of the age
distribution in the sample.

## Three Methods for Population Definition

[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
offers three ways to define the population of interest:

### 1. Value Method: Specific Covariate Level

Estimate effect at a **single covariate value** (e.g., age = 60):

``` r
result_age60 <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 60,  # Effect for 60-year-olds
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Effect for 60-year-olds:\n")
#> Effect for 60-year-olds:
cat("  Mean effect:", result_age60$delta$mean_delta, "\n")
#>   Mean effect: 0.145 -0.008
cat("  Posterior probability:", result_age60$delta$pop, "\n")
#>   Posterior probability: 0.448
```

### 2. Empirical Method: Observed Covariate Range

Average effect over **observed covariate values** in a specified range:

``` r
result_young <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(-Inf, 55),  # Younger patients (age <= 55)
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)

result_old <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(55, Inf),  # Older patients (age > 55)
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)

result_all_ages <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf),  # Older patients (age > 55)
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Effect in younger patients (age <= 55):\n")
#> Effect in younger patients (age <= 55):
cat("  Mean effect:", result_young$delta$mean_delta, "\n")
#>   Mean effect: 0.202 0.002
cat("  Posterior probability:", result_young$delta$pop, "\n\n")
#>   Posterior probability: 0.514

cat("Effect in older patients (age > 55):\n")
#> Effect in older patients (age > 55):
cat("  Mean effect:", result_old$delta$mean_delta, "\n")
#>   Mean effect: 0.107 -0.015
cat("  Posterior probability:", result_old$delta$pop, "\n")
#>   Posterior probability: 0.36

cat("Effect in all patients:\n")
#> Effect in all patients:
cat("  Mean effect:", result_all_ages$delta$mean_delta, "\n")
#>   Mean effect: 0.161 -0.008
cat("  Posterior probability:", result_all_ages$delta$pop, "\n\n")
#>   Posterior probability: 0.46
```

The empirical method uses the **actual distribution** of ages in the
data within each range.

### 3. Analytical Method: Theoretical Covariate Distribution

Average effect assuming covariate follows a **normal distribution**:

``` r
result_analytical <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Analytical",
  x_def = c(40, 70),  # Age range 40-70
  test = "right_sided",
  rule = "All",
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Effect for age 40-70 (analytical integration):\n")
#> Effect for age 40-70 (analytical integration):
cat("  Mean effect:", result_analytical$delta$mean_delta, "\n")
#>   Mean effect: 0.167 -0.006
cat("  Posterior probability:", result_analytical$delta$pop, "\n")
#>   Posterior probability: 0.449
```

The analytical method integrates over a **truncated normal
distribution** fitted to the data.

### Choosing a Method

| Method | When to use | Advantages | Disadvantages |
|----|----|----|----|
| **Value** | Specific subgroup (e.g., typical patient) | Simple, interpretable | Single point |
| **Empirical** | Generalize to observed range | Uses actual data distribution | Limited to observed values, not suitable for small subgroups (low *n*) |
| **Analytical** | Smooth extrapolation | Can extrapolate beyond data, no minimum subgroup size | Assumes normality |

## Comparing Results Across Methods

``` r
cat("Mean treatment difference by method:\n")
#> Mean treatment difference by method:
cat("  Value (age=60):     ", result_age60$delta$mean_delta, "\n")
#>   Value (age=60):      0.145 -0.008
cat("  Empirical (young):  ", result_young$delta$mean_delta, "\n")
#>   Empirical (young):   0.202 0.002
cat("  Empirical (old):    ", result_old$delta$mean_delta, "\n")
#>   Empirical (old):     0.107 -0.015
cat("  Empirical (all):    ", result_all_ages$delta$mean_delta, "\n")
#>   Empirical (all):     0.161 -0.008
cat("  Analytical (40-70): ", result_analytical$delta$mean_delta, "\n")
#>   Analytical (40-70):  0.167 -0.006

cat("Posterior Probabilities by method:\n")
#> Posterior Probabilities by method:
cat("  Value (age=60):     ", result_age60$delta$pop, "\n")
#>   Value (age=60):      0.448
cat("  Empirical (young):  ", result_young$delta$pop, "\n")
#>   Empirical (young):   0.514
cat("  Empirical (old):    ", result_old$delta$pop, "\n")
#>   Empirical (old):     0.36
cat("  Empirical (all):    ", result_all_ages$delta$pop, "\n")
#>   Empirical (all):     0.46
cat("  Analytical (40-70): ", result_analytical$delta$pop, "\n")
#>   Analytical (40-70):  0.449
```

The differences reflect:

1.  **Heterogeneous treatment effects**: More evidence in favor of the
    new treatment in younger patients
2.  **Different target populations**: Each method answers a slightly
    different question
3.  **Sample composition**: Empirical methods reflect the actual age
    distribution

## Understanding Regression Coefficients

The model includes:

- **Intercept**: Baseline log-odds for reference group (standard
  treatment, mean age)
- **Group effect**: Main effect of treatment
- **Age effect**: Effect of age (centered)
- **Interaction**: How treatment effect varies with age

``` r
# Regression parameter estimates
result_all_ages$estimates$b
#>             [,1]        [,2]         [,3] [,4]
#> [1,]  1.09885712 -0.37254242 -0.662861082    0
#> [2,]  2.42652889  3.06187763  1.192947060    0
#> [3,] -0.02765016 -0.00447535  0.008619568    0
#> [4,] -0.03181387 -0.03426866 -0.013661140    0
```

These multinomial regression coefficients have no straightforward
interpretation in terms of treatment effects. Therefore, (posterior
samples of) regression coefficients are transformed to (posterior
samples of) marginal success probabilities ($`\mathbf{\theta}_{A}`$,
$`\mathbf{\theta}_{B}`$) and marginal probability differences
($`\mathbf{\delta} = \mathbf{\theta}_{B} - \mathbf{\theta}_{A}`$) to be
used for further analysis and decision-making via the three
abovementioned methods.

## Decision Rules

All three decision rules work:

### All Rule

``` r
result_all <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 55,
  rule = "All",  # Both outcomes must improve
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Mean delta:", result_all$delta$mean_delta, "\n")
#> Mean delta: 0.174 -0.001
cat("P(new treatment better on BOTH outcomes | age=55):", result_all$delta$pop, "\n")
#> P(new treatment better on BOTH outcomes | age=55): 0.481
```

### Any Rule

``` r
result_any <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 55,
  rule = "Any",  # At least one outcome must improve
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Mean delta:", result_any$delta$mean_delta, "\n")
#> Mean delta: 0.168 -0.006
cat("P(new treatment better on AT LEAST ONE outcome | age=55):", result_any$delta$pop, "\n")
#> P(new treatment better on AT LEAST ONE outcome | age=55): 0.993
```

### Compensatory Rule

``` r
result_comp <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 55,
  rule = "Comp",
  w = c(0.6, 0.4),  # Pain relief weighted more heavily
  n_burn = 200,
  n_it = 500
)
```

``` r
cat("Mean delta:", result_comp$delta$mean_delta, "\n")
#> Mean delta: 0.171 -0.008
cat("P(weighted improvement | age=55):", result_comp$delta$pop, "\n")
#> P(weighted improvement | age=55): 0.961
cat("Weighted effect:", result_comp$delta$w_delta, "\n")
#> Weighted effect: 0.1
```

## Practical Example: Subgroup Analysis

Analyze treatment effects across age quartiles:

``` r
age_quartiles <- quantile(clinical_data$age, probs = c(0, 0.25, 0.5, 0.75, 1))

subgroup_results <- data.frame(
  age_group = c("Q1 (youngest)", "Q2", "Q3", "Q4 (oldest)"),
  age_range = paste0(
    round(age_quartiles[1:4]), "-", 
    round(age_quartiles[2:5])
  ),
  mean_delta1 = NA,
  mean_delta2 = NA,
  post_prob = NA
)

for (i in 1:4) {
  result_q <- bglm(
    data = clinical_data,
    grp = "treatment",
    grp_a = "standard",
    grp_b = "new",
    x_var = "age",
    y_vars = c("pain_relief", "mobility_improved"),
    x_method = "Empirical",
    x_def = c(age_quartiles[i], age_quartiles[i + 1]),
    rule = "All",
    n_burn = 1000,
    n_it = 3000
  )
  subgroup_results$mean_delta1[i] <- result_q$delta$mean_delta[1]
  subgroup_results$mean_delta2[i] <- result_q$delta$mean_delta[2]
  subgroup_results$post_prob[i] <- result_q$delta$pop
}

print(subgroup_results)
#>       age_group age_range mean_delta1 mean_delta2 post_prob
#> 1 Q1 (youngest)     16-46       0.209      -0.009     0.462
#> 2            Q2     46-55       0.195      -0.005     0.469
#> 3            Q3     55-63       0.150      -0.011     0.434
#> 4   Q4 (oldest)     63-82       0.080      -0.014     0.324
```

This reveals how treatment effectiveness varies across age groups.

## Discrete Covariates

Covariates can also be categorical (automatically detected):

``` r

# Add disease severity
clinical_data$severity <- factor(
  sample(c("mild", "severe"), 
         nrow(clinical_data), 
         replace = TRUE)
)

result_discrete <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "severity",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value", # For discrete covariates, use Value method
  x_def = 0,  
  rule = "All",
  n_burn = 200,
  n_it = 500
)

print(result_discrete)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>     Group mean y1 mean y2
#>  standard   0.459   0.593
#>       new   0.550   0.501
#>   n(standard) = 75    n(new) = 75
#> 
#> Posterior probability P(new > standard) [All rule]: 0.151
#> Marginalization: Value over [0]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
```

## Sample Size Considerations

Subgroup analysis typically requires larger samples than
[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md),
especially when the age distribution in the sample is used (method
`Empirical`:

Too small samples may lead to:

- Unstable coefficient estimates
- Wide credible intervals
- Poor convergence

For more information on sample size computations, see @Kavelaars2024.

## Specifying Prior Distributions

Regression coefficients have a **multivariate normal prior**:

``` math
\beta \sim N(\mu_0, \Sigma_0)
```

### Default Priors

``` r
# Default for bglm/bglmm
p <- 4
b_mu0 <- rep(0, p)                 # Zero prior mean
b_sigma0 <- solve(diag(100, p))    # Small prior precision (large prior variance). Weakly informative. 
                                   # solve() generates precision matrix from variance matrix
```

Where `p` is the number of fixed effects:

- Intercept
- Group effect
- Covariate effect(s)
- Interaction term(s)

#### Custom Fixed Effects Priors

**Example 1**: Informative prior favoring positive treatment effect

``` r
# Assume treatment improves outcomes (positive coefficient on group)
# For binary covariate with 4 parameters: Intercept, grp, x, grp:x

custom_b_mu0 <- c(0, 0.5, 0, 0)      # Expect positive group effect
custom_b_sigma0 <- solve(diag(c(10, 5, 10, 10)))  # More certainty on group effect
result_informative <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 55,
  b_mu0 = custom_b_mu0,
  b_sigma0 = custom_b_sigma0,
  n_burn = 200,
  n_it = 500
)
```

**Example 2**: Skeptical prior (no effect expected)

``` r
# Skeptical about treatment effect - tight prior around zero
skeptical_b_mu0 <- c(0, 0, 0, 0)        # No effect expected
skeptical_b_sigma0 <- solve(diag(c(1, 1, 1, 1)))  # Small variance

result_skeptical <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 55,
  b_mu0 = skeptical_b_mu0,
  b_sigma0 = skeptical_b_sigma0,
  n_burn = 200,
  n_it = 500
)
```

### Prior Sensitivity Analysis

Compare results under different priors:

``` r
# Weakly informative (default)
result_weak <- bglm(
  data = clinical_data,
  grp = "treatment", grp_a = "standard", grp_b = "new",
  x_var = "age", y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value", x_def = 55,
  b_mu0 = c(0, 0, 0, 0),
  b_sigma0 = solve(diag(10, 4)),
  test = "right_sided",
  rule = "All",
  n_burn = 200, n_it = 500
)

# More informative
result_more_informative <- bglm(
  data = clinical_data,
  grp = "treatment", grp_a = "standard", grp_b = "new",
  x_var = "age", y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value", x_def = 55,
  b_mu0 = c(0, 0.3, 0, 0.5),
  b_sigma0 = solve(diag(c(1, 5, 1, 1))),
  test = "right_sided",
  rule = "All",
  n_burn = 200, n_it = 500
)
```

``` r
cat("Prior Sensitivity:\n")
#> Prior Sensitivity:
cat("  Weakly informative:\n")
#>   Weakly informative:
cat("    Mean effect:", result_weak$delta$mean_delta, "\n\n")
#>     Mean effect: 0.166 -0.01
cat("    Posterior prob:", result_weak$delta$pop, "\n")
#>     Posterior prob: 0.431

cat("  More informative (favoring treatment):\n")
#>   More informative (favoring treatment):
cat("    Mean effect:", result_more_informative$delta$mean_delta, "\n")
#>     Mean effect: 0.174 -0.006
cat("    Posterior prob:", result_more_informative$delta$pop, "\n")
#>     Posterior prob: 0.468

cat("  Informative (favoring treatment):\n")
#>   Informative (favoring treatment):
cat("    Mean effect:", result_informative$delta$mean_delta, "\n")
#>     Mean effect: 0.172 -0.007
cat("    Posterior prob:", result_informative$delta$pop, "\n")
#>     Posterior prob: 0.454

cat("  Skeptical:\n")
#>   Skeptical:
cat("    Mean effect:", result_skeptical$delta$mean_delta, "\n")
#>     Mean effect: 0.154 -0.014
cat("    Posterior prob:", result_skeptical$delta$pop, "\n")
#>     Posterior prob: 0.428
```

**Interpretation**:

- If results are **similar** → data dominate, prior has little influence
- If results **differ substantially** → prior is influential

### Further Reading

For Bayesian logistic regression priors: - Gelman, A. Jakulin, A.,
Pittau, M. G. & Su, Y-S (2008). A weakly informative default prior
distribution for logistic and other regression models. *Annals of
Applied Statistics, 2*(4), 1360-1383. .

## MCMC Diagnostics

The function internally checks the multivariate potential scale
reduction factor (MPSRF) and warns if \> 1.1. Additional MCMC
diagnostics will be returned When `return_diagnostics = TRUE`. When
`return_diagnostic_plots = TRUE`,
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
returns diagnostic plots as well.

``` r
result <- bglm(
    # ... arguments ...
    n_burn = 200,  # Increase if convergence issues
    n_it = 500,
    return_diagnostics = TRUE,
    return_diagnostic_plots = TRUE
  )
```

Signs of poor convergence:

- Warning: “MCMC chains may not have converged”
- Very wide standard errors
- Unexpected posterior probabilities

Solutions:

1.  Increase `n_burn` (e.g., 5000)
2.  Increase `n_it` (e.g., 10000)
3.  Check for data issues (perfect separation, outliers)

## Comparison: `bmvb()` and `bglm()`

``` r
# Multivariate Bernoulli-distribution:
result_bmvb <- bmvb(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  y_vars = c("pain_relief", "mobility_improved"),
  n_it = 500
)

# Logistic regression, full range:
result_bglm <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Empirical",
  x_def = c(-Inf, Inf),
  n_burn = 200,
  n_it = 500
)

cat("Comparison:\n")
#> Comparison:
cat("  Multivariate Bernoulli:\n") 
#>   Multivariate Bernoulli:
cat(" Mean treatment difference: ", result_bmvb$delta$mean_delta, "\n")
#>  Mean treatment difference:  0.165 0
cat(" Posterior probability: ", result_bmvb$delta$pop, "\n")
#>  Posterior probability:  0.472

cat("  Logistic regression:\n   ")
#>   Logistic regression:
#> 
cat(" Mean treatment difference: ", result_bglm$delta$mean_delta, "\n")
#>  Mean treatment difference:  0.159 -0.008
cat(" Posterior probability: ", result_bglm$delta$pop, "\n")
#>  Posterior probability:  0.45

cat("  Difference:\n            ")
#>   Difference:
#> 
cat(" Treatment difference:", abs(result_bglm$delta$mean_delta - result_bmvb$delta$mean_delta), "\n")
#>  Treatment difference: 0.006 0.008
cat(" Posterior probability:", abs(result_bglm$delta$pop - result_bmvb$delta$pop), "\n")
#>  Posterior probability: 0.022
```

When the full sample range of the covariate is used, results should be
similar.

## Advanced: Extracting Predictions

``` r
result_samples <- bglm(
  data = clinical_data,
  grp = "treatment",
  grp_a = "standard",
  grp_b = "new",
  x_var = "age",
  y_vars = c("pain_relief", "mobility_improved"),
  x_method = "Value",
  x_def = 65,
  n_burn = 200,
  n_it = 500,
  return_samples = TRUE
)
```

``` r
# Treatment effect samples for age = 65
delta_samples <- result_samples$samples$delta

# Custom probability: effect on pain relief > 0.2
cat("P(effect on pain > 0.2 | age=65):", mean(delta_samples[, 1] > 0.2), "\n")
#> P(effect on pain > 0.2 | age=65): 0.193

# Credible interval
ci <- apply(delta_samples, 2, quantile, probs = c(0.025, 0.975))
cat("\n95% Credible Intervals for age=65:\n")
#> 
#> 95% Credible Intervals for age=65:
print(ci)
#>             [,1]       [,2]
#> 2.5%  -0.1035631 -0.2166126
#> 97.5%  0.3162754  0.1985520
```

## Summary

[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
extends Bayesian multivariate analysis by:

1.  **Estimating subgroup effects** (Value, Empirical, Analytical
    methods)
2.  **Handling interactions** between treatment and covariates
3.  **Maintaining multivariate inference** (All, Any, Compensatory
    rules)

Choose your analysis:

- **[`bmvb()`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md)**:
  No covariates → Simple comparison
- **[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)**:
  Covariates → Subgroup comparison
- **[`bglmm()`](https://xynthiakavelaars.github.io/bmco/reference/bglmm.md)**:
  Covariates and/or clustering → Full hierarchical model

## References

For more details on statistical methodology or when using the
[`bglm()`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
function, please cite:

Kavelaars, X., J. Mulder, and M. Kaptein (2024). “Bayesian Multivariate,
Logistic Regression for Superiority and Inferiority Decision-Making,
under Observable Treatment Heterogeneity”. In: *Multivariate Behavioral,
Research* 59.4, pp. 859-882. DOI:,
[10.1080/00273171.2024.2337340](https://doi.org/10.1080%2F00273171.2024.2337340).
