## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----echo=FALSE---------------------------------------------------------------
biblio <- RefManageR::ReadBib("../inst/REFERENCES.bib", check = "error")

## ----setup, eval=TRUE---------------------------------------------------------
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
summary(clinical_data$age)
table(clinical_data$treatment)

## ----full_sample, eval=TRUE---------------------------------------------------
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
summary(result_full_sample)

## ----adjusted_analysis, results='hide', message=FALSE, eval=TRUE--------------
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

## ----show_adjusted, eval=TRUE-------------------------------------------------
print(result_subgroup)

## ----value_method, results='hide', message=FALSE, eval=TRUE-------------------
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

## ----show_value, eval=TRUE----------------------------------------------------
cat("Effect for 60-year-olds:\n")
cat("  Mean effect:", result_age60$delta$mean_delta, "\n")
cat("  Posterior probability:", result_age60$delta$pop, "\n")


## ----empirical_method, results='hide', message=FALSE, eval=TRUE---------------
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

## ----show_empirical, eval=TRUE------------------------------------------------
cat("Effect in younger patients (age <= 55):\n")
cat("  Mean effect:", result_young$delta$mean_delta, "\n")
cat("  Posterior probability:", result_young$delta$pop, "\n\n")

cat("Effect in older patients (age > 55):\n")
cat("  Mean effect:", result_old$delta$mean_delta, "\n")
cat("  Posterior probability:", result_old$delta$pop, "\n")

cat("Effect in all patients:\n")
cat("  Mean effect:", result_all_ages$delta$mean_delta, "\n")
cat("  Posterior probability:", result_all_ages$delta$pop, "\n\n")


## ----analytical_method, results='hide', message=FALSE, eval=TRUE--------------
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

## ----show_analytical, eval=TRUE-----------------------------------------------
cat("Effect for age 40-70 (analytical integration):\n")
cat("  Mean effect:", result_analytical$delta$mean_delta, "\n")
cat("  Posterior probability:", result_analytical$delta$pop, "\n")

## ----compare_methods,eval=TRUE------------------------------------------------
cat("Mean treatment difference by method:\n")
cat("  Value (age=60):     ", result_age60$delta$mean_delta, "\n")
cat("  Empirical (young):  ", result_young$delta$mean_delta, "\n")
cat("  Empirical (old):    ", result_old$delta$mean_delta, "\n")
cat("  Empirical (all):    ", result_all_ages$delta$mean_delta, "\n")
cat("  Analytical (40-70): ", result_analytical$delta$mean_delta, "\n")

cat("Posterior Probabilities by method:\n")
cat("  Value (age=60):     ", result_age60$delta$pop, "\n")
cat("  Empirical (young):  ", result_young$delta$pop, "\n")
cat("  Empirical (old):    ", result_old$delta$pop, "\n")
cat("  Empirical (all):    ", result_all_ages$delta$pop, "\n")
cat("  Analytical (40-70): ", result_analytical$delta$pop, "\n")

## ----show_coefficients, eval=TRUE---------------------------------------------
# Regression parameter estimates
result_all_ages$estimates$b

## ----rule_all, results='hide', message=FALSE,eval=TRUE------------------------
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

## ----show_all, eval=TRUE------------------------------------------------------
cat("Mean delta:", result_all$delta$mean_delta, "\n")
cat("P(new treatment better on BOTH outcomes | age=55):", result_all$delta$pop, "\n")

## ----rule_any, results='hide', message=FALSE, eval=TRUE-----------------------
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

## ----show_any, eval=TRUE------------------------------------------------------
cat("Mean delta:", result_any$delta$mean_delta, "\n")
cat("P(new treatment better on AT LEAST ONE outcome | age=55):", result_any$delta$pop, "\n")

## ----rule_comp, results='hide', message=FALSE, eval=TRUE----------------------
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

## ----show_comp, eval=TRUE-----------------------------------------------------
cat("Mean delta:", result_comp$delta$mean_delta, "\n")
cat("P(weighted improvement | age=55):", result_comp$delta$pop, "\n")
cat("Weighted effect:", result_comp$delta$w_delta, "\n")

## ----subgroup_analysis, eval=TRUE---------------------------------------------
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

## ----discrete_covariate, eval=TRUE--------------------------------------------

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

## ----default_priors, eval=TRUE------------------------------------------------
# Default for bglm/bglmm
p <- 4
b_mu0 <- rep(0, p)                 # Zero prior mean
b_sigma0 <- solve(diag(100, p))    # Small prior precision (large prior variance). Weakly informative. 
                                   # solve() generates precision matrix from variance matrix

## ----custom_prior_1, results='hide', message=FALSE, warning=FALSE, eval=TRUE----
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

## ----custom_prior_2, results='hide', message=FALSE, warning=FALSE, eval=TRUE----
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

## ----prior_sensitivity, results='hide', message=FALSE, warning=FALSE, eval=TRUE----
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

## ----compare_priors, eval=TRUE------------------------------------------------
cat("Prior Sensitivity:\n")
cat("  Weakly informative:\n")
cat("    Mean effect:", result_weak$delta$mean_delta, "\n\n")
cat("    Posterior prob:", result_weak$delta$pop, "\n")

cat("  More informative (favoring treatment):\n")
cat("    Mean effect:", result_more_informative$delta$mean_delta, "\n")
cat("    Posterior prob:", result_more_informative$delta$pop, "\n")

cat("  Informative (favoring treatment):\n")
cat("    Mean effect:", result_informative$delta$mean_delta, "\n")
cat("    Posterior prob:", result_informative$delta$pop, "\n")

cat("  Skeptical:\n")
cat("    Mean effect:", result_skeptical$delta$mean_delta, "\n")
cat("    Posterior prob:", result_skeptical$delta$pop, "\n")

## ----convergence, eval = FALSE------------------------------------------------
# result <- bglm(
#     # ... arguments ...
#     n_burn = 200,  # Increase if convergence issues
#     n_it = 500,
#     return_diagnostics = TRUE,
#     return_diagnostic_plots = TRUE
#   )

## ----comparison, eval=TRUE----------------------------------------------------
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
cat("  Multivariate Bernoulli:\n") 
cat(" Mean treatment difference: ", result_bmvb$delta$mean_delta, "\n")
cat(" Posterior probability: ", result_bmvb$delta$pop, "\n")

cat("  Logistic regression:\n   ")
cat(" Mean treatment difference: ", result_bglm$delta$mean_delta, "\n")
cat(" Posterior probability: ", result_bglm$delta$pop, "\n")

cat("  Difference:\n            ")
cat(" Treatment difference:", abs(result_bglm$delta$mean_delta - result_bmvb$delta$mean_delta), "\n")
cat(" Posterior probability:", abs(result_bglm$delta$pop - result_bmvb$delta$pop), "\n")

## ----predictions, results='hide', message=FALSE, eval=TRUE--------------------
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

## ----use_predictions, eval=TRUE-----------------------------------------------
# Treatment effect samples for age = 65
delta_samples <- result_samples$samples$delta

# Custom probability: effect on pain relief > 0.2
cat("P(effect on pain > 0.2 | age=65):", mean(delta_samples[, 1] > 0.2), "\n")

# Credible interval
ci <- apply(delta_samples, 2, quantile, probs = c(0.025, 0.975))
cat("\n95% Credible Intervals for age=65:\n")
print(ci)

