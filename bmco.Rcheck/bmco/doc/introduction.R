## ----echo=FALSE---------------------------------------------------------------
biblio <- RefManageR::ReadBib("../inst/REFERENCES.bib", check = "error")

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
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
summary(result)


## ----bmvb_example, eval=TRUE--------------------------------------------------
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

## ----bglm_example, eval=TRUE--------------------------------------------------
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
summary(result_bglm) 

## ----bglmm_example, eval=TRUE-------------------------------------------------
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
summary(result_bglmm)

## ----rule_all-----------------------------------------------------------------
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

## ----rule_any-----------------------------------------------------------------
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

## ----rule_comp----------------------------------------------------------------
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
cat("Weighted effect:", result_comp$delta$w_delta, "\n")

## ----test_right---------------------------------------------------------------
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",      
  grp_b = "drug", # Put the group you expect to be better here
  y_vars = c("symptom_relief", "quality_of_life"),
  test = "right_sided",  # P(drug > placebo)
  n_it = 5000
)

## ----test_left----------------------------------------------------------------
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  test = "left_sided",  # P(placebo > drug)
  n_it = 5000
)

## ----output_structure---------------------------------------------------------
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

# Posterior estimates for each group
result$estimates

# Sample sizes
result$sample_sizes

# Treatment effect
result$delta

# Test information
result$info

## ----samples------------------------------------------------------------------
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

# Custom probability calculations
samples <- result_samples$samples$delta

# P(effect on symptom relief > 0.1)
cat("P(delta[symptom] > 0.1):", mean(samples[, 1] > 0.1), "\n")

# P(effect on QoL > 0.15)
cat("P(delta[QoL] > 0.15):", mean(samples[, 2] > 0.15), "\n")

# Credible intervals
cat("\n95% Credible Intervals:\n")
print(apply(samples, 2, quantile, probs = c(0.025, 0.975)))

## ----mcmc_settings, eval=TRUE-------------------------------------------------
bmvb(
  data = trial_data,
  grp = "treatment",
  grp_a = "placebo",
  grp_b = "drug",
  y_vars = c("symptom_relief", "quality_of_life"),
  n_it = 10000  # More iterations for smoother estimates
)

## ----mcmc_settings_reg, eval = FALSE------------------------------------------
# bglm(
#   # ... arguments ...
#   n_burn = 2000,  # Burn-in / warm-up iterations
#   n_it = 5000,    # Sampling iterations
#   n_thin = 1      # Keep every nth sample (reduces memory)
# )

## ----missing_data, eval=TRUE--------------------------------------------------
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

# Sample sizes are reduced
result_missing$sample_sizes

## ----frequentist, eval=TRUE---------------------------------------------------
# Separate tests for each outcome
chisq_symptom <- chisq.test(table(trial_data$treatment, trial_data$symptom_relief))
chisq_qol <- chisq.test(table(trial_data$treatment, trial_data$quality_of_life))

cat("Frequentist p-values:\n")
cat("  Symptom relief:", chisq_symptom$p.value, "\n")
cat("  Quality of life:", chisq_qol$p.value, "\n")
cat("  Bonferroni-adjusted alpha: 0.025 (needed for Any-rule only; not needed for All-rule)\n\n")

## ----bayesian, eval=TRUE------------------------------------------------------
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
cat("  P(drug better on BOTH):", result_bayes$delta$pop, "\n")

