## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----echo=FALSE---------------------------------------------------------------
biblio <- RefManageR::ReadBib("../inst/REFERENCES.bib", check = "error")

## ----setup--------------------------------------------------------------------
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


## ----fit_bglmm, results='hide', message=FALSE, warning=FALSE, eval=TRUE-------
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

## ----show_results, eval=TRUE--------------------------------------------------
print(result)
summary(result)

## ----value_method, results='hide', message=FALSE, warning=FALSE, eval = FALSE----
# set.seed(2027)
# result_value <- bglmm(
#   # ... other arguments ...
#   x_method = "Value",
#   x_def = 50,  # Average ability
# )

## ----show_value, eval = FALSE-------------------------------------------------
# cat("Effect for students with baseline ability = 50:\n")
# cat("  Treatment difference:", result_value$delta$mean_delta, "\n")
# cat("  Posterior probability:", result_value$delta$pop, "\n")

## ----empirical_method, results='hide', message=FALSE, warning=FALSE, eval = FALSE----
# set.seed(2028)
# result_empirical <- bglmm(
#   # ... other arguments ...
#   x_method = "Empirical",
#   x_def = c(40,60) # Ability range 40-60
# )

## ----analytical_method, results='hide', message=FALSE, warning=FALSE, eval = FALSE----
# set.seed(2029)
# result_analytical <- bglmm(
#   # ... other arguments ...
#     x_method = "Analytical",
#     x_def = c(40,60)
#  )

## ----rule_all, eval = FALSE---------------------------------------------------
# result_all <- bglmm(
#   data = study_data,
#   # ... other arguments ...
#   rule = "All"  # P(new_method better on BOTH outcomes)
# )

## ----rule_any, eval = FALSE---------------------------------------------------
# result_any <- bglmm(
#   data = study_data,
#   # ... other arguments ...
#   rule = "Any"  # P(new_method better on AT LEAST ONE outcome)
# )

## ----rule_comp, eval = FALSE--------------------------------------------------
# result_comp <- bglmm(
#   data = study_data,
#   # ... other arguments ...
#   rule = "Comp",
#   w = c(0.6, 0.4)  # Weights for math and reading
# )

## ----default_priors, eval=TRUE------------------------------------------------
# Default
p_fixed <- 2
b_mu0 <- rep(0, p_fixed)                 # Zero prior mean
b_sigma0 <- solve(diag(1e1, p_fixed))    # Small prior precision (large prior variance on the log-odds scale). 
                                         # Weakly informative. 
                                         # solve() generates precision matrix needed as input

## ----custom_prior_1, results='hide', message=FALSE, warning=FALSE, eval=TRUE----
# Assume treatment improves outcomes (positive coefficient on group)

custom_b_mu0 <- c(0, 0.5)                 # Expect positive group effect 
custom_b_sigma0 <- solve(diag(c(1e1, 5))) # More certainty on group effect

## ----random_default, eval=TRUE------------------------------------------------
p_random <- 2
g_mu0 <- rep(0, p_random)               # Prior mean
g_sigma0 <- solve(diag(1e1, p_random))  # Weakly informative prior variance; 
                                        # solve() transforms variance matrix to precision matrix needed as input

## ----iwishart_default, eval=TRUE----------------------------------------------
nu0 <- p_random                # Degrees of freedom (dimension)
tau0 <- diag(1e-1, p_random)   # Scale matrix

## ----prior_sensitivity, results='hide', message=FALSE, warning=FALSE, eval = FALSE----
# # -----------------------------------------------------------------
# # Three priors that pull in clearly different directions.
# # -----------------------------------------------------------------
# 
# # 1. Weakly informative (neutral starting point)
# #    Diffuse on fixed effects, diffuse IW for random effects.
# result_weak <- bglmm(
#   # ... other arguments ...
#   b_mu0    = rep(0, p_fixed),
#   b_sigma0 = solve(diag(10, p_fixed)),
#   g_mu0    = rep(0, p_random),
#   g_sigma0 = solve(diag(10, p_random)),
#   nu0  = p_random,
#   tau0 = diag(0.1, p_random),
# )
# 
# # 2. Skeptical prior
# #    Tight variance so the prior resists positive or negative evidence.
# #    Also a small tau0, implying homogeneity between schools.
# result_skeptical <- bglmm(
#   # ... other arguments ...
#   b_mu0    = rep(0, p_fixed),
#   b_sigma0 = solve(diag(c(1, 1), p_fixed)),  # Tight — hard(er) to overcome
#   g_mu0    = rep(0, p_random),
#   g_sigma0 = solve(diag(c(1, 1), p_random)), # Tight — hard(er) to overcome
#   nu0  = p_random + 4,                 # More informative IW
#   tau0 = diag(0.05, p_random),         # Expects small school variation
# )
# 

## ----prior_check, eval = FALSE------------------------------------------------
#    # Run with default priors
#    # Run with informative priors
#    # Compare posterior probabilities and effect sizes
#    # If similar → robust; if different → prior-dependent

## ----diagnose, results='hide', message=FALSE, warning=FALSE, eval=TRUE--------
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

## ----show_diag, eval=TRUE-----------------------------------------------------
# Check convergence
cat("Random effects convergence:\n")
print(result_diag$diags$g$convergence)

cat("\nRandom effects covariance convergence:\n")
print(result_diag$diags$tau$convergence)


## ----plot_convergence, eval = FALSE-------------------------------------------
# plot(result_diag)                                    # All
# plot(result_diag, type = "trace")                    # Trace only
# plot(result_diag, type = "density")                  # Density only
# 
# plot(result_diag, which = "fixed")                   # Fixed effects
# plot(result_diag, which = "random")                  # Random effects
# plot(result_diag, which = "variance")                # Variance components
# plot(result_diag, which = "all")                     # All
# 
# # Combine type + which
# plot(result_diag, type = "trace", which = "all")          # Trace for all parameters
# plot(result_diag, type = "density", which = "variance")   # Variance densities

## ----samples, results='hide', message=FALSE, warning=FALSE, eval = FALSE------
# set.seed(2031)
# result_samples <- bglmm(
# # ... other arguments ...
#   return_samples = TRUE
# )

## ----use_samples, eval = FALSE------------------------------------------------
# # Treatment effect samples
# str(result_samples$samples$delta)
# 
# # Compute custom summaries
# cat("Treatment effect quantiles:\n")
# print(apply(result_samples$samples$delta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
# 
# # Probability that effect on math > 0.1
# cat("\nP(delta_math > 0.1):", mean(result_samples$samples$delta[, 1] > 0.1), "\n")

