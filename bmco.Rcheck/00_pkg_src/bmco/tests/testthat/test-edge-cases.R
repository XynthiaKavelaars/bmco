#### Comprehensive Edge Case Tests for bmco ####

library(testthat)
library(bmco)

# ---------------------------------------------------------------------------
# Helper: checks whether all warnings originating from expr are given.
# use this when a call generates multiple warnings and you wish to verify presence of all of them.
# ---------------------------------------------------------------------------
expect_all_warnings <- function(expr, patterns) {
  warns <- testthat::capture_warnings(expr)
  expect_true(
    all(colSums(sapply(patterns, grepl, warns)) > 0),
    info = paste("Not all expected warnings found.\nFound:", paste(warns, collapse = "\n"))
  )
}

# =============================================================================
# PART 1: DATA STRUCTURE EDGE CASES
# =============================================================================

test_that("error when data is NULL", {
  expect_error(
    bmvb(data = NULL, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "data frame"
  )
})

test_that("error when data is not a data frame", {
  expect_error(
    bmvb(data = list(a = 1, b = 2), grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "data frame"
  )
})

test_that("error when data has zero rows", {
  empty_data <- data.frame(grp = character(0), y1 = integer(0), y2 = integer(0))
  expect_error(
    bmvb(data = empty_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "Group.*not found in grouping variable"
  )
})

# =============================================================================
# PART 2: GROUP VARIABLE EDGE CASES
# =============================================================================

test_that("error when group variable doesn't exist", {
  test_data <- data.frame(
    treatment = rep(c("A", "B"), each = 20),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "nonexistent", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "not found in data"
  )
})

test_that("error when grp_a value not found", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "X", grp_b = "B",
         y_vars = c("y1", "y2")),
    "not found in grouping variable"
  )
})

test_that("error when grp_b value not found", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "Y",
         y_vars = c("y1", "y2")),
    "not found in grouping variable"
  )
})

test_that("error when grp has only one unique value", {
  test_data <- data.frame(
    grp = rep("A", 40),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "not found in grouping variable"
  )
})

test_that("error when one group is empty after NA removal", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = c(rep(NA, 20), rbinom(20, 1, 0.5)),  # All A are NA
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    suppressWarnings(bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"))),
    "No observations in group"
  )
})

# =============================================================================
# PART 3: OUTCOME VARIABLE EDGE CASES
# =============================================================================

test_that("error when outcome variable doesn't exist", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y_missing")),
    "not found in data"
  )
})

test_that("error for more than 2 outcomes", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rbinom(40, 1, 0.5),
    y2 = rbinom(40, 1, 0.5),
    y3 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2", "y3")),
    "only supported for 2 outcomes"
  )
})

test_that("error for non-binary outcome values", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = sample(0:2, 40, replace = TRUE),  # 0, 1, 2
    y2 = rbinom(40, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2")),
    "must be binary"
  )
})

test_that("warning for outcome with no variation (all zeros) — three warnings expected", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rep(0, 40),
    y2 = rbinom(40, 1, 0.5)
  )
  # check_input geeft: globale warning (sectie 2) + per-groep warning (sectie 3)
  expect_all_warnings(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), n_it = 100),
    c(
      "no variation \\(all 0\\)",
      "no variation in group A \\(all 0\\)",
      "no variation in group B \\(all 0\\)"
    )
  )
})

test_that("warning for outcome with no variation (all ones) — drie warnings verwacht", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 20),
    y1 = rep(1, 40),
    y2 = rbinom(40, 1, 0.5)
  )
  expect_all_warnings(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), n_it = 100),
    c(
      "no variation \\(all 1\\)",
      "no variation in group A \\(all 1\\)",
      "no variation in group B \\(all 1\\)"
    )
  )
})

# =============================================================================
# PART 4: COVARIATE EDGE CASES (bglm)
# =============================================================================

test_that("error when covariate doesn't exist", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "nonexistent", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 0, n_burn = 10, n_it = 20),
    "not found in data"
  )
})

test_that("warning when covariate has no variation", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rep(50, 60),  # No variation
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_all_warnings(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 50, n_burn = 10, n_it = 20),
    c("no variation", "sample size", "converged")
  )
})

test_that("warning when covariate is highly collinear with group", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    # x nearly perfectly predicts group
    x = rep(c(1, 2), each = 30) + rnorm(60, 0, 0.01),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_all_warnings(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = c(-Inf, Inf), n_burn = 10, n_it = 20),
    c("High correlation", "sample size", "converged")
  )
})

test_that("error for discrete covariate with more than 2 levels", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    severity = factor(sample(c("mild", "moderate", "severe"), 60, replace = TRUE)),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "severity", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 0, n_burn = 10, n_it = 20),
    "only binary covariates.*2 categories.*supported"
  )
})

test_that("error for discrete covariate with only 1 level", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    constant = factor(rep("same", 60)),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "constant", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 0, n_burn = 10, n_it = 20),
    "only 1 category.*At least 2 categories"
  )
})

test_that("error when x_method = 'Value' but x_def is a vector", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60, 50, 10),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = c(40, 60), n_burn = 10, n_it = 20),
    "x_method"
  )
})

test_that("error when x_method = 'Empirical' but x_def is a scalar", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60, 50, 10),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = 50, n_burn = 10, n_it = 20),
    "x_method"
  )
})

test_that("error when x_def range is inverted", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60, 50, 10),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = c(60, 40), n_burn = 10, n_it = 20),
    "increasing order"
  )
})

test_that("error when x_def is below minimum", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60, 50, 10),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = c(-100, min(test_data$age) - 1),
         n_burn = 10, n_it = 20),
    "below minimum observed value"
  )
})

# =============================================================================
# PART 5: CLUSTER EDGE CASES (bglmm)
# =============================================================================

test_that("error when id_var is NULL for bglmm", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglmm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
          id_var = NULL, x_var = "age", y_vars = c("y1", "y2"),
          x_method = "Value", x_def = 0, n_burn = 10, n_it = 20),
    "id_var must be specified"
  )
})

test_that("error when id_var doesn't exist", {
  test_data <- data.frame(
    cluster = factor(rep(1:6, each = 10)),
    grp = rep(c("A", "B"), times = 30),
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglmm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
          id_var = "nonexistent", x_var = "age", y_vars = c("y1", "y2"),
          x_method = "Value", x_def = 0, n_burn = 10, n_it = 20),
    "not found in data"
  )
})

test_that("warning for very few clusters (J < 5)", {
  test_data <- data.frame(
    cluster = factor(rep(1:4, each = 15)),
    grp = rep(c("A", "B"), times = 30),
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )

    warns <- testthat::capture_warnings(
    bglmm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
          id_var = "cluster", x_var = "age", y_vars = c("y1", "y2"),
          x_method = "Value", x_def = 0, n_burn = 20, n_it = 1000,
          return_diagnostics = FALSE)
  )
  expect_true(any(grepl("Very few clusters", warns)))
})

test_that("warning for small cluster size", {
  test_data <- data.frame(
    cluster = factor(c(rep(1, 2), rep(2:10, each = 10))),
    grp = rep(c("A", "B"), length.out = 92),
    age = rnorm(92),
    y1 = rbinom(92, 1, 0.5),
    y2 = rbinom(92, 1, 0.5)
  )
  warns <- testthat::capture_warnings(
    bglmm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
          id_var = "cluster", x_var = "age", y_vars = c("y1", "y2"),
          x_method = "Value", x_def = 0, n_burn = 100, n_it = 1000,
          return_diagnostics = FALSE)
  )
  expect_true(any(grepl("very few observations", warns)))
})

test_that("warning for cluster with only one group", {
  test_data <- data.frame(
    cluster = factor(rep(1:6, each = 10)),
    grp = c(rep("A", 10), rep(c("A", "B"), times = 25)),  # Cluster 1 all A
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  warns <- testthat::capture_warnings(
    bglmm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
          id_var = "cluster", x_var = "age", y_vars = c("y1", "y2"),
          x_method = "Value", x_def = 0, n_burn = 100, n_it = 1000,
          return_diagnostics = FALSE)
  )
  expect_true(any(grepl("only one group", warns)))
})

# =============================================================================
# PART 6: DECISION RULE EDGE CASES
# =============================================================================

test_that("error when rule = 'Comp' but w is NULL", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_message(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), rule = "Comp", w = NULL, n_it = 100),
    "equal weights"
  )
})

test_that("error when w has wrong length", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), rule = "Comp", w = c(0.5), n_it = 100),
    "not equal to number of outcomes"
  )
})

test_that("error when w contains negative values", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), rule = "Comp", w = c(-0.5, 0.5), n_it = 100),
    "non-negative"
  )
})

# =============================================================================
# PART 7: MCMC SETTINGS EDGE CASES
# =============================================================================

test_that("error when n_it = 0", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), n_it = 0),
    "iterations"
  )
})

test_that("warning when n_burn = 0", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60),
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
     bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 0, n_burn = 0, n_it = 100, n_thin = 1)
    ,
    "burnin"
  )
})

# NOT RUN (computation time)
#test_that("message when n_it is very large without thinning", {
#  test_data <- data.frame(
#    grp = rep(c("A", "B"), each = 30),
#    age = rnorm(60),
#    y1 = rbinom(60, 1, 0.5),
#    y2 = rbinom(60, 1, 0.5)
#  )
#  expect_message(
#    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
#         x_var = "age", y_vars = c("y1", "y2"),
#         x_method = "Value", x_def = 0, n_burn = 10, n_it = 1+1e6, n_thin = 1),
#    "thinning"
#  )
#})

# =============================================================================
# PART 8: PERFECT SEPARATION CASES
# =============================================================================

test_that("model runs with perfect group separation in outcomes", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    y1 = c(rep(0, 30), rep(1, 30)),  # Perfect separation
    y2 = rbinom(60, 1, 0.5)
  )
  # Should run but may have poor diagnostics
  expect_s3_class(
    suppressWarnings(bmvb(data = test_data, grp = "grp", grp_a = "A",
                          grp_b = "B", y_vars = c("y1", "y2"), n_it = 100)),
    "bmvb"
  )
})


# =============================================================================
# PART 9: NUMERICAL STABILITY TESTS
# =============================================================================

test_that("error when covariate has extreme values", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = c(rnorm(30, 50000, 100), rnorm(30, 50000, 100)),  # ~50k (extreme)
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = 50000,
         n_burn = 10, n_it = 20),
    "extreme values.*numerical instability"
  )
})

test_that("error when covariate range is very wide", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = c(runif(30, 0, 20000), runif(30, 0, 20000)),  # Range = 20k (too wide)
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "age", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = c(-Inf, Inf),
         n_burn = 10, n_it = 20),
    "very wide range.*numerical instability"
  )
})

test_that("model works with moderate values and range", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age = rnorm(60, 50, 10),  # Moderate: mean=50, range ~20-80
    y1 = rbinom(60, 1, 0.5),
    y2 = rbinom(60, 1, 0.5)
  )
  # Should run without error
  expect_s3_class(
    suppressWarnings(bglm(data = test_data, grp = "grp", grp_a = "A",
                          grp_b = "B", x_var = "age", y_vars = c("y1", "y2"),
                          x_method = "Value", x_def = 50,
                          n_burn = 50, n_it = 100)),
    "bglm"
  )
})

test_that("standardized covariates work well", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 30),
    age_raw = rnorm(60, 50, 10)
  )
  # Standardize (mean=0, sd=1)
  test_data$age <- scale(test_data$age_raw)[, 1]
  test_data$y1 <- rbinom(60, 1, 0.5)
  test_data$y2 <- rbinom(60, 1, 0.5)

  expect_s3_class(
    suppressWarnings(bglm(data = test_data, grp = "grp", grp_a = "A",
                          grp_b = "B", x_var = "age", y_vars = c("y1", "y2"),
                          x_method = "Value", x_def = 0,
                          n_burn = 50, n_it = 100)),
    "bglm"
  )
})
# =============================================================================
# PART 10: EXTREME IMBALANCE
# =============================================================================

test_that("model handles extreme group imbalance", {
  test_data <- data.frame(
    grp = c(rep("A", 10), rep("B", 100)),
    y1 = rbinom(110, 1, 0.5),
    y2 = rbinom(110, 1, 0.5)
  )
  expect_s3_class(
    suppressWarnings(bmvb(data = test_data, grp = "grp", grp_a = "A",
                          grp_b = "B", y_vars = c("y1", "y2"), n_it = 100)),
    "bmvb"
  )
})
