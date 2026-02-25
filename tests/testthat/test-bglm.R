#### Unit Tests for bglm() ####
# Covers all checks in check_input() sections 1-6 as well as the
# design-matrix creation (section 8) that now lives inside check_input().

library(testthat)
library(bmco)

# ---------------------------------------------------------------------------
# Helper: controleer dat ALLE patronen voorkomen in de warnings van expr.
# ---------------------------------------------------------------------------
expect_all_warnings <- function(expr, patterns) {
  warns <- testthat::capture_warnings(expr)
  expect_true(
    all(colSums(sapply(patterns, grepl, warns)) > 0),
    info = paste("Niet alle verwachte warnings gevonden.\nGevonden:", paste(warns, collapse = "\n"))
  )
}

# ---------------------------------------------------------------------------
# Shared helper: minimal valid dataset with a continuous covariate
# ---------------------------------------------------------------------------
make_bglm_data <- function(n = 80, seed = 123) {
  set.seed(seed)
  data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    x   = rnorm(n, 0, 1),
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
}

#### Test 1: Basic Functionality ####
test_that("bglm runs with valid input", {
  test_data <- make_bglm_data()
  result <- bglm(
    data   = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var  = "x", y_vars = c("y1", "y2"), n_it = 1000
  )
  expect_s3_class(result, "bglm")
  expect_type(result, "list")
  expect_true(all(c("estimates", "sample_sizes", "delta", "info") %in% names(result)))
})

#### Test 2: Return Samples ####
test_that("bglm returns samples when requested", {
  test_data <- make_bglm_data(seed = 124)
  result <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), n_it = 1000, return_samples = TRUE
  )
  expect_true("samples" %in% names(result))
  expect_type(result$samples$theta_a, "double")
  expect_type(result$samples$theta_b, "double")
  expect_type(result$samples$delta,   "double")
  expect_equal(nrow(result$samples$theta_a), 2000)
})

#### Test 3: Posterior Probability Range ####
test_that("posterior probability is between 0 and 1", {
  test_data <- make_bglm_data(seed = 125)
  result <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), n_it = 1000
  )
  expect_gte(result$delta$pop, 0)
  expect_lte(result$delta$pop, 1)
  expect_type(result$delta$pop, "double")
})

#### Test 4: Decision Rules ####
test_that("all decision rules work", {
  test_data <- make_bglm_data(seed = 126)

  result_all <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), rule = "All", n_it = 1000
  )
  expect_equal(result_all$info$rule, "All")

  result_any <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), rule = "Any", n_it = 1000
  )
  expect_equal(result_any$info$rule, "Any")

  result_comp <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"),
    rule = "Comp", w = c(0.6, 0.4), n_it = 1000
  )
  expect_equal(result_comp$info$rule, "Comp")
  expect_equal(result_comp$info$w, c(0.6, 0.4))
  expect_true("w_delta" %in% names(result_comp$delta))
})

#### Test 5: Test Directions ####
test_that("test directions work correctly", {
  test_data <- make_bglm_data(seed = 127)

  result_right <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), test = "right_sided", n_it = 1000
  )
  expect_equal(result_right$info$test, "right_sided")
  expect_match(result_right$info$test_label, "B > A")

  result_left <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), test = "left_sided", n_it = 1000
  )
  expect_equal(result_left$info$test, "left_sided")
  expect_match(result_left$info$test_label, "A > B")
})

#### Test 6: Equal Weights Default ####
test_that("equal weights are used by default for Comp rule", {
  test_data <- make_bglm_data(seed = 128)
  result <- suppressMessages(
    bglm(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      x_var = "x", y_vars = c("y1", "y2"),
      rule = "Comp", w = NULL, n_it = 1000
    )
  )
  expect_equal(result$info$w, c(0.5, 0.5))
})

#### Test 7: Prior Parameters ####
test_that("different prior parameters can be specified", {
  test_data <- make_bglm_data(seed = 129)
  result <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"),
    b_mu0 = 1, b_sigma0 = 1e-1, n_it = 1000
  )
  expect_s3_class(result, "bglm")
})

#### Test 8: Sample Sizes are Reported ####
test_that("sample sizes are correctly reported", {
  set.seed(130)
  n_a <- 30; n_b <- 40
  test_data <- data.frame(
    grp = c(rep("A", n_a), rep("B", n_b)),
    x   = rnorm(n_a + n_b, 0, 1),
    y1  = rbinom(n_a + n_b, 1, 0.4),
    y2  = rbinom(n_a + n_b, 1, 0.6)
  )
  result <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), n_it = 1000
  )
  expect_equal(result$sample_sizes$n_a, n_a)
  expect_equal(result$sample_sizes$n_b, n_b)
})

#### Test 9a: Missing Data on y ####
test_that("missing data on y triggers warning and reduces sample size", {
  test_data <- make_bglm_data(n = 100, seed = 131)
  test_data$y1[c(1, 3, 5)] <- NA
  test_data$y2[c(2, 4)]    <- NA

  expect_warning(
    result <- bglm(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      x_var = "x", y_vars = c("y1", "y2"), n_it = 1000
    ),
    "Missing data detected"
  )
  expect_true(result$sample_sizes$n_a < 50)
  expect_s3_class(result, "bglm")
})

#### Test 9b: Missing Data on x ####
test_that("missing data on x triggers warning and reduces total sample size", {
  test_data <- make_bglm_data(n = 100, seed = 135)
  test_data$x[c(5, 10, 15)] <- NA   # Three missing covariate values

  expect_warning(
    result <- bglm(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      x_var = "x", y_vars = c("y1", "y2"), n_it = 1000
    ),
    "Missing data detected"
  )
  # Total observations used must be < 100
  expect_true(result$sample_sizes$n_a + result$sample_sizes$n_b < 100)
  expect_s3_class(result, "bglm")
})

#### Test 9c: Design Matrix Returned by check_input ####
test_that("check_input returns x_data and x_names for bglm", {
  test_data <- make_bglm_data(n = 60, seed = 136)
  validated <- bmco:::check_input(
    data     = test_data,
    grp      = "grp", grp_a = "A", grp_b = "B",
    y_vars   = c("y1", "y2"),
    test     = "right_sided", rule = "All", w = NULL,
    analysis = "bglm",
    x_var    = "x", x_method = "Empirical", x_def = c(-Inf, Inf),
    n_it     = 100
  )
  # x_data is a matrix with 4 columns (Intercept, grp, x, grp:x)
  expect_true(is.matrix(validated$x_data))
  expect_equal(ncol(validated$x_data), 4)
  expect_equal(length(validated$x_names), 4)
  # Row count must equal n_a + n_b
  expect_equal(nrow(validated$x_data), validated$n_a + validated$n_b)
})

#### Test 9d: x_data and y_a / y_b are consistent after joint NA removal ####
test_that("x_data rows match n_a + n_b after joint x/y NA removal", {
  test_data <- make_bglm_data(n = 100, seed = 137)
  test_data$x[1:5]  <- NA   # missing x in group A
  test_data$y1[51]  <- NA   # missing y in group B
  validated <- suppressWarnings(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), n_it = 1000
    )
  )
  expect_equal(nrow(validated$x_data), validated$n_a + validated$n_b)
  expect_equal(nrow(validated$y_a), validated$n_a)
  expect_equal(nrow(validated$y_b), validated$n_b)
})

#### Test 10: Error Handling – General Input Checks ####

test_that("error for missing group variable", {
  test_data <- data.frame(x = rnorm(2), y1 = c(0, 1), y2 = c(1, 0))
  expect_error(
    bglm(data = test_data, grp = "missing_var", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2")),
    "not found in data"
  )
})

test_that("error for missing outcome variable", {
  test_data <- data.frame(grp = c("A", "B"), x = rnorm(2), y1 = c(0, 1))
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y_missing")),
    "not found in data"
  )
})

test_that("error for invalid group value", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "C", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2")),
    "not found in grouping variable"
  )
})

test_that("error for more than 2 outcomes", {
  set.seed(1)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25), x = rnorm(50),
    y1 = rbinom(50, 1, 0.5), y2 = rbinom(50, 1, 0.5), y3 = rbinom(50, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2", "y3")),
    "only supported for 2 outcomes"
  )
})

test_that("error for wrong weight vector length", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), rule = "Comp", w = c(0.5)),
    "not equal to number of outcomes"
  )
})

test_that("error for negative weights", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), rule = "Comp", w = c(-0.5, 0.5)),
    "non-negative"
  )
})

test_that("error for non-binary outcomes", {
  set.seed(1)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25), x = rnorm(50),
    y1 = sample(0:2, 50, replace = TRUE), y2 = rbinom(50, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2")),
    "must be binary"
  )
})

#### Test 11: Error Handling – bglm-Specific Checks (section 6 of check_input) ####

test_that("error when x_var is not found in data", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "nonexistent", y_vars = c("y1", "y2")),
    "not found in data"
  )
})

test_that("error for unknown measurement level (list column)", {
  test_data <- make_bglm_data()
  test_data$x <- as.list(test_data$x)   # neither numeric nor factor
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2")),
    "unknown measurement level"
  )
})

test_that("error when x_def length 1 but x_method is not 'Value'", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = 0),
    "x_method"
  )
})

test_that("error when x_def length 2 but x_method is 'Value'", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"),
         x_method = "Value", x_def = c(-1, 1)),
    "x_method"
  )
})

test_that("error when x_def range is not in increasing order", {
  test_data <- make_bglm_data()
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"),
         x_method = "Empirical", x_def = c(1, -1)),
    "increasing order"
  )
})

test_that("warning for covariate with no variation", {
  set.seed(140)
  n <- 60
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    x   = rep(5, n),           # Constant covariate
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), n_it = 1000
    ),
    "no variation"
  )
})

test_that("warning for high collinearity between group and covariate", {
  set.seed(141)
  n <- 100
  grp_ind <- rep(c(0, 1), each = n / 2)
  test_data <- data.frame(
    grp = ifelse(grp_ind == 1, "B", "A"),
    x   = grp_ind + rnorm(n, 0, 0.01),  # r approximately 0.999
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), n_it = 1000
    ),
    "High correlation"
  )
})

test_that("warning when x_def lower bound is below observed minimum", {
  test_data <- make_bglm_data(n = 100, seed = 142)
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-100, Inf), n_it = 1000
    ),
    "below minimum observed value"
  )
})

test_that("warning when x_def upper bound is above observed maximum", {
  test_data <- make_bglm_data(n = 100, seed = 143)
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, 100), n_it = 1000
    ),
    "above maximum observed value"
  )
})

test_that("warning for few unique values in continuous covariate with Empirical method", {
  set.seed(144)
  n <- 60
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    x   = sample(1:3, n, replace = TRUE),   # Only 3 unique values
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), n_it = 1000
    ),
    "Few unique values"
  )
})

test_that("message for large n_it without thinning", {
  test_data <- make_bglm_data()
  expect_message(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), n_it = 2e6, n_thin = 1
    ),
    "thinning"
  )
})

test_that("error for discrete covariate with more than 2 categories", {
  set.seed(150)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 40),
    x   = factor(sample(c("low", "medium", "high"), 80, replace = TRUE)),  # 3 categories
    y1  = rbinom(80, 1, 0.5),
    y2  = rbinom(80, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), x_method = "Value", x_def = 0,
         n_burn = 100, n_it = 100),
    "only binary covariates.*2 categories.*supported"
  )
})

test_that("error for discrete covariate with only 1 category", {
  set.seed(151)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 40),
    x   = factor(rep("only_one_level", 80)),  # 1 category
    y1  = rbinom(80, 1, 0.5),
    y2  = rbinom(80, 1, 0.5)
  )
  expect_error(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), x_method = "Value", x_def = 0,
         n_burn = 100, n_it = 100),
    "only 1 category.*At least 2 categories"
  )
})
#### Test 12: Print Method ####
test_that("print method works", {
  test_data <- make_bglm_data(n = 60, seed = 132)
  result <- bglm(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    x_var = "x", y_vars = c("y1", "y2"), n_it = 1000
  )
  expect_output(print(result), "Bayesian Multivariate Logistic Regression")
  expect_output(print(result), "Group Estimates")
  expect_output(print(result), "Posterior probability")
})

#### Test 13: Warning for Small Sample Size ####
# check_input geeft één warning per groep: "Small sample size in group A (n = 8)"
# EN "Small sample size in group B (n = 8)". De oude constructie
# suppressWarnings(..., classes = ...) BINNENIN expect_warning onderdrukt
# mogelijk alle warnings vóór expect_warning ze kan zien (R < 4.2). Fix:
# gebruik capture_warnings() om beide groep-warnings te verifiëren.
test_that("warning for small sample size in both groups", {
  set.seed(133)
  test_data <- data.frame(
    grp = c(rep("A", 8), rep("B", 8)),
    x   = rnorm(16),
    y1  = rbinom(16, 1, 0.5),
    y2  = rbinom(16, 1, 0.5)
  )
  expect_all_warnings(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), n_it = 1000),
    c("Small sample size.*A", "Small sample size.*B")
  )
})

#### Test 14: Warning for No Variation in Outcome ####
# check_input geeft drie warnings wanneer y1 constant is:
#  1. Globaal (sectie 2): "has no variation (all 1)"
#  2. Binnen groep A (sectie 3): "has no variation in group A (all 1)"
#  3. Binnen groep B (sectie 3): "has no variation in group B (all 1)"
test_that("warning for outcome with no variation — drie warnings verwacht", {
  set.seed(134)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    x   = rnorm(50),
    y1  = rep(1, 50),
    y2  = rbinom(50, 1, 0.5)
  )
  expect_all_warnings(
    bglm(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2"), n_it = 1000),
    c(
      "no variation \\(all 1\\)",
      "no variation in group A \\(all 1\\)",
      "no variation in group B \\(all 1\\)"
    )
  )
})

#### Test 15: Discrete (Factor) Covariate ####
test_that("discrete (factor) covariate is accepted and measurement level is 'discrete'", {
  set.seed(145)
  n <- 80
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n / 2),
    x   = factor(sample(c("low", "high"), n, replace = TRUE)),
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglm", x_var = "x", x_method = "Value", x_def = 0,
    n_it = 1000
  )
  expect_equal(validated$ml, "discrete")
  expect_true(is.matrix(validated$x_data))
})
