#### Unit Tests for bmvb() ####

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

#### Test 1: Basic Functionality ####
test_that("bmvb runs with valid input", {
  # Simulate data
  set.seed(123)
  n <- 100
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 100  # Low for speed
  )

  # Check result structure
  expect_s3_class(result, "bmvb")
  expect_type(result, "list")

  # Check required components
  expect_true("estimates" %in% names(result))
  expect_true("sample_sizes" %in% names(result))
  expect_true("delta" %in% names(result))
  expect_true("info" %in% names(result))
})

#### Test 2: Return Samples ####
test_that("bmvb returns samples when requested", {
  set.seed(124)
  n <- 50
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.4),
    y2 = rbinom(n, 1, 0.6)
  )

  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 100,
    return_samples = TRUE
  )

  expect_true("samples" %in% names(result))
  expect_type(result$samples$theta_a, "double")
  expect_type(result$samples$theta_b, "double")
  expect_type(result$samples$delta, "double")
  expect_equal(nrow(result$samples$theta_a), 100)
})

#### Test 3: Posterior Probability Range ####
test_that("posterior probability is between 0 and 1", {
  set.seed(125)
  n <- 60
  test_data <- data.frame(
    grp = rep(c("X", "Y"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "X",
    grp_b = "Y",
    y_vars = c("y1", "y2"),
    n_it = 100
  )

  expect_true(result$delta$pop >= 0)
  expect_true(result$delta$pop <= 1)
  expect_type(result$delta$pop, "double")
})

#### Test 4: Decision Rules ####
test_that("all decision rules work", {
  set.seed(126)
  n <- 80
  test_data <- data.frame(
    grp = factor(rep(c("A", "B"), each = n/2)),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Test All rule
  result_all <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    rule = "All",
    n_it = 100
  )
  expect_equal(result_all$info$rule, "All")

  # Test Any rule
  result_any <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    rule = "Any",
    n_it = 100
  )
  expect_equal(result_any$info$rule, "Any")

  # Test Comp rule with weights
  result_comp <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    rule = "Comp",
    w = c(0.6, 0.4),
    n_it = 100
  )
  expect_equal(result_comp$info$rule, "Comp")
  expect_equal(result_comp$info$w, c(0.6, 0.4))
  expect_true("w_delta" %in% names(result_comp$delta))
})

#### Test 5: Test Directions ####
test_that("test directions work correctly", {
  set.seed(127)
  n <- 70
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Right-sided test
  result_right <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    test = "right_sided",
    n_it = 100
  )
  expect_equal(result_right$info$test, "right_sided")
  expect_match(result_right$info$test_label, "B > A")

  # Left-sided test
  result_left <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    test = "left_sided",
    n_it = 100
  )
  expect_equal(result_left$info$test, "left_sided")
  expect_match(result_left$info$test_label, "A > B")
})

#### Test 6: Equal Weights Default ####
test_that("equal weights are used by default for Comp rule", {
  set.seed(128)
  n <- 50
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Suppress message about using equal weights
  result <- suppressMessages(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2"),
      rule = "Comp",
      w = NULL,  # Not specified
      n_it = 100
    )
  )

  expect_equal(result$info$w, c(0.5, 0.5))
})

#### Test 7: Prior Parameters ####
test_that("different priors can be specified", {
  set.seed(129)
  n <- 60
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Should work with different priors
  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    prior_a = 1.0,
    prior_b = 1.0,
    n_it = 100
  )

  expect_s3_class(result, "bmvb")
})

#### Test 8: Sample Sizes are Reported ####
test_that("sample sizes are correctly reported", {
  set.seed(130)
  n_a <- 30
  n_b <- 40
  test_data <- data.frame(
    grp = c(rep("A", n_a), rep("B", n_b)),
    y1 = rbinom(n_a + n_b, 1, 0.5),
    y2 = rbinom(n_a + n_b, 1, 0.5)
  )

  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 100
  )

  expect_equal(result$sample_sizes$n_a, n_a)
  expect_equal(result$sample_sizes$n_b, n_b)
})

#### Test 9: Missing Data Handling ####
test_that("missing data is handled correctly", {
  set.seed(131)
  n <- 100
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Add some missing data
  test_data$y1[c(1, 3, 5)] <- NA
  test_data$y2[c(2, 4)] <- NA

  # Should give warning but still work
  expect_warning(
    result <- bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2"),
      n_it = 100
    ),
    "Missing data detected"
  )

  # Sample size should be reduced
  expect_true(result$sample_sizes$n_a < 50)
  expect_s3_class(result, "bmvb")
})

#### Test 10: Error Handling - Invalid Input ####

test_that("error for missing group variable", {
  test_data <- data.frame(y1 = c(0, 1), y2 = c(1, 0))

  expect_error(
    bmvb(
      data = test_data,
      grp = "missing_var",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2")
    ),
    "not found in data"
  )
})

test_that("error for missing outcome variable", {
  test_data <- data.frame(grp = c("A", "B"), y1 = c(0, 1))

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y_missing")
    ),
    "not found in data"
  )
})

test_that("error for invalid group value", {
  test_data <- data.frame(
    grp = c("A", "B"),
    y1 = c(0, 1),
    y2 = c(1, 0)
  )

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "C",  # Doesn't exist
      grp_b = "B",
      y_vars = c("y1", "y2")
    ),
    "not found in grouping variable"
  )
})

test_that("error for more than 2 outcomes", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    y1 = rbinom(50, 1, 0.5),
    y2 = rbinom(50, 1, 0.5),
    y3 = rbinom(50, 1, 0.5)
  )

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2", "y3")
    ),
    "only supported for 2 outcomes"
  )
})

test_that("error for wrong weight vector length", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    y1 = rbinom(50, 1, 0.5),
    y2 = rbinom(50, 1, 0.5)
  )

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2"),
      rule = "Comp",
      w = c(0.5)  # Should be length 2
    ),
    "not equal to number of outcomes"
  )
})

test_that("error for non-binary outcomes", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    y1 = sample(0:2, 50, replace = TRUE),  # Has values 0, 1, 2
    y2 = rbinom(50, 1, 0.5)
  )

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2")
    ),
    "must be binary"
  )
})

test_that("warning for negative prior", {
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    y1 = rbinom(50, 1, 0.5),
    y2 = rbinom(50, 1, 0.5)
  )

  expect_error(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      y_vars = c("y1", "y2"),
      prior_a = -1  # Invalid
    ),
    "must be positive"
  )
})

#### Test 11: Print Method ####
test_that("print method works", {
  set.seed(132)
  n <- 50
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  result <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 100
  )

  # Should not error
  expect_output(print(result), "Multivariate Bernoulli")
  expect_output(print(result), "Group")
  expect_output(print(result), "Posterior probability")
})

#### Test 12: Warnings for Small Sample Size ####
# check_input geeft één warning per groep: "Small sample size in group A" EN
# "Small sample size in group B". Controleer beide met capture_warnings().
test_that("warning for small sample size in both groups", {
  set.seed(133)
  test_data <- data.frame(
    grp = c(rep("A", 8), rep("B", 8)),  # n < 10 per groep
    y1 = rbinom(16, 1, 0.5),
    y2 = rbinom(16, 1, 0.5)
  )
  expect_all_warnings(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), n_it = 100),
    c("Small sample size.*A", "Small sample size.*B")
  )
})

#### Test 13: Warning for No Variation ####
# check_input geeft drie warnings wanneer y1 constant is:
#  1. Globaal (sectie 2): "has no variation (all 1)"
#  2. Binnen groep A (sectie 3): "has no variation in group A (all 1)"
#  3. Binnen groep B (sectie 3): "has no variation in group B (all 1)"
# Controleer alle drie met capture_warnings().
test_that("warning for outcome with no variation — drie warnings verwacht", {
  set.seed(134)
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = 25),
    y1 = rep(1, 50),  # All 1s
    y2 = rbinom(50, 1, 0.5)
  )
  expect_all_warnings(
    bmvb(data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2"), n_it = 100),
    c(
      "no variation \\(all 1\\)",          # globale check (sectie 2)
      "no variation in group A \\(all 1\\)", # per-groep check (sectie 3)
      "no variation in group B \\(all 1\\)"
    )
  )
})
