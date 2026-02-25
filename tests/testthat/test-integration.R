#### Integration Tests for bmco Package ####

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

#### Test 1: Same Data Through All Three Methods ####
test_that("all three methods give consistent results on same data", {
  set.seed(200)
  n_clusters <- 10
  n_per_cluster <- 10

  # Generate hierarchical data
  test_data <- data.frame(
    cluster = rep(1:n_clusters, each = n_per_cluster),
    grp = rep(rep(c("control", "treatment"), each = n_per_cluster/2), n_clusters),
    age = rnorm(n_clusters * n_per_cluster, mean = 50, sd = 10),
    y1 = rbinom(n_clusters * n_per_cluster, 1, 0.5),
    y2 = rbinom(n_clusters * n_per_cluster, 1, 0.5)
  )

  # Test bmvb (ignoring clustering and covariate)
  result_bmvb <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "control",
    grp_b = "treatment",
    y_vars = c("y1", "y2"),
    n_it = 1000
  )

  # Test bglm (ignoring clustering, adjusting for age)
  result_bglm <- suppressWarnings(
    bglm(
      data = test_data,
      grp = "grp",
      grp_a = "control",
      grp_b = "treatment",
      x_var = "age",
      y_vars = c("y1", "y2"),
      x_method = "Empirical",
      x_def = c(-Inf, Inf),
      n_burn = 50,
      n_it = 500
  )
  )


  # Test bglmm (accounting for clustering and age)
  result_bglmm <- suppressWarnings(  # Suppress warnings about few clusters
    bglmm(
      data = test_data,
      grp = "grp",
      grp_a = "control",
      grp_b = "treatment",
      id_var = "cluster",
      x_var = "age",
      y_vars = c("y1", "y2"),
      x_method = "Empirical",
      x_def = c(-Inf, Inf),
      n_burn = 50,
      n_it = 200,
      n_thin = 2,
      return_diagnostics = FALSE
    )
  )

  # All should produce valid results
  expect_s3_class(result_bmvb, "bmvb")
  expect_s3_class(result_bglm, "bglm")
  expect_s3_class(result_bglmm, "bglmm")

  # All should have posterior probabilities
  expect_true(result_bmvb$delta$pop >= 0 && result_bmvb$delta$pop <= 1)
  expect_true(result_bglm$delta$pop >= 0 && result_bglm$delta$pop <= 1)
  expect_true(result_bglmm$delta$pop >= 0 && result_bglmm$delta$pop <= 1)
})

#### Test 2: Samples Compatibility ####
test_that("posterior samples have correct structure across methods", {
  set.seed(201)
  n <- 100
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    age = rnorm(n, 50, 10),
    cluster = rep(1:10, each = n/10),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # bmvb with samples
  result_bmvb <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 1000,
    return_samples = TRUE
  )

  # bglm with samples
  result_bglm <- suppressWarnings(
    bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "age",
    y_vars = c("y1", "y2"),
    x_method = "Value",
    x_def = 50,
    n_burn = 50,
    n_it = 500,
    return_samples = TRUE
  )
  )

  # Both should have samples
  expect_true("samples" %in% names(result_bmvb))
  expect_true("samples" %in% names(result_bglm))

  # Samples should have correct dimensions
  expect_equal(nrow(result_bmvb$samples$delta), 1000)
  expect_equal(ncol(result_bmvb$samples$delta), 2)
  expect_equal(nrow(result_bglm$samples$delta), 1000)
  expect_equal(ncol(result_bglm$samples$delta), 2)
})

#### Test 3: Decision Rules Consistency ####
test_that("decision rules give consistent results across methods", {
  set.seed(202)
  n <- 80
  test_data <- data.frame(
    grp = rep(c("X", "Y"), each = n/2),
    x = rnorm(n),
    y1 = rbinom(n, 1, 0.6),
    y2 = rbinom(n, 1, 0.4)
  )

  # Test All rule in both bmvb and bglm
  bmvb_all <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "X",
    grp_b = "Y",
    y_vars = c("y1", "y2"),
    rule = "All",
    n_it = 1000
  )

  bglm_all <- bglm(
    data = test_data,
    grp = "grp",
    grp_a = "X",
    grp_b = "Y",
    x_var = "x",
    y_vars = c("y1", "y2"),
    x_method = "Empirical",
    x_def = c(-Inf, Inf),
    rule = "All",
    n_burn = 50,
    n_it = 500
  )

  # Both should use All rule
  expect_equal(bmvb_all$info$rule, "All")
  expect_equal(bglm_all$info$rule, "All")

  # Test Comp rule in both
  bmvb_comp <- suppressMessages(
    bmvb(
      data = test_data,
      grp = "grp",
      grp_a = "X",
      grp_b = "Y",
      y_vars = c("y1", "y2"),
      rule = "Comp",
      n_it = 1000
    )
  )

  bglm_comp <- suppressMessages(
    bglm(
      data = test_data,
      grp = "grp",
      grp_a = "X",
      grp_b = "Y",
      x_var = "x",
      y_vars = c("y1", "y2"),
      x_method = "Value",
      x_def = 0,
      rule = "Comp",
      n_burn = 50,
      n_it = 500
    )
  )

  # Both should use Comp rule and have weights
  expect_equal(bmvb_comp$info$rule, "Comp")
  expect_equal(bglm_comp$info$rule, "Comp")
  expect_true("w_delta" %in% names(bmvb_comp$delta))
  expect_true("w_delta" %in% names(bglm_comp$delta))
})

#### Test 4: Test Label Consistency ####
test_that("test labels are correctly formed across all methods", {
  set.seed(203)
  n <- 60
  test_data <- data.frame(
    treatment = factor(c(rep("placebo", n/2), rep("drug", n/2))),
    cluster = rep(1:6, each = n/6),
    age = rnorm(n, 50, 10),
    outcome1 = rbinom(n, 1, 0.5),
    outcome2 = rbinom(n, 1, 0.5)
  )

  # bmvb
  result_bmvb <- bmvb(
    data = test_data,
    grp = "treatment",
    grp_a = "placebo",
    grp_b = "drug",
    y_vars = c("outcome1", "outcome2"),
    test = "right_sided",
    n_it = 1000
  )

  # bglm
  result_bglm <- suppressWarnings(
    bglm(
    data = test_data,
    grp = "treatment",
    grp_a = "placebo",
    grp_b = "drug",
    x_var = "age",
    y_vars = c("outcome1", "outcome2"),
    x_method = "Value",
    x_def = 50,
    test = "right_sided",
    n_burn = 50,
    n_it = 500
  )
  )

  # Both should have same test label
  expect_equal(result_bmvb$info$test_label, "P(drug > placebo)")
  expect_equal(result_bglm$info$test_label, "P(drug > placebo)")

  # Test left-sided
  result_bmvb_left <- bmvb(
    data = test_data,
    grp = "treatment",
    grp_a = "placebo",
    grp_b = "drug",
    y_vars = c("outcome1", "outcome2"),
    test = "left_sided",
    n_it = 1000
  )

  expect_equal(result_bmvb_left$info$test_label, "P(placebo > drug)")
})

#### Test 5: Covariate Methods in bglm ####
test_that("different covariate methods work in bglm", {
  set.seed(204)
  n <- 100
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    age = rnorm(n, 50, 10),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  # Empirical method
  result_empirical <- suppressWarnings(
    bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "age",
    y_vars = c("y1", "y2"),
    x_method = "Empirical",
    x_def = c(40, 60),
    n_burn = 50,
    n_it = 500
  )
  )

  # Value method
  result_value <- suppressWarnings(
    bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "age",
    y_vars = c("y1", "y2"),
    x_method = "Value",
    x_def = 50,
    n_burn = 50,
    n_it = 500
  )
  )

  # Analytical method
  result_analytical <- suppressWarnings(
    bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "age",
    y_vars = c("y1", "y2"),
    x_method = "Analytical",
    x_def = c(40, 60),
    n_burn = 50,
    n_it = 500
  )
  )

  # All should work
  expect_s3_class(result_empirical, "bglm")
  expect_s3_class(result_value, "bglm")
  expect_s3_class(result_analytical, "bglm")

  # Should have correct method recorded
  expect_equal(result_empirical$info$marg_method, "Empirical")
  expect_equal(result_value$info$marg_method, "Value")
  expect_equal(result_analytical$info$marg_method, "Analytical")
})

#### Test 6: Output Structure Consistency ####
test_that("all methods return consistent output structure", {
  set.seed(205)
  n <- 80
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    cluster = rep(1:8, each = n/8),
    x = rnorm(n),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  result_bmvb <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 1000
  )

  result_bglm <- bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "x",
    y_vars = c("y1", "y2"),
    x_method = "Value",
    x_def = 0,
    n_burn = 50,
    n_it = 500
  )

  result_bglmm <- suppressWarnings(
    bglmm(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      id_var = "cluster",
      x_var = "x",
      y_vars = c("y1", "y2"),
      x_method = "Value",
      x_def = 0,
      n_burn = 50,
      n_it = 200,
      n_thin = 2,
      return_diagnostics = FALSE
    )
  )

  # All should have these components
  for (result in list(result_bmvb, result_bglm, result_bglmm)) {
    expect_true("estimates" %in% names(result))
    expect_true("sample_sizes" %in% names(result))
    expect_true("delta" %in% names(result))
    expect_true("info" %in% names(result))

    # Estimates should have these
    expect_true("mean_a" %in% names(result$estimates))
    expect_true("mean_b" %in% names(result$estimates))
    expect_true("sd_a" %in% names(result$estimates))
    expect_true("sd_b" %in% names(result$estimates))

    # Delta should have these
    expect_true("mean_delta" %in% names(result$delta))
    expect_true("se_delta" %in% names(result$delta))
    expect_true("pop" %in% names(result$delta))

    # Info should have these
    expect_true("grp_a" %in% names(result$info))
    expect_true("grp_b" %in% names(result$info))
    expect_true("test_label" %in% names(result$info))
  }

  # bglmm should additionally have J
  expect_true("J" %in% names(result_bglmm$sample_sizes))
})

#### Test 7: Error Consistency ####
test_that("all methods give consistent errors for invalid input", {
  test_data <- data.frame(
    grp = c("A", "B", "A", "B"),
    y1 = c(0, 1, 0, 1),
    y2 = c(1, 0, 1, 0)
  )

  # Missing group variable
  expect_error(
    bmvb(test_data, grp = "missing", grp_a = "A", grp_b = "B", y_vars = c("y1", "y2")),
    "not found"
  )

  expect_error(
    bglm(test_data, grp = "missing", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2")),
    "not found"
  )

  # Too many outcomes
  test_data$y3 <- c(0, 1, 1, 0)

  expect_error(
    bmvb(test_data, grp = "grp", grp_a = "A", grp_b = "B",
         y_vars = c("y1", "y2", "y3")),
    "only supported for 2 outcomes"
  )

  expect_error(
    bglm(test_data, grp = "grp", grp_a = "A", grp_b = "B",
         x_var = "x", y_vars = c("y1", "y2", "y3")),
    "only supported for 2 outcomes"
  )
})

#### Test 8: Print Methods Work for All ####
test_that("print methods work for all result objects", {
  set.seed(206)
  n <- 60
  test_data <- data.frame(
    grp = rep(c("A", "B"), each = n/2),
    cluster = rep(1:6, each = n/6),
    x = rnorm(n),
    y1 = rbinom(n, 1, 0.5),
    y2 = rbinom(n, 1, 0.5)
  )

  result_bmvb <- bmvb(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    y_vars = c("y1", "y2"),
    n_it = 1000
  )

  result_bglm <- bglm(
    data = test_data,
    grp = "grp",
    grp_a = "A",
    grp_b = "B",
    x_var = "x",
    y_vars = c("y1", "y2"),
    x_method = "Value",
    x_def = 0,
    n_burn = 50,
    n_it = 500
  )

  result_bglmm <- suppressWarnings(
    bglmm(
      data = test_data,
      grp = "grp",
      grp_a = "A",
      grp_b = "B",
      id_var = "cluster",
      x_var = "x",
      y_vars = c("y1", "y2"),
      x_method = "Value",
      x_def = 0,
      n_burn = 50,
      n_it = 200,
      n_thin = 2,
      return_diagnostics = FALSE
    )
  )

  # All should print without error
  expect_output(print(result_bmvb), "Bernoulli")
  expect_output(print(result_bglm), "Logistic")
  expect_output(print(result_bglmm), "Multilevel")
})

#### Test 9: Credible Intervals Helper ####
test_that("credible interval function works", {
  samples <- rnorm(1000, mean = 0, sd = 1)

  ci <- credible_interval(samples, prob = 0.95)

  expect_length(ci, 2)
  expect_true(ci[1] < ci[2])
  expect_true(ci[1] < 0)  # Mean is 0, so lower bound should be negative
  expect_true(ci[2] > 0)  # Upper bound should be positive

  # 90% CI should be narrower than 95% CI
  ci_90 <- credible_interval(samples, prob = 0.90)
  ci_95 <- credible_interval(samples, prob = 0.95)

  expect_true(ci_90[1] > ci_95[1])
  expect_true(ci_90[2] < ci_95[2])
})

#### Test 10: Workflow Test ####
test_that("complete workflow works", {
  set.seed(207)

  # 1. Generate data
  n <- 100
  test_data <- data.frame(
    treatment = factor(rep(c("placebo", "drug"), each = n/2)),
    cluster = rep(1:10, each = n/10),
    age = rnorm(n, 50, 10),
    response1 = rbinom(n, 1, 0.5),
    response2 = rbinom(n, 1, 0.5)
  )

  # 2. Run analysis
  result <- bmvb(
    data = test_data,
    grp = "treatment",
    grp_a = "placebo",
    grp_b = "drug",
    y_vars = c("response1", "response2"),
    test = "right_sided",
    rule = "All",
    n_it = 1000,
    return_samples = TRUE
  )

  # 3. Print results
  expect_output(print(result))

  # 4. Extract key information
  posterior_prob <- result$delta$pop
  expect_true(is.numeric(posterior_prob))
  expect_true(posterior_prob >= 0 && posterior_prob <= 1)

  # 5. Compute credible intervals
  ci <- apply(result$samples$delta, 2, credible_interval)
  expect_equal(dim(ci), c(2, 2))

  # Workflow should complete without errors
  expect_true(TRUE)
})
