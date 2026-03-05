#### Unit Tests for bglmm() ####
# Covers all checks in check_input() sections 1-7 and the design-matrix
# creation (section 8) that now lives inside check_input().
#
# All input-validation tests call check_input() directly to avoid running MCMC.

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

# ---------------------------------------------------------------------------
# Shared helper: multilevel dataset with J clusters of n_per_cluster rows.
# Groups alternate within each cluster so every cluster contains both A and B.
# ---------------------------------------------------------------------------
make_bglmm_data <- function(J = 10, n_per_cluster = 20, seed = 200) {
  set.seed(seed)
  n <- J * n_per_cluster
  data.frame(
    id  = factor(rep(seq_len(J), each = n_per_cluster)),
    grp = rep(c("A", "B"), times = n / 2),   # alternates: A,B,A,B,... within every cluster
    x   = rnorm(n, 0, 1),
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
}

# ============================================================================
# TEST GROUP 1: Basic Functionality (structure of validated output)
# ============================================================================

test_that("check_input succeeds with valid bglmm input", {
  test_data <- make_bglmm_data()
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_type(validated, "list")
  expect_true(all(c("y_a", "y_b", "n_a", "n_b", "w",
                    "grp_a", "grp_b", "J", "ml",
                    "fixed", "random", "x_names",
                    "data_list", "x_data", "y_data") %in% names(validated)))
})

test_that("J is correctly determined from id_var", {
  test_data <- make_bglmm_data(J = 8)
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_equal(validated$J, 8)
})

test_that("fixed, random, and x_names have correct values", {
  test_data <- make_bglmm_data()
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_equal(validated$fixed,   c("x", "grp_x"))
  expect_equal(validated$random,  c("Intercept", "grp"))
  expect_equal(validated$x_names, c("x", "grp_x", "Intercept", "grp"))
})

test_that("data_list has J non-empty data frames", {
  J <- 10
  test_data <- make_bglmm_data(J = J)
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_length(validated$data_list, J)
  expect_true(all(sapply(validated$data_list, is.data.frame)))
  expect_true(all(sapply(validated$data_list, nrow) > 0))
})

test_that("x_data is a list of J matrices each with 4 columns", {
  J <- 10
  test_data <- make_bglmm_data(J = J)
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_length(validated$x_data, J)
  expect_true(all(sapply(validated$x_data, is.matrix)))
  expect_true(all(sapply(validated$x_data, ncol) == 4))
})

test_that("y_data is a list of J matrices with 4 columns (Q=4 for K=2)", {
  J <- 10
  test_data <- make_bglmm_data(J = J)
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  expect_length(validated$y_data, J)
  expect_true(all(sapply(validated$y_data, is.matrix)))
  expect_true(all(sapply(validated$y_data, ncol) == 4))
})

test_that("row counts of x_data, y_data, and data_list match per cluster", {
  test_data <- make_bglmm_data(J = 6, n_per_cluster = 10)
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Empirical",
    x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
  )
  for (j in seq_len(validated$J)) {
    expect_equal(nrow(validated$x_data[[j]]),    nrow(validated$data_list[[j]]))
    expect_equal(nrow(validated$y_data[[j]]),    nrow(validated$data_list[[j]]))
  }
})

# ============================================================================
# TEST GROUP 2: General Input Checks (sections 1-3 of check_input)
# ============================================================================

test_that("error when data is not a data frame", {
  expect_error(
    bmco:::check_input(
      data = list(a = 1), grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "data frame"
  )
})

test_that("error for missing group variable", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "missing_grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "not found in data"
  )
})

test_that("error for missing outcome variable", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y_missing"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "not found in data"
  )
})

test_that("error for invalid group value (grp_a)", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "X", grp_b = "B",  # X does not exist
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "not found in grouping variable"
  )
})

test_that("error for more than 2 outcomes", {
  test_data <- make_bglmm_data()
  test_data$y3 <- rbinom(nrow(test_data), 1, 0.5)
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2", "y3"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "only supported for 2 outcomes"
  )
})

test_that("error for non-binary outcome", {
  test_data <- make_bglmm_data()
  test_data$y1 <- sample(0:2, nrow(test_data), replace = TRUE)
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "must be binary"
  )
})

test_that("warning for outcome with no variation", {
  test_data <- make_bglmm_data()
  test_data$y1 <- rep(1, nrow(test_data))
  expect_all_warnings(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    c("no variation.*A", "no variation.*B")
  )
})

test_that("warning for small total sample size in a group", {
  set.seed(210)
  # 5 clusters of 4 obs each; only 1 A per cluster = 5 A total (< 10)
  test_data <- data.frame(
    id  = factor(rep(1:5, each = 4)),
    grp = rep(c("A", "B", "B", "B"), times = 5),   # 1 A and 3 B per cluster
    x   = rnorm(20),
    y1  = rbinom(20, 1, 0.5),
    y2  = rbinom(20, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "Small sample size"
  )
})

# ============================================================================
# TEST GROUP 3: Decision Rule Checks (section 4 of check_input)
# ============================================================================

test_that("error for wrong weight vector length", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "Comp",
      w = c(0.5), analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "not equal to number of outcomes"
  )
})

test_that("error for negative weights", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "Comp",
      w = c(-0.3, 0.7), analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "non-negative"
  )
})

test_that("equal weights used by default for Comp rule", {
  test_data <- make_bglmm_data()
  validated <- suppressMessages(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "Comp",
      w = NULL, analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    )
  )
  expect_equal(validated$w, c(0.5, 0.5))
})

# ============================================================================
# TEST GROUP 4: bglm/bglmm-Specific Checks (section 6 of check_input)
# ============================================================================

test_that("error when x_var is not found in data", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "z_nonexistent", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "not found in data"
  )
})

test_that("error for unknown measurement level (list column as covariate)", {
  test_data <- make_bglmm_data()
  test_data$x <- as.list(test_data$x)
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "unknown measurement level"
  )
})

test_that("error when x_def length 1 but x_method is not 'Value'", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = 0, id_var = "id", n_it = 100, n_thin = 10
    ),
    "x_method"
  )
})

test_that("error when x_def length 2 but x_method is 'Value'", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Value",
      x_def = c(-1, 1), id_var = "id", n_it = 100, n_thin = 10
    ),
    "x_method"
  )
})

test_that("error when x_def range is not in increasing order", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(1, -1), id_var = "id", n_it = 100, n_thin = 10
    ),
    "increasing order"
  )
})

test_that("warning for covariate with no variation", {
  test_data <- make_bglmm_data()
  test_data$x <- rep(0, nrow(test_data))
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "no variation"
  )
})

test_that("warning for high collinearity between group and covariate", {
  # Groups alternate within every cluster so each cluster has both A and B.
  # x is nearly identical to the group indicator, producing r ≈ 0.999 overall.
  set.seed(220)
  n <- 200
  test_data <- make_bglmm_data(J = 10, n_per_cluster = 20, seed = 220)
  grp_numeric <- ifelse(test_data$grp == "B", 1, 0)
  test_data$x <- grp_numeric + rnorm(n, 0, 0.01)
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "High correlation"
  )
})

test_that("warning when x_def lower bound is below observed minimum", {
  test_data <- make_bglmm_data(seed = 230)
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-100, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "below minimum observed value"
  )
})

test_that("warning when x_def upper bound is above observed maximum", {
  test_data <- make_bglmm_data(seed = 231)
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, 100), id_var = "id", n_it = 100, n_thin = 10
    ),
    "above maximum observed value"
  )
})

test_that("warning for few unique values in continuous covariate", {
  set.seed(232)
  n <- 200
  test_data <- data.frame(
    id  = factor(rep(seq_len(10), each = 20)),
    grp = rep(c("A", "B"), times = 100),
    x   = sample(1:3, n, replace = TRUE),   # Only 3 unique values
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "Few unique values"
  )
})
#not run due to computation time
#test_that("message for large n_it without thinning", {
#  test_data <- make_bglmm_data()
#  expect_message(
#    bmco:::check_input(
#      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
#      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
#      analysis = "bglmm", x_var = "x", x_method = "Empirical",
#      x_def = c(-Inf, Inf), id_var = "id", n_it = 1+1e6, n_thin = 1
#    ),
#    "thinning"
#  )
#})

# ============================================================================
# TEST GROUP 5: bglmm-Specific Checks (section 7 of check_input)
# ============================================================================

test_that("error when id_var is NULL", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = NULL, n_it = 100, n_thin = 10
    ),
    "id_var must be specified"
  )
})

test_that("error when id_var is not found in data", {
  test_data <- make_bglmm_data()
  expect_error(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "missing_id", n_it = 100, n_thin = 10
    ),
    "not found in data"
  )
})

test_that("warning for very few clusters (J < 5)", {
  set.seed(240)
  test_data <- data.frame(
    id  = factor(rep(1:4, each = 15)),           # Only 4 clusters
    grp = rep(c("A", "B"), times = 30),
    x   = rnorm(60),
    y1  = rbinom(60, 1, 0.5),
    y2  = rbinom(60, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "Very few clusters"
  )
})

test_that("warning for small minimum cluster size (< 3)", {
  set.seed(241)
  # Cluster 1 has only 2 observations
  test_data <- data.frame(
    id  = factor(c(rep(1, 2), rep(2:10, each = 15))),
    grp = rep(c("A", "B"), length.out = 2 + 9 * 15),
    x   = rnorm(2 + 9 * 15),
    y1  = rbinom(2 + 9 * 15, 1, 0.5),
    y2  = rbinom(2 + 9 * 15, 1, 0.5)
  )
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "very few observations"
  )
})

test_that("warning for clusters with observations from only one group", {
  set.seed(242)
  n <- 10 * 10
  test_data <- data.frame(
    id  = factor(rep(seq_len(10), each = 10)),
    grp = rep(c("A", "B"), times = 50),
    x   = rnorm(n),
    y1  = rbinom(n, 1, 0.5),
    y2  = rbinom(n, 1, 0.5)
  )
  # Force cluster 1 to be all group A
  test_data$grp[test_data$id == "1"] <- "A"
  expect_warning(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "only one group"
  )
})

# ============================================================================
# TEST GROUP 6: Missing Data Handling (section 8 / design-matrix creation)
# ============================================================================

test_that("missing y data reduces rows in data_list", {
  test_data <- make_bglmm_data(J = 8, n_per_cluster = 20)
  test_data$y1[c(1, 5, 21)] <- NA   # A few missing y values

  validated <- suppressWarnings(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    )
  )
  total_rows <- sum(sapply(validated$data_list, nrow))
  expect_true(total_rows < nrow(test_data))
})

test_that("missing x data triggers warning and rows are dropped from data_list", {
  test_data <- make_bglmm_data(J = 8, n_per_cluster = 20)
  test_data$x[c(2, 40, 100)] <- NA   # Three missing covariate values

  expect_warning(
    validated <- bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    ),
    "Missing data"
  )
  total_rows <- sum(sapply(validated$data_list, nrow))
  expect_true(total_rows < nrow(test_data))
})

test_that("x_data and y_data rows match data_list per cluster after NA removal", {
  test_data <- make_bglmm_data(J = 6, n_per_cluster = 15)
  test_data$x[c(3, 20)]  <- NA
  test_data$y2[c(7, 35)] <- NA

  validated <- suppressWarnings(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    )
  )
  for (j in seq_len(validated$J)) {
    expect_equal(nrow(validated$x_data[[j]]), nrow(validated$data_list[[j]]))
    expect_equal(nrow(validated$y_data[[j]]), nrow(validated$data_list[[j]]))
  }
})

test_that("n_a and n_b reflect joint complete-case filtering; y_a and y_b are consistent", {
  test_data <- make_bglmm_data(J = 10, n_per_cluster = 20)
  test_data$x[1:10] <- NA    # 10 rows with missing x

  validated <- suppressWarnings(
    bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Empirical",
      x_def = c(-Inf, Inf), id_var = "id", n_it = 100, n_thin = 10
    )
  )
  expect_true(validated$n_a + validated$n_b < nrow(test_data))
  expect_equal(nrow(validated$y_a), validated$n_a)
  expect_equal(nrow(validated$y_b), validated$n_b)
})

# ============================================================================
# TEST GROUP 7: Discrete (Factor) Covariate
# ============================================================================

test_that("factor covariate is accepted and measurement level is 'discrete'", {
  set.seed(250)
  test_data <- make_bglmm_data()
  test_data$x <- factor(sample(c("low", "high"), nrow(test_data), replace = TRUE))
  validated <- bmco:::check_input(
    data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
    y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
    analysis = "bglmm", x_var = "x", x_method = "Value", x_def = 0,
    id_var = "id", n_it = 100, n_thin = 10
  )
  expect_equal(validated$ml, "discrete")
  expect_length(validated$x_data, validated$J)
})

test_that("error for discrete covariate with more than 2 categories", {
  set.seed(150)
  test_data <- make_bglmm_data()
  test_data$x <- factor(sample(c("low", "medium", "high"), nrow(test_data), replace = TRUE))

  expect_error(
     validated <- bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Value", x_def = 0,
      id_var = "id", n_it = 100, n_thin = 10
    ),
    "only binary covariates.*2 categories.*supported"
  )
})

test_that("error for discrete covariate with only 1 category", {
  set.seed(150)
  test_data <- make_bglmm_data()
  test_data$x <- factor(sample(c("low"), nrow(test_data), replace = TRUE))

  expect_error(
    validated <- bmco:::check_input(
      data = test_data, grp = "grp", grp_a = "A", grp_b = "B",
      y_vars = c("y1", "y2"), test = "right_sided", rule = "All", w = NULL,
      analysis = "bglmm", x_var = "x", x_method = "Value", x_def = 0,
      id_var = "id", n_it = 100, n_thin = 10
    ),
  "only 1 category.*At least 2 categories"
  )
})
