## data-raw/generate_examples.R
##
## Run this script ONCE (from the package root) to (re)generate the
## pre-computed example objects that are shipped with the package.
##
##   source("data-raw/generate_examples.R")
##
## The resulting .rda files are placed in data/ and documented in R/data.R.
## They allow print.bglm, print.bglmm, plot.bglm, and plot.bglmm to have
## runnable examples that finish well under 5 seconds.

library(bmco)

## ── bglm_fit ────────────────────────────────────────────────────────────────
## Single-level dataset: 2 groups × 100 subjects, 1 continuous covariate, 2
## binary outcomes.

set.seed(2024)
n <- 200

bglm_data <- data.frame(
  group = rep(c("placebo", "drug"), each = n / 2),
  age   = rnorm(n, mean = 50, sd = 10),
  stringsAsFactors = FALSE
)

bglm_data$y1 <- bglm_data$y2 <- NA_integer_

for (i in seq_len(n)) {
  grpB <- as.integer(bglm_data$group[i] == "drug")
  bglm_data$y1[i] <- rbinom(
    1, 1,
    plogis(-0.50 + 0.75 * grpB + 0.10 * bglm_data$age[i] / 10)
  )
  bglm_data$y2[i] <- rbinom(
    1, 1,
    plogis(-0.50 + 0.80 * grpB + 0.05 * bglm_data$age[i] / 10)
  )
}

bglm_fit <- bglm(
  data    = bglm_data,
  grp     = "group",
  grp_a   = "placebo",
  grp_b   = "drug",
  x_var   = "age",
  y_vars  = c("y1", "y2"),
  x_method = "Empirical",
  x_def   = c(-Inf, Inf),
  test    = "right_sided",
  rule    = "All",
  n_burn  = 10000,
  n_it    = 10000,
  n_chain = 2,
  return_diagnostics = TRUE,
  return_samples     = TRUE
)

## ── bglmm_fit ────────────────────────────────────────────────────────────────
## Multilevel dataset: 20 clusters × 15 subjects, random intercepts.

set.seed(2024)
J  <- 20   # number of clusters
nJ <- 15   # subjects per cluster

uj_1 <- rnorm(J, sd = 0.5)
uj_2 <- rnorm(J, sd = 0.5)

bglmm_data <- data.frame(
  id    = factor(rep(seq_len(J), each = nJ)),
  group = rep(rep(c("placebo", "drug"), each = J / 2), each = nJ),
  age   = rnorm(J * nJ, mean = 50, sd = 10),
  stringsAsFactors = FALSE
)

bglmm_data$y1 <- bglmm_data$y2 <- NA_integer_

for (i in seq_len(J * nJ)) {
  j    <- as.integer(bglmm_data$id[i])
  grpB <- as.integer(bglmm_data$group[i] == "drug")
  bglmm_data$y1[i] <- rbinom(
    1, 1,
    plogis(-0.50 + 0.75 * grpB + 0.10 * bglmm_data$age[i] / 10 + uj_1[j])
  )
  bglmm_data$y2[i] <- rbinom(
    1, 1,
    plogis(-0.50 + 0.80 * grpB + 0.05 * bglmm_data$age[i] / 10 + uj_2[j])
  )
}

bglmm_fit <- bglmm(
  data    = bglmm_data,
  grp     = "group",
  grp_a   = "placebo",
  grp_b   = "drug",
  id_var  = "id",
  x_var   = "age",
  y_vars  = c("y1", "y2"),
  fixed   = c("group", "age", "group_age"),
  random  = "Intercept",
  x_method = "Empirical",
  x_def   = c(-Inf, Inf),
  test    = "right_sided",
  rule    = "All",
  n_burn  = 10000,
  n_it    = 20000,
  n_thin  = 10,
  n_chain = 2,
  return_diagnostics = TRUE,
  return_samples     = TRUE
)
