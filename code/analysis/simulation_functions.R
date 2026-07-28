# Shared functions for the alpha-if-item-removed simulation study.
#
# Sourced by both:
#   - validation.qmd  (validates the fast alpha functions against psych::alpha)
#   - simulation.qmd  (runs the simulation)
#
# Kept in one file rather than duplicated into each .qmd so that the functions
# validated in the first file are provably the same objects used by the second.

library(lavaan)
library(dplyr)
library(tibble)


# Data-generating process -------------------------------------------------

# draw a named vector of standardized single-factor loadings for one item pool, from a
# scaled Beta fitted to n = 3783 single-factor CFA loadings pooled across trait scales from
# Bainbridge et al. (2022), the AIID battery, and Brysbaert et al. (2024). The shape values
# are the fitted parameters and are not manipulated in this simulation; they are estimated
# elsewhere, in analysis_factor_loading_distribution.Rmd (ported from Hussey's granary repo).
draw_lambda <- function(k_indicators, shape1 = 11.7748, shape2 = 3.8375) {
  lambda <- 2 * rbeta(n = k_indicators, shape1 = shape1, shape2 = shape2) - 1
  names(lambda) <- paste0("item", seq_len(k_indicators))
  return(lambda)
}
# # testing
# draw_lambda(k_indicators = 5)


# build lavaan single-factor model syntax from a named lambda vector.
#
# The residual variances are stated EXPLICITLY as psi_i = 1 - lambda_i^2, and the factor
# variance explicitly as 1. This is essential and not cosmetic: given syntax that specifies
# loadings only, lavaan leaves the residual variances free and simulateData() then defaults
# them to 1, producing items with variance lambda_i^2 + 1 rather than unit variance. The
# closed-form population targets below assume unit item variance (they are written for
# standardized loadings), so with the loadings-only syntax the DGP and the population
# targets describe different populations -- at k = 8, N = 200,000 the generated data's alpha
# was 0.733 against a stated population alpha of 0.812. Stating psi explicitly makes the
# generated covariance structure exactly Cov(X_i, X_j) = lambda_i * lambda_j, Var(X_i) = 1,
# which is the model documented in the DGP section of simulation.qmd.
#
# Numbers are written in fixed notation at full double precision rather than via format()
# (which rounds) or %g (which can emit scientific notation such as "1e-05" for a near-zero
# residual variance, which lavaan's syntax parser does not accept).
build_cfa_model <- function(lambda) {

  as_lavaan_number <- function(x) formatC(x, format = "f", digits = 17)

  loadings <- paste0("factor =~ ",
                     paste0(as_lavaan_number(lambda), "*", names(lambda), collapse = " + "))

  # psi_i = 1 - lambda_i^2, so that Var(X_i) = lambda_i^2 + psi_i = 1
  residual_variances <- paste0(names(lambda), " ~~ ",
                               as_lavaan_number(1 - lambda^2), "*", names(lambda))

  # factor variance fixed at 1 (also simulateData()'s std.lv = TRUE default; stated for clarity)
  factor_variance <- "factor ~~ 1.0*factor"

  paste(c(loadings, residual_variances, factor_variance), collapse = "\n")
}
# # testing
# cat(build_cfa_model(draw_lambda(5)))


# simulate item-level data from a lambda vector and a participant N
generate_data <- function(lambda, n_participants) {
  simulateData(model = build_cfa_model(lambda), sample.nobs = n_participants)
}
# # testing
# generate_data(lambda = draw_lambda(5), n_participants = 50)


# Population targets ------------------------------------------------------

# true population Cronbach's alpha implied by a vector of standardized single-factor
# loadings, under the model used by generate_data() (unit total item variance, error
# variance psi_i = 1 - lambda_i^2, so Cov(item_i, item_j) = lambda_i * lambda_j). Ported
# from alpha_from_loadings() in Hussey's granary repo
# (code/analysis/analysis_factor_loading_distribution.Rmd), which derives this in closed
# form from that covariance structure (Cronbach, 1951).
population_alpha_from_loadings <- function(lambda) {
  k  <- length(lambda)
  S1 <- sum(lambda)
  S2 <- sum(lambda^2)
  (k / (k - 1)) * ((S1^2 - S2) / (k + S1^2 - S2))
}
# # testing
# population_alpha_from_loadings(rep(0.5, 10))


# true reliability (McDonald's omega) of the unit-weighted sum score implied by the same
# loadings: the proportion of sum-score variance due to the common factor,
# omega = S1^2 / (S1^2 + sum(psi_i)) with psi_i = 1 - lambda_i^2. Under this congeneric
# DGP alpha_pop(S) <= omega(S) by Cauchy-Schwarz, with equality only when all loadings in
# S are equal (essential tau-equivalence) -- this gap is estimand 1b and is why the
# reported alpha is evaluated against two distinct targets.
population_omega_from_loadings <- function(lambda) {
  S1 <- sum(lambda)
  S1^2 / (S1^2 + sum(1 - lambda^2))
}
# # testing
# population_omega_from_loadings(rep(0.5, 10))   # equals alpha under equal loadings
# population_omega_from_loadings(c(0.2, 0.5, 0.8, 0.5, 0.6))  # exceeds alpha otherwise


# Sample alpha ------------------------------------------------------------

# Sample (raw) Cronbach's alpha from a covariance matrix:
#   alpha = k/(k-1) * (1 - sum(diag(S)) / sum(S))
# This is the same quantity as psych::alpha(data)$total$raw_alpha (raw = covariance-based,
# not correlation-based) but computed directly, since psych::alpha() builds a large object
# of item statistics that this simulation never uses. Validated against psych::alpha() in
# validation.qmd. Note that, like psych::alpha(check.keys = FALSE), no
# automatic reverse-keying is applied: the fitted loading distribution has a small negative
# tail, and items with negative loadings must be left as-is for the item-dropping rules to
# be able to find and drop them.
alpha_from_cov <- function(S) {
  k <- ncol(S)
  (k / (k - 1)) * (1 - sum(diag(S)) / sum(S))
}
# # testing
# alpha_from_cov(cov(generate_data(draw_lambda(5), 100)))


# Sample alpha of the scale with each item removed in turn: element j is the alpha of the
# scale without item j. Derived from a single covariance matrix rather than by refitting,
# using the fact that deleting item j removes row j and column j, so
#   sum(diag(S_-j)) = sum(diag(S)) - S[j, j]
#   sum(S_-j)       = sum(S) - 2 * colSums(S)[j] + S[j, j]
# (the +S[j, j] corrects for S[j, j] having been subtracted twice). This makes all k
# leave-one-out alphas cost one covariance matrix rather than k calls to psych::alpha().
alphas_if_each_item_removed_from_cov <- function(S) {
  m <- ncol(S) - 1
  diag_S <- diag(S)
  total_removed <- sum(S) - 2 * colSums(S) + diag_S
  diag_removed  <- sum(diag_S) - diag_S
  setNames((m / (m - 1)) * (1 - diag_removed / total_removed), colnames(S))
}
# # testing
# alphas_if_each_item_removed_from_cov(cov(generate_data(draw_lambda(5), 100)))


# Analysis ----------------------------------------------------------------

# Extract, from one generated dataset and its known population loadings, every quantity
# the estimands require (see "Estimands and Targets" and "Methods and extracted
# quantities" in simulation.qmd). Returns a one-row tibble.
#
# Both strategies nominate the same candidate item J* (the argmax of the alpha-if-removed
# table) and differ only in whether acting on it is gated on an observed improvement, so
# J* is found once here and each strategy's outcome is derived downstream.
analyse <- function(data, lambda) {

  covariance_matrix <- cov(data)

  alpha_full_sample <- alpha_from_cov(covariance_matrix)
  alphas_if_removed <- alphas_if_each_item_removed_from_cov(covariance_matrix)

  # the candidate item J*: the item whose removal maximises sample alpha
  candidate_index <- which.max(alphas_if_removed)
  alpha_dropped_sample <- unname(alphas_if_removed[candidate_index])

  results <-
    tibble(# sample estimates
           alpha_full_sample    = alpha_full_sample,
           alpha_dropped_sample = alpha_dropped_sample,
           candidate_item       = names(lambda)[candidate_index],
           # population targets, full pool
           alpha_full_pop = population_alpha_from_loadings(lambda),
           omega_full_pop = population_omega_from_loadings(lambda),
           # population targets, pool minus the candidate item
           alpha_dropped_pop = population_alpha_from_loadings(lambda[-candidate_index]),
           omega_dropped_pop = population_omega_from_loadings(lambda[-candidate_index]),
           # did the conditional strategy act? (the forced strategy always does)
           dropped_conditional = alpha_dropped_sample > alpha_full_sample)

  return(results)
}
# # testing
# lambda <- draw_lambda(5)
# analyse(data = generate_data(lambda = lambda, n_participants = 100), lambda = lambda)


# Wrapper: draw loadings, generate data, analyse, return only the extracted quantities.
# The generated dataset and the loading vector exist only inside this function's frame, so
# neither is ever stored in the simulation data frame -- with 25,000 replications and up to
# 39 items x 1000 participants each, retaining them as list-columns would cost several GB
# of memory and slow every downstream operation.
generate_and_analyse <- function(k_indicators, n_participants) {

  lambda <- draw_lambda(k_indicators = k_indicators)

  data <- generate_data(lambda = lambda,
                        n_participants = n_participants)

  analyse(data = data, lambda = lambda)
}
# # testing
# generate_and_analyse(k_indicators = 5, n_participants = 100)
