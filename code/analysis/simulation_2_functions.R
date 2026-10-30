# Functions for simulation 2 (replicability of item-dropping decisions).
#
# Sourced by simulation_2_replicability.qmd. Wraps simulation_functions.R rather than
# editing it: the validation caches read by the manuscript embed a hash of that file, so any
# edit to it would (deliberately) break the manuscript render. Everything here reuses the
# validated DGP, population-target, and sample-alpha functions unchanged.
#
# Design: one item pool, two independent samples of the same size (discovery D, replication
# R). See simulation_2_replicability_plan.md for the full specification.

source("simulation_functions.R")


# Analysis ----------------------------------------------------------------

# Extract, from one discovery sample, one replication sample, and their shared known
# loadings, every quantity the Study 2 estimands require. Returns a one-row tibble.
#
# The discovery columns are named exactly as simulation 1's analyse() names the same
# quantities (alpha_full_sample, alpha_dropped_sample, candidate_item, dropped_conditional,
# and the four population columns), so the Study 1 identity check is a direct column-by-
# column comparison. The replication columns carry an _r suffix.
#
# alpha_dropped_r_at_d is the replication sample's alpha with the DISCOVERY sample's item
# removed: the item is held fixed from D and only re-evaluated in R (RQ1). The replication
# sample's own nomination (item_r, dropped_r) is used only for item agreement (RQ2).
analyse_pair <- function(data_d, data_r, lambda) {

  # discovery sample: identical computation to simulation 1's analyse()
  cov_d <- cov(data_d)
  alpha_full_d <- alpha_from_cov(cov_d)
  alphas_removed_d <- alphas_if_each_item_removed_from_cov(cov_d)
  index_d <- which.max(alphas_removed_d)
  alpha_dropped_d <- unname(alphas_removed_d[index_d])

  # replication sample: its own nomination, and the discovery item held fixed
  cov_r <- cov(data_r)
  alpha_full_r <- alpha_from_cov(cov_r)
  alphas_removed_r <- alphas_if_each_item_removed_from_cov(cov_r)
  index_r <- which.max(alphas_removed_r)

  tibble(
    # discovery sample (names match simulation 1's analyse())
    alpha_full_sample    = alpha_full_d,
    alpha_dropped_sample = alpha_dropped_d,
    candidate_item       = names(lambda)[index_d],
    dropped_conditional  = alpha_dropped_d > alpha_full_d,
    # replication sample
    alpha_full_r         = alpha_full_r,
    alpha_dropped_r_at_d = unname(alphas_removed_r[index_d]),
    alpha_dropped_r_own  = unname(alphas_removed_r[index_r]),
    item_r               = names(lambda)[index_r],
    dropped_r            = unname(alphas_removed_r[index_r]) > alpha_full_r,
    # population targets, full pool and pool minus the DISCOVERY item
    alpha_full_pop    = population_alpha_from_loadings(lambda),
    omega_full_pop    = population_omega_from_loadings(lambda),
    alpha_dropped_pop = population_alpha_from_loadings(lambda[-index_d]),
    omega_dropped_pop = population_omega_from_loadings(lambda[-index_d])
  )
}
# # testing
# lambda <- draw_lambda(5)
# analyse_pair(generate_data(lambda, 50), generate_data(lambda, 50), lambda)


# Wrapper: draw loadings, generate the discovery then the replication sample, analyse.
#
# The ORDER of the random draws is load-bearing. Simulation 1's generate_and_analyse() calls
# draw_lambda() and then generate_data() once, and nothing else consumes random numbers. By
# making the same two calls first, in the same order, and drawing the replication sample
# only afterwards, each task's L'Ecuyer stream produces exactly simulation 1's loadings and
# sample before it produces anything new. Run with simulation 1's task grid, shuffle and
# master seed, every row's discovery sample is therefore simulation 1's sample for that row
# (on the same platform; see the reuse check in simulation_2_replicability.qmd). Do not
# reorder these lines or add random draws before the replication sample.
generate_and_analyse_pair <- function(k_indicators, n_participants) {

  lambda <- draw_lambda(k_indicators = k_indicators)

  data_d <- generate_data(lambda = lambda, n_participants = n_participants)
  data_r <- generate_data(lambda = lambda, n_participants = n_participants)

  analyse_pair(data_d = data_d, data_r = data_r, lambda = lambda)
}
# # testing
# generate_and_analyse_pair(k_indicators = 5, n_participants = 50)


# Exchangeable calibration cell: as above, but with every loading fixed at the same value
# instead of drawn. Under exchangeability every item is equally likely to be nominated, so
# item agreement between the two samples must equal 1/k (a bug check with a known answer).
generate_and_analyse_pair_exchangeable <- function(k_indicators, n_participants,
                                                   loading = 0.5) {

  lambda <- setNames(rep(loading, k_indicators), paste0("item", seq_len(k_indicators)))

  data_d <- generate_data(lambda = lambda, n_participants = n_participants)
  data_r <- generate_data(lambda = lambda, n_participants = n_participants)

  analyse_pair(data_d = data_d, data_r = data_r, lambda = lambda)
}
