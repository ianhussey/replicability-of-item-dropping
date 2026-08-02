# Starter for simulation 2 (see "ideas for sim2 re replicability of item dropping
# choices.md"). NOT used by the manuscript or by simulation_1_bias.qmd.
#
# Simulation 1 draws a fresh item pool AND a fresh sample in every replication, so its
# Aim-4 prediction intervals marginalise over two sources of variation at once and cannot
# say how much of the spread is "which scale you have" versus "which sample you drew".
# This script runs the same DGP, strategies, and extracted quantities in a NESTED design
# (S pools per condition, M independent samples per pool), which is also the structure
# sim 2's item-choice replicability question needs: with several samples per pool, both
# the between/within variance decomposition and same-item agreement across samples are
# computable from one run. In particular it can decompose sim 1's marginal harm rate:
# ~50% harm is compatible both with "every short scale is a coin flip" and with a mixture
# of scales where dropping is near-certain to harm (no weak item) and near-certain to
# help (one genuinely weak item) -- opposite practical advice, indistinguishable in the
# marginal design.
#
# Design choices (all provisional; revisit against the sim-2 ADEMP spec):
#   - Conditions: k = 5 at N = 20 (sim 1's headline harm-rate cell) and N = 1000 (its
#     large-sample contrast). Extend to the full grid for the real sim-2 run.
#   - S = 500 pools x M = 20 samples = 10,000 replications per condition, the same Monte
#     Carlo budget per condition as the main simulation. M = 20 gives per-pool harm-rate
#     estimates with SE <= 0.11, coarse per pool but ample for the distribution ACROSS
#     pools and for the ANOVA variance decomposition, both of which pool over S = 500.
#   - Loadings still come from the same marginal scaled Beta (draw_lambda()), NOT from a
#     hierarchical scale-level model, so pools differ only through their k draws from the
#     marginal. Because the marginal fit ignores between-scale clustering of loadings,
#     the between-pool variance component estimated here is, if anything, an
#     underestimate. See the manuscript's Limitations.
#
# Output: results/simulation_1_supplement_nested_S500_M20.rds -- one row per
# (condition, pool, sample) with the same extracted quantities as the main simulation,
# plus pool and sample ids. Summaries are computed in the manuscript at render time.

library(furrr)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

source("simulation_functions.R")

n_workers <- max(1L, unname(parallelly::availableCores()))
plan(multisession, workers = n_workers)
message(sprintf("Parallel plan: %d workers.", n_workers))

S <- 500   # pools per condition
M <- 20    # samples per pool

nested_path <- sprintf("results/simulation_1_supplement_nested_S%d_M%d.rds", S, M)

# one task = one pool: draw lambda once, then analyse M independent samples from it
one_pool <- function(k_indicators, n_participants, pool) {
  lambda <- draw_lambda(k_indicators = k_indicators)
  map_dfr(seq_len(M), function(s) {
    analyse(data = generate_data(lambda = lambda, n_participants = n_participants),
            lambda = lambda) |>
      mutate(sample = s)
  }) |>
    mutate(k_indicators = k_indicators, n_participants = n_participants, pool = pool)
}

if (file.exists(nested_path)) {

  message("cache exists, nothing to do: ", nested_path)

} else {

  tasks <- expand_grid(k_indicators = 5,
                       n_participants = c(20, 1000),
                       pool = seq_len(S))

  # shuffle for load balancing (N = 1000 pools cost ~50x an N = 20 pool); fixed seed so
  # the run is reproducible. furrr assigns one L'Ecuyer stream per task, so results are
  # identical regardless of worker count.
  set.seed(101)
  tasks <- tasks[sample(nrow(tasks)), ]

  set.seed(42)
  t0 <- Sys.time()
  nested <- tasks |>
    future_pmap(.f = one_pool,
                .progress = TRUE,
                .options = furrr_options(seed = TRUE)) |>
    bind_rows() |>
    arrange(n_participants, pool, sample)
  message(sprintf("run time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  dir.create("results", showWarnings = FALSE)
  write_rds(nested, nested_path, compress = "gz")
  message("wrote ", nested_path)
}
