# Item-dropping to improve in-sample Cronbach’s alpha can harm true reliability: A Monte Carlo simulation

Drafts of [Hussey et al. (2025)](https://doi.org/10.1177/25152459241287123) grappled with the idea that researchers often argue that poorly performing items should be dropped before calculating scale mean scores. I would counter by arguing that the our ability to precisely detect poorly performing items is poor, and that the recommendations of item-drop methods have poor replicability in small samples. It turned that we had misread and miscited some of the studies in this area, and so removed it from the final manuscript. My impression was that there is not much evidence either way regarding the replicability of item dropping in realistic sample sizes. This repo therefore examines this in real and simluated data. This includes a replication of Kopalle and Lehman (1997). 



## Data source

Hussey, I. (2026). granary: A very large item-level open dataset of psychological traits and attitudes (Version 1.0.1) [Data set]. https://doi.org/10.5281/zenodo.21620120 

- All data processing and beta regression fits are done in that repo, not here.



## Repository contents

```
code/analysis/
  simulation_functions.R          Data-generating process, closed-form population targets
                                  (population alpha and omega from loadings), fast sample
                                  alpha and leave-one-out alphas from one covariance
                                  matrix. Sourced by BOTH files below, so the functions
                                  validated are provably the ones the simulation runs.
  validation.qmd / .html          Validates those functions before use: agreement with
                                  psych::alpha(), the closed-form targets against the
                                  model-implied covariance matrix, and that
                                  lavaan::simulateData() actually samples from the
                                  structure the targets assume. Halts if any check fails.
  simulation_1_bias.qmd / .html   The simulation itself, with the full ADEMP
                                  specification. 25 conditions (item-pool size k in
                                  {5, 10, 15, 20, 40} x sample size N in
                                  {20, 30, 50, 200, 1000}) x 10,000 replications =
                                  250,000 generated datasets. Reports every estimand
                                  per condition at more precision, and over more
                                  quantities, than the manuscript has room for.
  simulation_2_nested_starter.R   Scaffold for the planned follow-up (nested design:
                                  several samples per item pool). NOT used by the
                                  manuscript or by simulation_1_bias.qmd.
  results/                        Cached simulation output (.rds). See Reproducibility.
  templates/                      Generic simulation templates, not specific to this study.
  illustrate_cronbachs_alpha.Rmd  Standalone teaching illustration.
  ideas for sim2 ....md           Design notes for the follow-up simulation.

communication/manuscript/
  manuscript.qmd                  The manuscript. Every number, table and figure in its
                                  Results is computed at render time from the caches in
                                  code/analysis/results/ -- nothing is hardcoded.
  manuscript.pdf / .tex           Rendered output.
  bibliography.bib, apa.csl       References and APA 7th citation style.
  _extensions/, elsarticle.cls    Quarto Elsevier template.

data/                             Real item-level datasets, for the planned analyses of
                                  empirical data. The simulation reported in the
                                  manuscript does not read from here: its data-generating
                                  process is parameterised by a loading distribution
                                  fitted in a separate repository (see below).
literature/                       Bibliography files for cited works (PDFs gitignored).
old/                              Superseded work, retained for provenance (gitignored).
```

**External dependency.** The scaled-Beta loading distribution and the item-count
percentiles that define the design grid were fitted in
[`ianhussey/granary`](https://github.com/ianhussey/granary)
(`code/analysis/analysis_factor_loading_distribution.Rmd`), from 3,783 single-factor
loadings across 206 (scale, source) units. Those fitted values enter this repository only
as fixed constants in `draw_lambda()`; the fit is not re-run here.

## Reproducibility

**Software.** R 4.5.2, Quarto 1.8.27, with `lavaan` 0.7.2, `psych` 2.5.3, `furrr` 0.3.1,
`simhelpers` 0.3.1, `dplyr` 1.2.1, `ggplot2` 4.0.3 and `kableExtra` 1.4.1. Full session
information is printed at the end of each rendered `.html`.

**Order of operations.** Render in this order, from within `code/analysis/`:

1. `validation.qmd` — must pass before the simulation is run. It `stop()`s on failure
   rather than warning, so a broken function cannot silently propagate.
2. `simulation_1_bias.qmd` — runs the simulation, or loads it from cache if present.
3. `communication/manuscript/manuscript.qmd` — reads both caches and rebuilds the whole
   Results section.

**Caching.** Every expensive step writes an `.rds` to `code/analysis/results/` and is
skipped if that file already exists, so the manuscript and reports can be re-rendered in
seconds without re-running anything. Delete a cache to force recomputation. The intended
workflow is to run steps 1--2 on an HPC cluster, download `results/`, and re-render
locally.

**Cache keys.** Validation caches embed an md5 hash of `simulation_functions.R` in their
filename (e.g. `..._fnsf2f42899.rds`). Editing those functions changes the hash, so the
old cache no longer matches and the validation re-runs instead of vouching for code that
no longer exists. `manuscript.qmd` recomputes the same hash when locating the validation
caches, and fails to render if they are missing — it cannot quote validation results for
functions that have since changed. Files carrying a different hash (e.g. `..._fns27c7fd3b.rds`)
are orphans from an earlier version of the functions and are safe to delete.

**Reproducible parallelism.** Both `validation.qmd` and `simulation_1_bias.qmd` use
`furrr` with `furrr_options(seed = TRUE)`, which advances an L'Ecuyer-CMRG stream per
task rather than per worker. Results are therefore **identical regardless of how many
cores are used**, so a cluster run and a laptop re-run produce the same numbers. Worker
count is auto-detected with `parallelly::availableCores()`, which honours SLURM
allocations and cgroup quotas; override it via `n_workers_override` at the top of each
file. Note that both files also shuffle their task grid (fixed seed) for load balancing;
this is reproducible but does change which stream each row draws, so results differ from
an unshuffled run.

## License

(c) Ian Hussey (2025)

All code is released under an MIT license. 

All text and images are released under a CC BY 4.0 license. 
