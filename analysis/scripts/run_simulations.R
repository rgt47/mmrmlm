#!/usr/bin/env Rscript
# Regenerate the Monte Carlo simulation output that report.Rmd reads
# from analysis/data/derived_data/. Run from the repository root:
#
#   Rscript analysis/scripts/run_simulations.R
#
# Per Morris, White, and Crowther (2019) section 4.1, the RNG seed is
# set ONCE here, before either simulation grid runs, and the two
# simulation kernels (R/sim_kernels.R) do not call set.seed()
# themselves.

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(nlme)

here_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) && d != dirname(d)) {
    d <- dirname(d)
  }
  d
}
root <- here_root()
source(file.path(root, "R", "sim_kernels.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260417)

out_dir <- file.path(root, "analysis", "data", "derived_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_sims <- as.integer(
  Sys.getenv("MMRMLM_N_SIMS", unset = "2000")
)
cat("Running simulations with n_sims =", n_sims, "\n")

# ---- Setting 1: two-timepoint ANCOVA vs MMRM -------------------
t0 <- Sys.time()
scenarios <- expand_grid(
  rho = c(0.3, 0.6, 0.9),
  sigma_ratio = c(0.5, 1.0, 2.0)
)

sim_results <- pmap_dfr(
  scenarios,
  function(rho, sigma_ratio) {
    sim_prepost(
      n_per_arm = 30,
      gamma_true = 3,
      sigma1 = 10,
      sigma2 = 10 * sigma_ratio,
      rho = rho,
      n_sims = n_sims
    ) |>
      mutate(
        rho = rho,
        sigma_ratio = sigma_ratio
      )
  }
)
cat(
  "Setting 1 complete:",
  round(as.numeric(Sys.time() - t0, units = "secs"), 1),
  "sec\n"
)
saveRDS(
  sim_results,
  file.path(out_dir, "sim_results_prepost.rds")
)

# ---- Setting 2: three-timepoint summary stat vs random slopes --
t0 <- Sys.time()
three_tp_results <- map_dfr(
  c(0.5, 2, 8),
  function(sb) {
    sim_three_tp(
      n_per_arm = 30,
      gamma_true = 2,
      sigma_b = sb,
      sigma_e = 4,
      n_sims = n_sims
    ) |>
      mutate(
        sigma_b = sb,
        sigma_e = 4,
        ratio = sb / 4
      )
  }
)
cat(
  "Setting 2 complete:",
  round(as.numeric(Sys.time() - t0, units = "secs"), 1),
  "sec\n"
)
saveRDS(
  three_tp_results,
  file.path(out_dir, "sim_results_three_tp.rds")
)

cat("Saved to", out_dir, "\n")
