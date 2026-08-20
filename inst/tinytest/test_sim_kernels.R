library(tinytest)

suppressMessages({
  library(nlme)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

pkg_root <- (function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "DESCRIPTION")) &&
    d != dirname(d)) {
    d <- dirname(d)
  }
  d
})()
source(file.path(pkg_root, "R", "sim_kernels.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260819)

# ---- sim_prepost(): ANCOVA/MMRM point-estimate equivalence ----
# This guards the paper's headline claim (Methods "Algebraic
# equivalence", Appendix A): the ANCOVA and MMRM point estimates
# should agree to within Monte Carlo noise on a small run, even
# though this test does not have the power to detect a subtle
# implementation bug in either fitting routine.
res_pp <- sim_prepost(
  n_per_arm = 25,
  gamma_true = 3,
  sigma1 = 10,
  sigma2 = 10,
  rho = 0.6,
  n_sims = 100
)

expect_true(
  all(c("sim", "method", "est", "se", "pval", "ci_lo", "ci_hi") %in%
    names(res_pp)),
  info = "sim_prepost() returns the documented columns"
)

pp_wide <- res_pp |>
  filter(!is.na(est)) |>
  select(sim, method, est) |>
  pivot_wider(names_from = method, values_from = est)

expect_true(
  nrow(pp_wide) >= 90,
  info = "at most a handful of gls() fits fail to converge at n=50"
)

max_gap <- max(abs(pp_wide$ANCOVA - pp_wide$MMRM))
expect_true(
  max_gap < 1e-3,
  info = paste(
    "sim_prepost(): ANCOVA and MMRM point estimates must agree to",
    "gls() numerical convergence tolerance per replication (max gap",
    signif(max_gap, 3), "), per the point-estimate identity",
    "proved in report.Rmd Appendix A. A gap near 1 or larger would",
    "indicate the baseline-constrained MMRM specification has",
    "regressed to the unconstrained cell-means bug this test guards",
    "against (see roxygen note on sim_prepost())."
  )
)

expect_true(
  all(is.finite(res_pp$se[!is.na(res_pp$se)])),
  info = "sim_prepost() standard errors are finite where estimated"
)

expect_true(
  all(res_pp$ci_hi[!is.na(res_pp$ci_lo)] >
    res_pp$ci_lo[!is.na(res_pp$ci_lo)]),
  info = "sim_prepost() confidence intervals have positive width"
)

# ---- sim_three_tp(): output shape and non-degenerate estimates ----
# Guards against the earlier DGM/fitted-model mismatch (random
# intercept + slope fitted against a slope-only DGM), which drove
# non-convergence and NA output silently absorbed by tryCatch.
res_tp <- sim_three_tp(
  n_per_arm = 25,
  gamma_true = 2,
  sigma_b = 2,
  sigma_e = 4,
  n_sims = 100
)

expect_true(
  all(c("sim", "method", "est", "se", "pval", "ci_lo", "ci_hi") %in%
    names(res_tp)),
  info = "sim_three_tp() returns the documented columns"
)

n_na_re <- sum(is.na(res_tp$est[res_tp$method == "Random slopes"]))
expect_true(
  n_na_re <= 5,
  info = paste(
    "at most a handful of lme() fits fail to converge out of 100",
    "replications with the matched slope-only DGM/fit (observed",
    n_na_re, "NA)"
  )
)

re_est <- res_tp$est[
  res_tp$method == "Random slopes" & !is.na(res_tp$est)
]
expect_true(
  length(re_est) > 0 && sd(re_est) > 0,
  info = paste(
    "random-slopes estimates are not degenerate (would indicate a",
    "silently broken fit, as in the pre-remediation",
    "random-intercept-and-slope model fit against a slope-only DGM)"
  )
)

expect_equal(
  sum(res_tp$method == "Summary stat" & is.na(res_tp$est)),
  0L,
  info = "the lm()-based summary-statistic estimator never fails"
)
