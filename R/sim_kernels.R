#' Simulate the two-timepoint (pre-post) ANCOVA vs MMRM comparison
#'
#' Generates `n_sims` replications of a balanced two-group trial
#' measured at baseline and follow-up, and fits both an ANCOVA (`lm`)
#' and an unstructured-covariance MMRM (`gls`) to each replication.
#'
#' **MMRM specification (baseline-constrained, per Lee 1974 and
#' Frost, Kenward, and Fox 2008, eq. 8).** The fitted MMRM does NOT
#' include a treatment main effect at baseline; treatment enters
#' only through an interaction with the follow-up visit indicator
#' (`grp_post <- group * (time_f == "1")`), so the baseline mean is
#' constrained equal across arms by construction. This matches
#' standard clinical-trial MMRM practice (no true treatment effect
#' is possible at a pre-randomization baseline visit) and is the
#' specification under which the GLS point estimate of the
#' follow-up treatment contrast is algebraically identical to the
#' ANCOVA estimate (see report.Rmd Appendix A for the proof).
#'
#' An earlier version of this function instead fit a fully
#' unconstrained "cell-means" MMRM (`y ~ 0 + tg`, four free means,
#' one per time-by-group cell) and extracted the follow-up contrast
#' `mu_11 - mu_10`. For a saturated per-group mean model, GLS is
#' algebraically identical to the sample mean vector regardless of
#' the covariance structure (a standard SUR/growth-curve result),
#' so that contrast reduces to the RAW, baseline-unadjusted
#' follow-up group difference and is NOT algebraically identical to
#' ANCOVA (confirmed both by direct computation and by derivation;
#' see `inst/tinytest/test_sim_kernels.R`). That version is
#' corrected here.
#'
#' Per Morris, White, and Crowther (2019) section 4.1, the RNG seed
#' is set ONCE by the caller (`RNGkind("L'Ecuyer-CMRG")` plus
#' `set.seed()` at the top of the calling script or document). This
#' function does not call `set.seed()` itself, so results depend on
#' the incoming RNG stream and the order in which scenarios are run.
#'
#' Both methods report Wald z-based p-values and 95% confidence
#' intervals (`est +/- 1.96 * se`), so the two methods are compared
#' on a common reference distribution. This is a disclosed design
#' choice: at n = 30 per arm the ANCOVA t-reference (57 df) and the
#' normal reference are numerically close, and `gls` does not supply
#' a small-sample denominator degrees-of-freedom correction (see
#' report.Rmd Limitations).
#'
#' @param n_per_arm Integer. Sample size per treatment arm.
#' @param gamma_true Numeric. True treatment effect at follow-up.
#' @param sigma1 Numeric. Baseline SD.
#' @param sigma2 Numeric. Follow-up SD.
#' @param rho Numeric. Baseline-follow-up correlation.
#' @param n_sims Integer. Number of Monte Carlo replications.
#'
#' @return A tibble with one row per (method, replication) with
#'   columns `sim`, `method`, `est`, `se`, `pval`, `ci_lo`, `ci_hi`.
#'   Rows for replications where the `gls` fit failed contain `NA`
#'   for the MMRM columns.
sim_prepost <- function(n_per_arm, gamma_true,
                         sigma1, sigma2, rho,
                         n_sims = 2000) {
  sigma12 <- rho * sigma1 * sigma2
  Sigma <- matrix(
    c(sigma1^2, sigma12, sigma12, sigma2^2),
    nrow = 2
  )
  L <- chol(Sigma)
  n <- 2 * n_per_arm
  g <- rep(0:1, each = n_per_arm)

  ancova_est <- ancova_se <- ancova_p <-
    numeric(n_sims)
  mmrm_est <- mmrm_se <- mmrm_p <-
    numeric(n_sims)

  for (s in seq_len(n_sims)) {
    Z <- matrix(rnorm(n * 2), ncol = 2) %*% L
    y0 <- 50 + Z[, 1]
    y1 <- 50 + gamma_true * g + Z[, 2]

    # -- ANCOVA --
    fit_a <- lm(y1 ~ factor(g) + y0)
    cf_a <- summary(fit_a)$coefficients
    ancova_est[s] <- cf_a[2, 1]
    ancova_se[s] <- cf_a[2, 2]
    # Wald z-based p-value, matching the MMRM reference
    # distribution below (common-footing fix; see roxygen note).
    ancova_p[s] <- 2 * pnorm(-abs(ancova_est[s] / ancova_se[s]))

    # -- MMRM (baseline-constrained, unstructured cov) --
    dat <- data.frame(
      id = factor(rep(seq_len(n), 2)),
      time_f = factor(rep(c(0, 1), each = n)),
      time_num = rep(c(0, 1), each = n),
      group = rep(g, 2),
      y = c(y0, y1)
    )
    dat$grp_post <- dat$group * (dat$time_f == "1")

    fit_m <- tryCatch(
      gls(y ~ time_f + grp_post,
        data = dat,
        correlation = corSymm(form = ~ 1 | id),
        weights = varIdent(form = ~ 1 | time_num)
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_m)) {
      cf_m <- summary(fit_m)$tTable
      mmrm_est[s] <- cf_m["grp_post", "Value"]
      mmrm_se[s] <- cf_m["grp_post", "Std.Error"]
      z <- mmrm_est[s] / mmrm_se[s]
      mmrm_p[s] <- 2 * pnorm(-abs(z))
    } else {
      mmrm_est[s] <- NA
      mmrm_se[s] <- NA
      mmrm_p[s] <- NA
    }
  }

  tibble::tibble(
    sim = rep(seq_len(n_sims), 2),
    method = rep(c("ANCOVA", "MMRM"),
      each = n_sims),
    est = c(ancova_est, mmrm_est),
    se = c(ancova_se, mmrm_se),
    pval = c(ancova_p, mmrm_p)
  ) |>
    dplyr::mutate(
      ci_lo = est - 1.96 * se,
      ci_hi = est + 1.96 * se
    )
}

#' Simulate the three-timepoint summary-statistic vs random-slopes
#' comparison
#'
#' Generates `n_sims` replications of a balanced two-group trial
#' measured at three evenly spaced timepoints under a **slope-only**
#' random-effects data-generating model (random slope `b_i`, no
#' random intercept), and compares a per-subject mean-change
#' estimator (`lm` on `y[last] - y[first]` averaged, via consecutive
#' differences) with a random-slopes LME fit as
#' `lme(y ~ time * group, random = ~ 0 + time | id)`.
#'
#' The DGM and the fitted random-effects structure are matched
#' deliberately (slope-only, `random = ~ 0 + time | id`): earlier
#' versions of this simulation generated data with no random
#' intercept but fit a model with both a random intercept and a
#' random slope (`random = ~ time | id`), which placed every fit on
#' the boundary of the intercept-variance parameter space. See
#' report.Rmd, Methods "Setting 2", for the disclosed rationale.
#'
#' Per Morris, White, and Crowther (2019) section 4.1, the RNG seed
#' is set ONCE by the caller; this function does not call
#' `set.seed()`.
#'
#' @param n_per_arm Integer. Sample size per treatment arm.
#' @param gamma_true Numeric. True treatment effect on the rate of
#'   change (slope contrast).
#' @param sigma_b Numeric. Between-subject SD of the random slope.
#' @param sigma_e Numeric. Residual (measurement) SD.
#' @param n_sims Integer. Number of Monte Carlo replications.
#'
#' @return A tibble with one row per (method, replication) with
#'   columns `sim`, `method`, `est`, `se`, `pval`, `ci_lo`, `ci_hi`.
#'   Rows for replications where the `lme` fit failed contain `NA`
#'   for the random-slopes columns.
sim_three_tp <- function(n_per_arm, gamma_true,
                          sigma_b, sigma_e,
                          n_sims = 2000) {
  n <- 2 * n_per_arm
  g <- rep(0:1, each = n_per_arm)
  t_vec <- c(0, 1, 2)
  fixed_est <- fixed_se <- fixed_p <- numeric(n_sims)
  re_est <- re_se <- re_p <- numeric(n_sims)

  for (s in seq_len(n_sims)) {
    b_i <- rnorm(n, 0, sigma_b)

    dat <- tidyr::expand_grid(
      subj = seq_len(n),
      time = t_vec
    ) |>
      dplyr::mutate(
        id = factor(subj),
        group = factor(g[subj]),
        b = b_i[subj],
        y = 50 +
          (2 + gamma_true *
            as.numeric(group == 1)) * time +
          b * time +
          rnorm(dplyr::n(), 0, sigma_e)
      )

    # -- Summary statistic: mean change --
    c_dat <- dat |>
      dplyr::arrange(subj, time) |>
      dplyr::group_by(subj) |>
      dplyr::mutate(change = y - dplyr::lag(y)) |>
      dplyr::filter(!is.na(change)) |>
      dplyr::summarise(
        mean_c = mean(change),
        group = dplyr::first(group),
        .groups = "drop"
      )

    fit_f <- lm(mean_c ~ group, data = c_dat)
    cf_f <- summary(fit_f)$coefficients
    fixed_est[s] <- cf_f[2, 1]
    fixed_se[s] <- cf_f[2, 2]
    fixed_p[s] <- 2 * pnorm(-abs(fixed_est[s] / fixed_se[s]))

    # -- Random slopes model (slope-only random effect, matching
    #    the DGM: no random intercept) --
    fit_r <- tryCatch(
      lme(y ~ time * group,
        random = ~ 0 + time | id,
        data = dat,
        control = lmeControl(
          opt = "optim",
          maxIter = 200
        )),
      error = function(e) NULL
    )

    if (!is.null(fit_r)) {
      cf_r <- summary(fit_r)$tTable
      rn <- grep("time:group", rownames(cf_r))
      re_est[s] <- cf_r[rn, "Value"]
      re_se[s] <- cf_r[rn, "Std.Error"]
      re_p[s] <- 2 * pnorm(-abs(re_est[s] / re_se[s]))
    } else {
      re_est[s] <- NA
      re_se[s] <- NA
      re_p[s] <- NA
    }
  }

  tibble::tibble(
    sim = rep(seq_len(n_sims), 2),
    method = rep(
      c("Summary stat", "Random slopes"),
      each = n_sims
    ),
    est = c(fixed_est, re_est),
    se = c(fixed_se, re_se),
    pval = c(fixed_p, re_p)
  ) |>
    dplyr::mutate(
      ci_lo = est - 1.96 * se,
      ci_hi = est + 1.96 * se
    )
}
