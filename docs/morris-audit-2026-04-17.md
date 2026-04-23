# Morris et al. (2019) ADEMP Audit: 18-mmrm-vs-lm-simple
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/report/report.Rmd` (the simulation lives entirely here)

Package `R/` and `analysis/scripts/` are empty. Top-level `archive/`
directory is deprecated and excluded from this audit.

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | comparison described in prose |
| DGMs documented | Met | pre-post and three-timepoint DGMs in Rmd chunks |
| Factors varied factorially | Partial | varied per-chunk, not one factorial grid |
| Estimand defined with true value | Met | treatment effect parameterised |
| Methods justified | Met | MMRM vs ANCOVA vs LM compared |
| Performance measures justified | Met | bias, empirical SE, coverage, power |
| n_sim stated | Met | `n_sims = 2000` at `report.Rmd:412` and `:651` |
| n_sim justified via MCSE | Not met | no derivation |
| MCSE reported per metric | Not met | `summary_prepost()` L529-539 returns no MCSE |
| Seed set once | Partial | seeds fixed per function via `seed = 42`, `seed = 123` — not once at program start |
| RNG states stored | Not met | not stored |
| Paired comparisons | Met | methods applied to same dataset per rep |
| Reproducibility | Partial | seeds present as function args; RNGkind not pinned |

## Overall verdict

**Partially compliant.**

## Gaps

- Simulation code lives entirely in the Rmd. Good for traceability,
  bad for reuse, testing, and isolation of seed management.
- `n_sims = 2000` is hardcoded in two places (`:412`, `:651`) with no
  MCSE justification.
- `summary_prepost()` at `report.Rmd:529-539` and `summary_3tp()` at
  `:749-757` compute bias / empirical SE / coverage / power without
  MCSE columns.
- Seed is set per function call via a `seed =` argument, not once at
  program start. Morris §4.1: one seed at the top.
- `RNGkind("L'Ecuyer-CMRG")` not pinned.

## Remediation plan

1. Migrate simulation kernels from `report.Rmd` into `R/` package
   code: `R/dgp.R` for data generation, `R/performance.R` for
   summaries with MCSE columns.
2. Add `mcse_*` columns to `summary_prepost()` and `summary_3tp()`
   covering bias, empirical SE, coverage, and power per Morris Table
   6.
3. Set a single seed at the top of `report.Rmd` (e.g.,
   `set.seed(20260310)`), pin `RNGkind("L'Ecuyer-CMRG")`, and remove
   the per-function `seed =` arguments.
4. Derive `n_sims` from a target MCSE and document. For coverage MCSE
   ≤ 0.5 pp at 0.95, need n ≥ 1900.
5. Store `.Random.seed` per rep in
   `analysis/data/derived_data/rng_states_*.rds`.
6. Add an ADEMP Methods section citing Morris Table 6.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/18-mmrm-vs-lm-simple/mmrmlm/docs/morris-audit-2026-04-17.md*
