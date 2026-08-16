# Publication Review Whitepaper: mmrmlm (18-mmrm-vs-lm-simple)
*Review date: 2026-08-16 10:10 PDT*

Referee-grade review of the research report in this workspace,
conducted to the standards of a statistical journal (Statistics in
Medicine, Biometrics, The American Statistician). Epistemic status is
marked throughout: verified (ran code), inspected (read source),
inferred, or unverified. No manuscript or code files were modified.

## 1. Summary of the work under review

One manuscript exists: `analysis/report/report.Rmd` (rendered PDF at
`analysis/report/report.pdf`; a stamped copy dated 2026-08-15 is
staged in `analysis/report/share/`, confirming the document currently
renders end to end; inspected, not re-rendered here). The paper,
titled "When Do Random Effects Models Reduce to Fixed-Effect
Alternatives? Conditions for Equivalence in Longitudinal Clinical
Trials," reviews the literature on equivalences between joint
longitudinal models and simpler conditional or summary-statistic
analyses (Brogan and Kutner 1980; Lee 1974; Frost, Kenward, and Fox
2008; Liang and Zeger 2000; Liu et al. 2009), then reports two Monte
Carlo studies: Setting 1 compares ANCOVA with an
unstructured-covariance MMRM (fit via `nlme::gls`) in a two-timepoint
design across a 3 x 3 grid of correlation and variance-ratio values;
Setting 2 compares a per-subject mean-change estimator with a
random-slopes LME (fit via `nlme::lme`) in a three-timepoint design
across three slope-variance values. A Discussion organizes findings
into three "conditions" for equivalence. Simulation code lives
entirely inside the Rmd; the package `R/`, `vignettes/`, and
`analysis/scripts/` directories are empty, and the test suite is a
single placeholder assertion (`inst/tinytest/test_basic.R` contains
only `expect_true(TRUE)`; inspected).

## 2. Major issues

1. **Promised third setting (cLDA) is never analyzed.**
   Location: `report.Rmd`, Introduction "Present study" (item 3),
   Discussion "Condition 3," Limitations. The Introduction states
   the paper investigates three canonical settings "through analytic
   argument and Monte Carlo simulation," including "Constrained
   longitudinal data analysis versus ANCOVA." No cLDA model is
   defined in Methods, no cLDA simulation exists in the code, and no
   cLDA results appear (inspected: the only simulation kernels are
   `sim_prepost()` and `sim_three_tp()`). Discussion Condition 3 and
   its "1--4% additional power" figure rest entirely on citations to
   Liu et al. (2009) and Lu et al. (2010), yet are presented in
   parallel with the paper's own simulated conditions. The
   Limitations section compounds this by referring to "the cLDA
   comparison" as if it were performed. A referee will treat this as
   a claim-evidence mismatch severe enough to question the whole
   paper. Remediation: either implement Setting 3 (model definition,
   DGM with and without missing baselines, simulation, table) or
   rewrite the Introduction, Discussion, and Limitations to scope
   the paper to two settings, with cLDA relegated to related work.

2. **The central "algebraic equivalence" is asserted, not proved,
   and the SE-identity claim is wrong as stated.**
   Location: Methods "Algebraic equivalence"; Discussion
   "Condition 1." The manuscript writes the MMRM estimator as
   difference-in-means minus r-hat times baseline imbalance, calls
   it "algebraically identical" to ANCOVA, and dismisses the proof
   in two sentences ("the proof follows from the block structure of
   the GLS normal equations"). No derivation is given, and the
   claimed identity of r-hat with the pooled ANCOVA slope
   S12/S11 requires care about which covariance estimator (ML,
   REML, group-pooled) is used. More seriously, Discussion
   Condition 1 claims "the treatment effect estimate and its
   standard error are algebraically identical to those from
   ANCOVA." This is incorrect as stated: the ANCOVA model-based SE
   is a conditional OLS quantity on n - 3 error degrees of freedom,
   while the `gls` SE is a marginal REML Wald quantity; they agree
   asymptotically and closely in practice but are not algebraically
   identical in finite samples (inferred from standard theory;
   consistent with the paper's own Table 1 showing near, not exact,
   agreement). A statistical journal will not accept "implicit in
   Lee (1974)" for the paper's headline result. Remediation: add a
   short appendix with a complete proof of point-estimate
   equivalence for the complete-data case, state precisely which SE
   quantities coincide and which differ (and at what order), and
   correct Condition 1 to distinguish estimate identity from SE
   near-identity.

3. **Internal contradiction and extrapolation on the Setting 2
   critical ratio.**
   Location: Results, paragraph following Table 2; Discussion
   "Condition 2." The Results paragraph states "the critical ratio
   is approximately sigma_b / sigma approx 4 (equivalently
   2 sigma / sigma_b = 1 for equally spaced designs)"; the
   parenthetical implies a ratio of 2, not 4, so the sentence
   contradicts itself. The Discussion then states the threshold is
   "below approximately 0.25." Three mutually inconsistent numbers
   (4, 2, 0.25) are given for the same quantity (inspected). In
   addition, the simulated grid is sigma_b in {0.5, 2, 8} with
   sigma = 4, i.e., ratios {0.125, 0.5, 2}; a claimed threshold of
   4 lies outside the simulated range and is unsupported.
   Remediation: derive or cite the actual threshold from Frost et
   al. (2008), simulate a grid that brackets it, and make Results
   and Discussion agree.

4. **Setting 2 data-generating model does not match the stated
   model.**
   Location: Methods "Setting 2" versus chunk `sim-three-tp`. The
   Methods specify random intercepts and slopes (a_i, b_i)' ~
   N(0, G) with variances sigma_a^2, sigma_b^2 and covariance
   sigma_ab. The code generates only b_i; there is no random
   intercept (sigma_a = 0) and no intercept-slope covariance
   (inspected). The fitted model `random = ~ time | id`
   nevertheless estimates an intercept variance whose true value is
   zero, placing every fit on the parameter-space boundary. This
   (a) misrepresents the DGM relative to the write-up, (b) invites
   convergence failures that are silently absorbed by `tryCatch`
   and dropped via `filter(!is.na(est))` without any reported
   failure rate, and (c) makes the SE comparison between the LME
   and the summary-statistic estimator hard to interpret, since the
   LME is penalized by estimating a null variance component.
   Remediation: generate a_i (and sigma_ab) per the stated G, or
   restate the DGM as slope-only and fit `random = ~ 0 + time | id`
   as a sensitivity arm; report the non-convergence count per
   scenario in the table (`n_valid` is already computed but never
   shown).

5. **Inference is not compared on equal footing, and Setting 2
   omits promised performance measures.**
   Location: chunks `sim-functions`, `results-prepost`,
   `sim-three-tp`; Methods "ADEMP structure." ANCOVA p-values come
   from the t reference (n - 3 df) while MMRM p-values are computed
   as 2 * pnorm(-|z|); both CIs use +/- 1.96 (inspected). At n = 30
   per arm the normal-versus-t choice alone shifts power and
   coverage by amounts comparable to the method differences under
   study, so the Table 1 power and coverage columns conflate the
   estimator comparison with an arbitrary reference-distribution
   choice; the ANCOVA p-value and its CI are also mutually
   inconsistent (t-based test, z-based interval). Separately, the
   ADEMP section promises "mean estimate, mean model SE, empirical
   SE, power, and coverage" for both settings, but Setting 2
   records no p-values and no CIs, so Table 2 reports neither power
   nor coverage. Remediation: use a common reference distribution
   per method (t with appropriate or Satterthwaite/Kenward-Roger
   df for both, given the small-sample framing that the paper
   itself cites via Kenward and Roger 1997), make each method's CI
   consistent with its test, and add power/coverage (with true
   gamma = 2) to Setting 2 or scope the ADEMP claims down.

6. **Monte Carlo standard errors are computed but never reported.**
   Location: chunks `results-prepost` and `table-three-tp` versus
   Methods ADEMP bullet "each with Monte Carlo SEs per Morris
   Table 6." The `mcse_*` columns exist in `summary_prepost` and
   `summary_3tp` but are dropped by the `select()` that builds each
   kable (inspected). The text then claims "the maximum discrepancy
   in any metric is attributable to Monte Carlo error" without
   displaying the MCSEs that would let a reader check this, and
   without any quantitative criterion. Remediation: report MCSEs in
   the tables (or parenthetically), and support the
   attributable-to-MC-error claim with an explicit comparison of
   observed discrepancies to their MCSEs.

7. **Stale, self-contradictory internal audit material is embedded
   in the manuscript, after the References heading.**
   Location: end of `report.Rmd`, section "Morris et al. (2019)
   ADEMP Compliance," placed after `# References`, so it renders as
   a subsection of the bibliography (inspected). Its listed "key
   gaps" state that seeds are set via `seed = 42` / `seed = 123`
   function arguments and that summaries lack MCSE columns; both
   statements are false for the current code, which pins
   `RNGkind("L'Ecuyer-CMRG")`, sets one seed in setup, and computes
   MCSE columns (inspected; the section reflects
   `docs/morris-audit-2026-04-17.md` before remediation). A
   manuscript that contradicts its own code in its final section is
   not submittable. Remediation: delete the section from the
   manuscript (internal audits belong in `docs/`), or convert the
   accurate parts into a brief reproducibility statement placed
   before References.

8. **The "Headline findings" section is a stub.**
   Location: Results, first subsection. It promises headline
   numbers "extracted via inline R from the `summary_prepost` and
   `summary_3tp` data frames" but contains no numbers and no inline
   R; it only tells the reader to look elsewhere (inspected). A
   referee reads this as scaffolding left in the submission.
   Remediation: either write the section with actual inline-R
   quantities (e.g., maximum absolute estimate discrepancy, SE
   ratio range) or remove it.

9. **Chunk caching undermines the single-seed reproducibility
   design.**
   Location: setup chunk (`cache = TRUE`, `cache.path = "cache/"`)
   together with the Morris single-stream seed discipline. Because
   all simulation chunks draw from one RNG stream set once in
   setup, the realized results depend on which chunks execute in a
   given render; with knitr caching, editing one chunk invalidates
   it but not others, so the RNG state entering a re-executed chunk
   differs from a clean run (inferred from knitr caching semantics;
   the presence of two divergent cache trees, `cache/` at the repo
   root and `analysis/report/cache/`, is consistent with this
   hazard, inspected). Results are therefore not a deterministic
   function of the source. Remediation: give each scenario an
   explicitly derived seed or substream (L'Ecuyer streams support
   this cleanly), or disable caching for simulation chunks; store
   per-run RNG state as the existing audit's remediation plan
   already prescribes.

10. **Missing submission apparatus.**
    Location: whole manuscript. There is no abstract, no keywords,
    no data/code availability statement, no software version or
    session information, no hardware/runtime disclosure, and no
    statement of the non-convergence handling rule (inspected).
    The compendium's reproducibility infrastructure (Dockerfile,
    `renv.lock`) exists but is never referenced in the paper.
    Remediation: add an abstract and keywords; add a reproducibility
    and code-availability paragraph citing the repository, R and
    package versions, and the container; report convergence failure
    counts.

## 3. Minor issues

1. Terminology: in Setting 1 the "MMRM" is fit by `gls` and
   contains no random effects; both comparators are fixed-effects
   models, and the title's "random effects models" therefore
   misdescribes half the paper. Recast as joint (marginal) versus
   conditional models, or motivate the random-intercept
   reparameterization explicitly (inspected).

2. The Frost near-equivalence conditions are quoted as
   "sigma_12 < sigma_1^2 and sigma_1^2 = sigma_2^2," yet the
   simulation varies sigma_2 / sigma_1 in {0.5, 1, 2} and finds
   equivalence throughout; the paper never reconciles the quoted
   homogeneity condition with its own evidence that the estimate
   equivalence does not require it (inspected).

3. Table 2 repeats the same `se_ratio` value on both method rows
   because the join key is `sigma_b` only; restructure so the ratio
   appears once per scenario (inspected).

4. `n_valid` (and thus the number of excluded non-converged fits)
   is computed but not displayed in either table (inspected;
   overlaps issue 4 but applies to Setting 1 as well).

5. Figure 2 uses `scales = "free"` across facets, which visually
   equates panels with very different SE magnitudes; fixed scales
   or an SE-ratio panel would be more honest (inspected).

6. The coverage target is hardcoded as `3` in
   `results-prepost` rather than referencing a named
   `gamma_true` object; fragile under scenario changes
   (inspected).

7. `rm(list = ls())` in the setup chunk is discouraged inside
   knitr documents and is redundant in a fresh render session
   (inspected).

8. US English: "towards" (Setting 2 results paragraph) and the
   `colour =` ggplot argument spelling in prose-adjacent code;
   the argument is valid API but `color =` matches house style
   (inspected).

9. The placeholder test suite (`expect_true(TRUE)`) provides no
   protection for the simulation kernels; at minimum, a test that
   `sim_prepost()` reproduces the ANCOVA/MMRM estimate identity to
   numerical tolerance on a tiny run would guard the paper's core
   claim (inspected).

10. Duplicate figure trees (`figure/`, `report_files/`, and a
    nested `report_files/report_files/`) and two cache trees
    suggest renders from inconsistent working directories; tidy
    and gitignore the byproducts (inspected).

11. The bibliography is adequate and all 23 cited keys resolve in
    `references.bib` (verified by key extraction). Missing
    references a referee may request: Laird (1983) or Holland and
    Rubin on change scores is optional, but Lu, Mehrotra, and Liu
    (2009) sample-size work is already present; add the FDA/EMA
    estimand framing (ICH E9(R1)) if the clinical-trials framing
    is kept, and Twisk et al. (2018) on baseline adjustment
    practice (unverified beyond key presence).

## 4. What remains to be done

Priority order for submission readiness.

(a) Required for correctness:

- Resolve the three-way contradiction on the Setting 2 critical
  ratio (major issue 3) and simulate a grid that brackets the
  corrected threshold.
- Align the Setting 2 DGM with the stated model, or restate the
  model; report non-convergence counts (major issue 4).
- Correct the Condition 1 claim of SE identity; prove
  point-estimate equivalence in an appendix (major issue 2).
- Remove or rewrite the stale ADEMP compliance section and the
  "cLDA comparison" phantom references (major issues 1, 7).
- Make each method's test and CI internally consistent and put
  both methods on a common reference distribution (major issue 5).

(b) Required for acceptance:

- Implement Setting 3 (cLDA) or rescope the paper to two settings
  throughout (major issue 1).
- Report MCSEs in tables and quantify the "attributable to Monte
  Carlo error" claim (major issue 6).
- Replace the stub "Headline findings" section with real inline
  quantities or delete it (major issue 8).
- Fix the seed/caching interaction so results are a deterministic
  function of the source (major issue 9).
- Add abstract, keywords, data/code availability, session and
  hardware information (major issue 10).
- Add power and coverage to Setting 2 per the ADEMP commitments
  (major issue 5).

(c) Desirable polish:

- Migrate simulation kernels to `R/` with tinytest coverage of
  the equivalence identity (minor issue 9; also the standing
  remediation plan in `docs/morris-audit-2026-04-17.md`).
- Table and figure refinements (minor issues 3, 4, 5), terminology
  pass (minor issue 1), US English pass (minor issue 8), byproduct
  cleanup (minor issue 10).
- A small worked real-data or realistic synthetic-trial example to
  ground the pedagogical claims.

## 5. Recommended framing

Plausible framings for this paper:

- (i) New methodological results: novel equivalence conditions
  between mixed and fixed-effect analyses.
- (ii) Methodological synthesis and tutorial: unify known
  equivalences (Lee 1974; Brogan-Kutner 1980; Frost et al. 2008;
  Liang-Zeger 2000; Liu et al. 2009) with confirmatory simulation
  and practical guidance on when the simpler analysis suffices.
- (iii) Software/tools paper: an R package that diagnoses when an
  MMRM reduces to a simpler model (currently only item 4 of Future
  research; no such software exists in the repo).

Recommendation: framing (ii). The literature survey in the paper's
own Introduction establishes that every equivalence it studies is
already published: the two-timepoint ANCOVA/joint-model identity
traces to Lee (1974) and is explicit in Frost et al. (2008); the
summary-statistic connection is Frison-Pocock and Senn; the cLDA
results are Liang-Zeger and Liu et al. Framing (i) would be
rejected for lack of novelty, and framing (iii) has no software to
describe. What the literature does not offer, and practitioners
demonstrably lack, is a single accessible account that states the
conditions side by side, verifies them by simulation with modern
ADEMP reporting, and tells the applied statistician when the
simpler analysis is exactly, approximately, or not at all
equivalent (inferred from the surveyed literature; consistent with
the paper's own "pedagogical and communicative" self-assessment in
the Discussion).

Implications of framing (ii):

- Title: drop the implication of new theory; something like
  "When Is the Mixed Model Doing Nothing More Than ANCOVA?
  A Practical Guide to Equivalences in Longitudinal Trials."
- Abstract and introduction: state explicitly that the
  contributions are synthesis, verification, and guidance; claim
  the decision rules (which condition applies, with thresholds) as
  the deliverable.
- Comparators: keep ANCOVA/MMRM and summary-statistic/LME; add
  cLDA to complete the promised triad, since a synthesis paper
  cannot omit the third strand it surveys; consider adding a
  change-score analysis arm, since that is what practitioners most
  often reach for.
- Emphasize: the three-condition organization of the Discussion
  (its strongest material), a proof appendix for Setting 1, and
  decision guidance with explicit thresholds.
- De-emphasize or move to supplement: the full 3 x 3 Setting 1
  grid (the identity makes most cells redundant; one table plus
  the proof suffices, with the grid as supplementary
  confirmation), per-simulation scatter figures, and all ADEMP
  bookkeeping beyond a compact reproducibility statement.
- Target journals: The American Statistician (Statistical Practice
  or Teacher's Corner) as first choice; Pharmaceutical Statistics
  or the Statistics in Medicine tutorial category as alternatives.
  The current Statistics in Medicine CSL styling is consistent
  with the alternatives.

## 6. Assessment

Verdict: major revision (were this a journal submission, the
realistic outcome is reject-and-resubmit at a methods journal, or
major revision at a practice-oriented outlet under framing (ii)).
The manuscript renders and its two implemented simulations appear
competently built at the code level, with genuine strengths in the
literature survey and the three-condition organization of the
Discussion. However, it promises three settings and delivers two;
its headline equivalence claim is stated without proof and is
partly incorrect as written; the Setting 2 conclusions contain
three mutually inconsistent threshold values, one of which lies
outside the simulated grid; the implemented DGM contradicts the
stated model; inference is compared across methods using different
reference distributions; and the document carries a stale internal
audit section that contradicts the current code. No individual
issue is unfixable, and most remediations are days rather than
months of work, but the paper is not reviewable in its present
state by an external journal.

## 7. Revision history

- 2026-08-16: Initial referee review. No prior whitepaper existed.
  Ten major and eleven minor issues identified; verdict major
  revision; recommended reframing as a methodological synthesis
  and tutorial targeting The American Statistician.
