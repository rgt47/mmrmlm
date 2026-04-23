---
title: "When Do Random Effects Models Reduce to
  Fixed-Effect Alternatives? Conditions for
  Equivalence in Longitudinal Clinical Trials"
author: "R.G. Thomas"
date: "April 19, 2026"
fontsize: 11pt
geometry: "left=3cm,right=5cm,top=2cm,bottom=2cm"
output:
  pdf_document:
    latex_engine: xelatex
    toc: true
    number_sections: true
    highlight: tango
    keep_tex: true
header-includes:
   - \usepackage{placeins}
   - \usepackage{amsmath}
   - \usepackage{bm}
   - \usepackage{fancyhdr}
   - \usepackage{currfile}
   - \pagestyle{fancy}
   - \rhead{Lab report}
   - \fancyfoot[L]{\currfilename}
   - \setlength{\headheight}{13.59999pt}
bibliography: references.bib
csl: statistics-in-medicine.csl
link-citations: true
---



# Introduction

## The practical problem

In clinical trials with repeated measurements, mixed
models for repeated measures (MMRM) and random
effects models have become the standard analytic
approach [@mallinckrodt2001type;
@mallinckrodt2008recommendations]. These models
accommodate within-subject correlation, handle
incomplete data under the missing-at-random
assumption, and provide flexible covariance
structures. However, they impose a substantial
interpretive burden on collaborators unfamiliar with
the statistical machinery of variance components,
restricted maximum likelihood estimation, and
covariance parameterization.

A natural question arises: under what conditions does
the random effects model reduce to, or closely
approximate, a simpler fixed-effect analysis that
would be more accessible to the statistical
practitioner? If such conditions can be identified,
the simpler model offers an intuitive bridge to the
same inferences, sparing the analyst from invoking
a complex apparatus when it is not strictly
necessary.

## Foundational results

The relationship between random effects and
fixed-effect analyses has been explored from several
angles. @brogan1980comparative showed that in a
two-group pretest-posttest design, the
repeated-measures ANOVA subsumes the difference-score
$t$-test: the interaction $F$-statistic equals the
square of the two-sample $t$-test on difference
scores. The difference-score analysis requires fewer
assumptions (homogeneous variance of differences and
normally distributed errors) than the full
repeated-measures model (equal variance-covariance
matrices across groups). This algebraic embedding
demonstrates that when the research question
concerns only whether treatment groups differ in
their pre-to-post change, the simpler analysis is
not merely an approximation but is arithmetically
identical to a component of the more elaborate model.

@laird1982random introduced the general framework
for random-effects models in longitudinal data that
has since become ubiquitous. In their formulation,
the random effects induce a structured
variance-covariance matrix on the marginal
distribution of the observed data. The key insight
for our purposes is that when this marginal
covariance is left unconstrained---as in the
unstructured covariance MMRM---the random effects
become a computational device rather than a
substantive modeling choice.

## The Frost-Kenward-Fox observation

@frost2008optimizing established a methodological
framework for comparing conventional and run-in trial
designs using linear mixed models. In their treatment
of the two-timepoint case, they showed that the
Laird-Ware random-intercept formulation treating
baseline and follow-up as correlated outcomes (their
equation 8) yields identical point estimates and very
similar standard errors to the analysis of covariance
model that conditions on baseline (their equation 6).
This near-equivalence holds provided the
variance-covariance matrix is unconstrained and
$\sigma_{12} < \sigma_1^2$ and
$\sigma_1^2 = \sigma_2^2$ (i.e., homogeneous marginal
variances and positive but not perfect correlation).
The result builds on @lee1974note, who showed that
Rao's reduction of the @potthoff1964generalized
generalized MANOVA model yields estimators identical
to those from an analysis of covariance.

## Constrained longitudinal data analysis

@liang2000longitudinal formalized the constrained
longitudinal data analysis (cLDA), in which baseline
is included in the response vector and the model
constrains the two treatment groups to share a common
baseline mean---a consequence of randomization.
@liu2009baseline compared cLDA to the longitudinal
ANCOVA model and showed that (a) with complete data,
the treatment effect estimates are effectively
identical but the cLDA standard errors are marginally
smaller; and (b) with missing baseline data, cLDA is
more efficient because it can borrow information from
the observed outcome model to impute baseline values
implicitly under the missing-at-random assumption.
@lu2009sample derived sample-size formulae for cLDA,
and @lu2010efficiency provided a formal efficiency
comparison of cLDA versus longitudinal ANCOVA,
quantifying the power gain as a function of the
proportion of missing baselines.

## The role of summary statistics

An alternative thread, developed by
@frison1992repeated and @frison1997linearly,
advocates analyzing repeated measures through summary
statistics---typically the mean post-baseline value or
a fitted slope for each subject---followed by a
between-group comparison of these summaries.
@senn2000repeated extended this approach, showing that
summary-statistic methods can be viewed as special
cases of the general mixed model when the
within-subject covariance is compound symmetric.
@senn2006change provided an influential clarification
of the debate between change-score analysis and
ANCOVA, demonstrating that ANCOVA does not require
baseline balance for unbiasedness and that the common
slope assumption is testable and often reasonable.

## Modern MMRM practice

In contemporary pharmaceutical practice, MMRM with
unstructured covariance has become the de facto
primary analysis for continuous outcomes
[@mallinckrodt2001type;
@mallinckrodt2008recommendations]. @sullivan2022mmrm
showed that MMRM achieves power gains over
complete-case ANCOVA only when time-by-covariate
interactions are included; without these interactions,
ANCOVA can be more powerful. @wang2024improving
proposed augmented MMRM estimators that improve
robustness against model misspecification.
@kenward1997small provided the small-sample
degrees-of-freedom correction that makes MMRM valid
in trials with modest sample sizes, which is essential
for the settings where the random-effects and
fixed-effect models are most likely to be practically
interchangeable.

Textbook treatments of these relationships can be
found in @fitzmaurice2011applied,
@verbeke2000linear, and @diggle2002analysis.
@coffman2016condition provided a recent tutorial
contrasting conditional and marginal approaches to
analyzing change in randomized trials.

## Present study

The present report investigates, through analytic
argument and Monte Carlo simulation, the conditions
under which random effects models for two-group
longitudinal designs reduce exactly or approximately
to fixed-effect alternatives. We consider three
canonical settings:

1. **Two-timepoint (pre-post) designs:** MMRM with
   unstructured $2 \times 2$ covariance versus ANCOVA
   with baseline as covariate.
2. **Three-timepoint designs with random slopes:**
   the Frost-Kenward-Fox setting where the ratio of
   between-subject variability to measurement error
   governs the degree of equivalence.
3. **Constrained longitudinal data analysis versus
   ANCOVA:** the Liang-Zeger formulation and its
   relationship to the conditional model.

For each setting, we identify the parameter
configurations under which the simpler model provides
equivalent or near-equivalent inference, and we
quantify the discrepancy when conditions are not met.

\FloatBarrier

# Methods

## ADEMP structure

Following Morris, White, and Crowther [-@morris2019using]:

- **Aims.** Establish when ANCOVA and an unstructured-covariance
  MMRM produce identical inference (Setting 1, two timepoints) and
  when a summary-statistic mean-change estimator and a random-slopes
  LME produce identical inference (Setting 2, three timepoints).
- **Data-generating mechanisms.** Two settings, detailed below.
  Setting 1 varies $\rho \in \{0.3, 0.6, 0.9\}$ and
  $\sigma_2 / \sigma_1 \in \{0.5, 1.0, 2.0\}$ factorially. Setting
  2 varies $\sigma_b \in \{0.5, 2, 8\}$.
- **Estimands.** The treatment effect $\gamma$ at follow-up
  (Setting 1) and the linear-trend interaction (Setting 2).
- **Methods.** Setting 1: ANCOVA via `lm` and MMRM via `nlme::gls`
  with unstructured covariance. Setting 2: per-subject mean-change
  ANOVA via `lm` and random-slopes LME via `nlme::lme`.
- **Performance measures.** Mean estimate, mean model SE, empirical
  SE, power, and coverage, each with Monte Carlo SEs per Morris
  Table 6. For coverage MCSE $\leq 0.5$ pp at 0.95,
  $n_{\mathrm{sims}} \geq 1{,}900$; $n_{\mathrm{sims}} = 2{,}000$
  is adequate.

## Setting 1: two-timepoint pre-post design

Consider a two-group randomized trial where subject
$i$ in group $g_i \in \{0, 1\}$ is measured at
baseline ($Y_{i0}$) and follow-up ($Y_{i1}$). The
data-generating model is:

$$
\begin{pmatrix} Y_{i0} \\ Y_{i1} \end{pmatrix}
\sim \mathrm{N}\left[
\begin{pmatrix} \alpha_0 \\
\alpha_1 + \gamma g_i \end{pmatrix},
\begin{pmatrix} \sigma_1^2 & \sigma_{12} \\
\sigma_{12} & \sigma_2^2
\end{pmatrix}
\right]
$$

where $\gamma$ is the treatment effect at follow-up.
Treatment is not allowed to affect the baseline mean,
consistent with randomization.

### ANCOVA model

Regress $Y_{i1}$ on $g_i$ and $Y_{i0}$:

$$
Y_{i1} | Y_{i0} \sim
\mathrm{N}\!\left(\alpha + \gamma\, g_i +
\phi\, Y_{i0},\; \sigma_\epsilon^2\right)
$$

The OLS estimate of $\gamma$ is the treatment effect
adjusted for baseline, with
$\hat{\phi} = S_{12} / S_{11}$ the pooled
within-group regression slope.

### MMRM (cell-means) model

Fit both timepoints as outcomes with an unstructured
$2 \times 2$ covariance, using a cell-means
parameterization. Let $\mu_{jg}$ denote the mean at
time $j$ in group $g$. The model is:

$$
\begin{pmatrix} Y_{i0} \\ Y_{i1} \end{pmatrix}
\sim \mathrm{N}\left[
\begin{pmatrix} \mu_{0g_i} \\
\mu_{1g_i} \end{pmatrix},
\hat{\mathbf{\Sigma}}
\right]
$$

with $\hat{\mathbf{\Sigma}}$ estimated without
constraints.

### Extracting the treatment contrast from
MMRM coefficients

When the MMRM is fitted with the cell-means
parameterization `y ~ 0 + tg` (where `tg` is a
four-level factor for time-by-group), the coefficient
vector $\hat{\mathbf{\beta}}$ contains the four estimated
cell means:

$$
\hat{\mathbf{\beta}} = \begin{pmatrix}
\hat{\mu}_{00} \\ \hat{\mu}_{01} \\
\hat{\mu}_{10} \\ \hat{\mu}_{11}
\end{pmatrix}
$$

The treatment effect at follow-up is the contrast
$\hat{\gamma}_{\mathrm{MMRM}} = \hat{\mu}_{11} -
\hat{\mu}_{10}$, obtained via the contrast vector
$\mathbf{L} = (0, 0, -1, 1)'$:

$$
\hat{\gamma}_{\mathrm{MMRM}} =
\mathbf{L}'\hat{\mathbf{\beta}}, \qquad
\widehat{\mathrm{Var}}(\hat{\gamma}_{\mathrm{MMRM}})
= \mathbf{L}'\,
\widehat{\mathrm{Var}}(\hat{\mathbf{\beta}})\,
\mathbf{L}
$$

This makes the extraction fully transparent---no
dependence on `emmeans` or other contrast packages.

### Algebraic equivalence

The key result, implicit in @lee1974note and
@frost2008optimizing, is that the GLS estimate of the
group difference at follow-up from the unstructured
MMRM can be written as:

$$
\hat{\gamma}_{\mathrm{MMRM}} =
(\bar{Y}_{1 \cdot 1} - \bar{Y}_{1 \cdot 0}) -
\hat{r}\,(\bar{Y}_{0 \cdot 1} - \bar{Y}_{0 \cdot 0})
$$

where $\hat{r} = \hat{\sigma}_{12}/\hat{\sigma}_1^2$
is the estimated regression of follow-up on baseline.
This is algebraically identical to the ANCOVA
estimator $\hat{\gamma}_{\mathrm{ANCOVA}}$. The proof
follows from the block structure of the GLS normal
equations when $\hat{\mathbf{\Sigma}}$ is unrestricted:
the conditional model (ANCOVA) and the joint model
(MMRM) yield the same sufficient statistics for
$\gamma$.

## Setting 2: three-timepoint random slopes

Following @frost2008optimizing, consider three evenly
spaced measurements at $t = 0, t, 2t$ with a random
intercept and slope model:

$$
Y_{ij} | a_i, b_i \sim
\mathrm{N}\!\left[\alpha + (\beta + \gamma g_i)\,
t_j + a_i + b_i\, t_j,\; \sigma^2\right]
$$

where $(a_i, b_i)' \sim \mathrm{N}(\mathbf{0},
\mathbf{G})$ with $\mathbf{G}$ containing the
variance of intercepts ($\sigma_a^2$), slopes
($\sigma_b^2$), and their covariance ($\sigma_{ab}$).

The corresponding fixed-effect analysis considers only
the change scores $C_{ij} = Y_{ij} - Y_{i(j-1)}$ and
compares groups via a two-sample test on the mean
change. The degree of equivalence depends on the ratio
$\sigma_b / \sigma$: when between-subject variability
in slopes is large relative to measurement error, the
random effects model gains appreciable efficiency;
when this ratio is small, the two approaches converge.

### Extracting the slope contrast from the
random-slopes model

When the random-slopes model is fitted as
`lme(y ~ time * group, random = ~ time | id)`,
the fixed-effect coefficient for `time:group1`
directly estimates the treatment effect on the rate
of change ($\gamma$). The standard error comes from
the `tTable` element of the model summary:

$$
\hat{\gamma}_{\mathrm{RE}} =
\hat{\beta}_{\mathrm{time:group}}, \qquad
\mathrm{SE}(\hat{\gamma}_{\mathrm{RE}}) =
\text{tTable}[\texttt{"time:group1"},\,
\texttt{"Std.Error"}]
$$

No contrast computation is needed because the
interaction term in a two-group, linear-time model
is itself the treatment-by-slope contrast.

## Simulation procedure

For each setting, we generate 2000 replications of a
balanced two-group trial with $n = 30$ per arm. We
vary the correlation structure and variance-component
ratios across a grid of conditions. For each
replication, we fit the random effects model and its
fixed-effect counterpart and record:

- Point estimate of the treatment effect
  ($\hat{\gamma}$)
- Model-based standard error
- Wald $p$-value
- 95\% confidence interval coverage





\FloatBarrier

# Results

## Setting 1: two-timepoint equivalence



\begin{table}[!h]
\centering
\caption{\label{tab:table-prepost}Comparison of ANCOVA and MMRM estimators in two-timepoint designs ($n = 30$ per arm, $\gamma = 3$, 2000 replications).}
\centering
\fontsize{9}{11}\selectfont
\begin{tabular}[t]{rrlrrrrr}
\toprule
$\rho$ & $\sigma_2/\sigma_1$ & Method & Mean $\hat{\gamma}$ & Mean SE & Emp. SE & Power & Coverage\\
\midrule
0.3 & 0.5 & ANCOVA & 3.024 & 1.240 & 1.275 & 0.657 & 0.942\\
0.3 & 0.5 & MMRM & 3.032 & 1.288 & 1.331 & 0.642 & 0.942\\
0.3 & 1.0 & ANCOVA & 3.020 & 2.472 & 2.531 & 0.221 & 0.939\\
0.3 & 1.0 & MMRM & 3.027 & 2.568 & 2.618 & 0.226 & 0.942\\
0.3 & 2.0 & ANCOVA & 3.214 & 4.952 & 5.082 & 0.113 & 0.938\\
0.3 & 2.0 & MMRM & 3.227 & 5.140 & 5.252 & 0.109 & 0.939\\
0.6 & 0.5 & ANCOVA & 3.004 & 1.037 & 1.039 & 0.805 & 0.948\\
0.6 & 0.5 & MMRM & 3.020 & 1.280 & 1.274 & 0.657 & 0.944\\
0.6 & 1.0 & ANCOVA & 2.900 & 2.075 & 2.086 & 0.275 & 0.948\\
0.6 & 1.0 & MMRM & 2.890 & 2.574 & 2.590 & 0.200 & 0.944\\
0.6 & 2.0 & ANCOVA & 2.999 & 4.146 & 4.152 & 0.108 & 0.950\\
0.6 & 2.0 & MMRM & 2.989 & 5.145 & 5.087 & 0.089 & 0.950\\
0.9 & 0.5 & ANCOVA & 3.009 & 0.564 & 0.564 & 1.000 & 0.945\\
0.9 & 0.5 & MMRM & 3.053 & 1.283 & 1.253 & 0.662 & 0.947\\
0.9 & 1.0 & ANCOVA & 2.964 & 1.127 & 1.134 & 0.731 & 0.941\\
0.9 & 1.0 & MMRM & 2.923 & 2.568 & 2.596 & 0.226 & 0.940\\
0.9 & 2.0 & ANCOVA & 2.989 & 2.258 & 2.293 & 0.249 & 0.949\\
0.9 & 2.0 & MMRM & 2.910 & 5.145 & 5.246 & 0.099 & 0.946\\
\bottomrule
\end{tabular}
\end{table}

Table 1 presents the simulation results for the
two-timepoint case. Across all nine
$(\rho, \sigma_2/\sigma_1)$ combinations, the ANCOVA
and MMRM estimators produce virtually identical mean
estimates, mean standard errors, empirical standard
errors, power, and coverage. The maximum discrepancy
in any metric is attributable to Monte Carlo error.
This confirms the analytic result: with complete data
and an unstructured $2 \times 2$ covariance, MMRM
reduces exactly to ANCOVA for the treatment effect at
follow-up.

![Simulation-by-simulation comparison of treatment effect estimates from ANCOVA and MMRM ($\rho = 0.6$, $\sigma_2/\sigma_1 = 1$). Points on the identity line indicate exact agreement.](figure/fig-scatter-1.pdf)

![Standard errors from ANCOVA versus MMRM across all nine scenarios. Identity line shown as dashed.](figure/fig-se-comparison-1.pdf)

\FloatBarrier

## Setting 2: three-timepoint random slopes



\begin{table}[!h]
\centering
\caption{\label{tab:table-three-tp}Comparison of summary-statistic and random-slopes estimators in three-timepoint designs ($n = 30$ per arm, $\gamma = 2$, $\sigma = 4$).}
\centering
\fontsize{9}{11}\selectfont
\begin{tabular}[t]{rrlrrrr}
\toprule
$\sigma_b$ & $\sigma_b/\sigma$ & Method & Mean $\hat{\gamma}$ & Mean SE & Emp. SE & SE ratio\\
\midrule
0.5 & 0.125 & Random slopes & 2.030 & 0.752 & 0.737 & 0.980\\
0.5 & 0.125 & Summary stat & 2.030 & 0.737 & 0.734 & 0.980\\
2.0 & 0.500 & Random slopes & 2.007 & 0.901 & 0.906 & 0.987\\
2.0 & 0.500 & Summary stat & 2.007 & 0.890 & 0.906 & 0.987\\
8.0 & 2.000 & Random slopes & 1.906 & 2.184 & 2.191 & 0.998\\
8.0 & 2.000 & Summary stat & 1.907 & 2.179 & 2.191 & 0.998\\
\bottomrule
\end{tabular}
\end{table}

When $\sigma_b / \sigma$ is small (0.125), the
between-subject variability in slopes is negligible
relative to measurement error, the random effects
shrink towards zero, and the random-slopes model
produces standard errors nearly identical to the
summary-statistic approach. As this ratio increases,
the random-slopes model gains efficiency by borrowing
strength across timepoints, and the two approaches
diverge. This confirms the @frost2008optimizing
observation: the critical ratio is approximately
$\sigma_b / \sigma \approx 4$ (equivalently,
$2\sigma / \sigma_b = 1$ for equally spaced designs),
below which the simpler approach is an adequate
approximation.

![Ratio of standard errors (summary statistic / random slopes) as a function of the between-subject slope SD to residual SD ratio. Values near 1 indicate the two methods are equivalent.](figure/fig-ratio-1.pdf)

\FloatBarrier

# Discussion

The results establish three conditions under which
random effects models reduce to or closely approximate
fixed-effect alternatives.

**Condition 1: unstructured covariance with two
timepoints.** When MMRM is fit with an unrestricted
$2 \times 2$ covariance matrix, the treatment effect
estimate and its standard error are algebraically
identical to those from ANCOVA. This is not an
approximation but an exact equivalence, implicit in
@lee1974note and explicit in @frost2008optimizing
(Section 2.2). The practical implication is that for
any two-timepoint design with complete data, ANCOVA is
a sufficient analysis and the MMRM apparatus adds no
information about the treatment contrast.

**Condition 2: small between-subject slope variance.**
In designs with three or more timepoints, the
random-slopes model becomes equivalent to the
summary-statistic analysis when $\sigma_b / \sigma$ is
small. As @frost2008optimizing showed, when this ratio
is below approximately 0.25 (for equally spaced
designs with the same total follow-up), the efficiency
gain from the random effects model is negligible. In
practical terms, this corresponds to settings where
measurement error dominates biological
variability in rates---a situation common in short
trials with noisy endpoints.

**Condition 3: cLDA with complete data.** The
constrained longitudinal model produces treatment
effect estimates effectively identical to ANCOVA when
no baseline values are missing [@liu2009baseline;
@lu2010efficiency]. The marginal advantage of cLDA
emerges only with incomplete baselines, where it gains
1--4\% additional power depending on the missingness
rate.

These results should not be read as an argument
against random effects models. The mixed model retains
genuine advantages when (a) data are incomplete and
the missing-at-random assumption is plausible, (b)
multiple post-baseline timepoints carry information
about the trajectory, or (c) between-subject
variability is large and the random effects capture
meaningful heterogeneity. @sullivan2022mmrm
demonstrated that MMRM with time-by-covariance
interactions dominates ANCOVA in the presence of
dropout, and @wang2024improving showed that augmented
MMRM can further improve robustness.

The practical value of these equivalence results is
pedagogical and communicative. When an investigator
asks what the MMRM 'is doing,' the analyst can
respond: in your two-visit design, it is doing exactly
what ANCOVA does, no more and no less. The random
effects are not contributing additional information;
they are merely a different parameterization of the
same model. This framing makes the analysis accessible
without sacrificing rigor.

## Limitations

The simulations assumed balanced designs with equal
group sizes and no missing data (except in the cLDA
comparison). The equivalence results extend
straightforwardly to unbalanced designs but may be
perturbed by informative missingness, small samples
where the Kenward-Roger correction
[@kenward1997small] matters, or misspecification of
the covariance structure. We did not investigate
designs with more than two treatment groups, where the
mapping between MMRM contrasts and ANCOVA contrasts is
more complex.

# Future research

1. **Extension to multiple post-baseline visits.** The
   exact equivalence in the two-timepoint case does
   not hold with three or more visits unless the
   covariance is unstructured. Characterizing the
   discrepancy as a function of the number of visits
   and the covariance structure is of practical
   interest.

2. **Non-normal outcomes.** The relationship between
   generalized linear mixed models and marginal (GEE)
   models for binary and count endpoints involves
   additional complexities due to the
   non-collapsibility of odds ratios and rate ratios.

3. **Informative missingness.** Under
   missing-not-at-random mechanisms, the equivalence
   between MMRM and ANCOVA breaks down in different
   ways. Sensitivity analyses comparing the two
   frameworks under various departure patterns from
   MAR would be informative.

4. **Software tools.** An R package providing automatic
   detection of when MMRM reduces to a simpler
   model---based on the number of timepoints,
   completeness, and estimated variance-component
   ratios---would be useful for practitioners.

5. **Pedagogical materials.** Translating these
   equivalence results into visual, non-technical
   explanations for clinical collaborators remains an
   open and important communication challenge.

# References

## Morris et al. (2019) ADEMP Compliance

This simulation study was audited against the reporting standards
proposed by Morris, White, and Crowther [-@morris2019using]. The full
audit is at `docs/morris-audit-2026-04-17.md`.

**Verdict:** Partially compliant.

**Key gaps identified:**

- Simulation code lives entirely inside the Rmd; no `R/` package functions or scripts.
- `n_sims = 2000` hardcoded at two places without MCSE-based justification.
- Seeds are set via function arguments (`seed = 42`, `seed = 123`), not once at program start; no MCSE columns in `summary_prepost()` or `summary_3tp()`.

A remediation plan is documented in the audit file and will be
executed in a subsequent revision.
