# Part 3: Statistical Comparison (Responders vs Non-Responders) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine which immune cell populations differ in relative frequency between miraclib responders and non-responders (melanoma, PBMC), with statistics rigorous enough to convince a skeptical reviewer, plus a boxplot per population.

**Architecture:** A new `analysis/comparison.py` module. It reuses Part 2's `sample_frequencies` (the "summary table") as the single source of relative frequencies, filters to the cohort, and applies a layered statistical strategy: a linear mixed-effects model (primary), a non-parametric per-subject Mann-Whitney U cross-check, Benjamini-Hochberg FDR correction, effect sizes with bootstrap CIs, and a boxplot figure. All functions are pure (return DataFrames / a Plotly figure); the pipeline and dashboard consume them.

**Tech Stack:** Python 3.11+, pandas, scipy (Mann-Whitney, bootstrap), statsmodels (mixed model), plotly (boxplots).

---

## Statistical Methodology & Justifications

This section is the scientific rationale the README/dashboard will cite. Every
test below exists to defend a specific threat to validity.

**The cohort.** Melanoma + miraclib + PBMC, response ∈ {yes, no}: 656 subjects
(331 responders / 325 non-responders), each contributing 3 PBMC samples at
timepoints 0/7/14. Response is a **between-subject** property (constant within a
subject).

**The central threat — pseudoreplication.** The 3 samples per subject are
repeated measures, not independent observations. Naively pooling all 1,968
samples and running a two-sample test treats the effective sample size as ~3×
larger than the number of independent units (subjects), deflating standard
errors and inflating false positives. Every design choice below is about
getting the *unit of analysis* right (the subject, not the sample).

1. **Primary test — linear mixed-effects model**, per population:
   `percentage ~ response + C(time)` with a random intercept per subject
   (`groups = subject_id`).
   - *Why:* A random intercept models the within-subject correlation of repeated
     timepoints, so the standard error on the `response` effect reflects 656
     independent subjects, not 1,968 correlated samples. It uses all the data
     (no discarding timepoints) while remaining valid under repeated measures —
     the standard approach for this design. `C(time)` absorbs any systematic
     drift in frequencies across the treatment course so it doesn't contaminate
     the response contrast.
   - *Report:* the `response` coefficient (mean %-point difference,
     responder − non-responder), its 95% CI, the p-value, and the **random-effect
     variance** (`group_var`). A near-zero `group_var` is informative: it means
     repeated samples behave almost independently, so pseudoreplication happened
     to be mild — but we only know that *because* we fit the model.
   - *Robustness note:* when `group_var` sits on the boundary (≈0), statsmodels
     emits a `ConvergenceWarning`. This is expected and benign; the code
     suppresses that specific warning and records `group_var` instead of hiding
     the diagnostic.

2. **Multiple-comparison control — Benjamini-Hochberg FDR** across the 5
   populations.
   - *Why:* Testing 5 populations at α=0.05 gives a ~23% chance of at least one
     false positive under the null. BH controls the false discovery rate (the
     expected fraction of flagged populations that are false) while keeping more
     power than Bonferroni — the right trade-off for exploratory biomarker
     screening, where the goal is a trustworthy shortlist for follow-up.
   - *Report:* both the raw p and the adjusted q for every population; flag
     q < 0.05.

3. **Distribution-free cross-check — Mann-Whitney U on per-subject means.**
   Collapse each subject to one value per population (mean % across their
   timepoints), then compare responders vs non-responders.
   - *Why:* The mixed model assumes approximately normal, homoscedastic
     residuals; relative frequencies are bounded [0,100] and can be skewed.
     Mann-Whitney U assumes neither normality nor equal variance and is robust to
     outliers, testing whether one group stochastically dominates. Collapsing to
     one value per subject **restores independence**, so the test is valid with
     no pseudoreplication. Agreement between the parametric and non-parametric
     routes is strong evidence; disagreement is a flag to report honestly (a
     population significant under one but not the other is a *candidate*, not a
     confirmed finding).
   - *Report:* U, p, BH-adjusted q, and the **rank-biserial correlation** as a
     non-parametric effect size, plus group medians.

4. **Effect size & uncertainty — bootstrap 95% CI** for the difference in
   per-subject means, per population.
   - *Why:* A p-value says whether an effect is distinguishable from zero, not
     whether it *matters*. For a response-predicting biomarker, a statistically
     significant but <1 %-point difference is clinically useless. Reporting the
     mean difference with a bootstrap CI (non-parametric, no normality
     assumption) makes magnitude and precision explicit, and lets the reader
     judge practical relevance. The bootstrap is seeded for reproducibility.

5. **Visualization — one boxplot per population, responders vs non-responders**,
   built from the per-subject means with individual points overlaid.
   - *Why:* Boxplots expose distribution shape, spread, overlap, and outliers —
     letting the reader see *why* a test did or didn't reach significance rather
     than trusting a number. Using the per-subject unit (matching the tests)
     keeps the visual honest about sample size and avoids visually implying
     1,968 independent points.

6. **Concordance-based significance rule.** A population is reported as a
   **primary hit** if it passes the mixed model (`q_mixed < α`), and as
   **concordant** only if it *also* passes the non-parametric test
   (`q_mw < α`). This two-tier rule prevents overselling a signal that hinges on
   one method's assumptions.

*(Baseline-only prediction — restricting to timepoint 0 — is the natural
extension for a pre-treatment predictive biomarker and is supported by the
`timepoint` parameter of `cohort_frequencies`; wiring it into the pipeline is a
Part 5 concern, not this plan.)*

---

## Global Constraints

- Analysis functions are pure: take a connection or a DataFrame (+ explicit
  params), **return** a DataFrame or a Plotly figure; never print, write files,
  or import Streamlit.
- Cohort filter for the headline analysis: `condition='melanoma'`,
  `treatment='miraclib'`, `sample_type='PBMC'`, `response ∈ {'yes','no'}`
  (exclude NULL/other responses).
- Relative frequencies come from Part 2's `analysis.frequencies.sample_frequencies`
  (single source of truth — "using the data reported in the summary table").
- Unit of analysis is the **subject** for independence-restoring steps
  (Mann-Whitney, bootstrap, boxplot); the mixed model uses all samples with a
  subject random effect.
- Population order is canonical: reuse `load_data.POPULATION_COLUMNS`
  (`b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte`).
- Responder coefficient/difference is oriented **responder − non-responder**
  (reference group `no`).
- Bootstrap is seeded (default `seed=0`) for reproducible CIs.
- The mixed model must suppress only the two `statsmodels` diagnostics that fire
  when the random-effect variance collapses to ≈0 — `ConvergenceWarning` and the
  `"Random effects covariance is singular"` `UserWarning` (both benign here;
  record `group_var` instead). Test output must otherwise be pristine.
- New dependencies pinned in `requirements.txt`: `scipy`, `statsmodels`,
  `plotly`.

---

## File Structure

- `requirements.txt` — add `scipy`, `statsmodels`, `plotly` (modify).
- `analysis/comparison.py` — cohort extraction, statistical primitives,
  orchestration, and the boxplot figure. Built up across Tasks 1–4.
- `tests/test_comparison_cohort.py` — cohort filtering + per-subject collapse (Task 1).
- `tests/test_comparison_stats.py` — the four statistical primitives (Task 2).
- `tests/test_comparison.py` — `compare_responders` orchestration + FDR + flags (Task 3).
- `tests/test_comparison_viz.py` — boxplot figure smoke tests (Task 4).

---

### Task 1: Cohort extraction + per-subject collapse

**Files:**
- Modify: `requirements.txt`
- Create: `analysis/comparison.py`
- Test: `tests/test_comparison_cohort.py`

**Interfaces:**
- Consumes: `analysis.frequencies.sample_frequencies(conn)`; `load_data.POPULATION_COLUMNS`; Part 1 tables.
- Produces:
  - `analysis.comparison.POPULATIONS: list[str]` (canonical order).
  - `analysis.comparison.cohort_frequencies(conn, condition='melanoma', treatment='miraclib', sample_type='PBMC', timepoint=None) -> DataFrame` with columns `sample, total_count, population, count, percentage, subject_id, time, response, sex` filtered to the cohort (response ∈ {yes,no}), optionally to one `timepoint`.
  - `analysis.comparison.subject_means(cohort) -> DataFrame` with columns `subject_id, response, population, percentage` (each subject's mean % per population).

- [ ] **Step 1: Add stats dependencies to `requirements.txt`**

Append these three lines to `requirements.txt` (keep the existing `pandas`/`pytest` lines):

```text
scipy==1.15.3
statsmodels==0.14.6
plotly==6.9.0
```

Then install them: `cd /Users/zoyagoel/teiko && python -m pip install -r requirements.txt`

- [ ] **Step 2: Write the failing test `tests/test_comparison_cohort.py`**

```python
import pandas as pd

from db.connection import get_connection
from load_data import build_database
from analysis.comparison import cohort_frequencies, subject_means, POPULATIONS

# 2 melanoma+miraclib+PBMC subjects (1 responder, 1 non-responder), each 2 timepoints,
# plus rows that must be EXCLUDED: wrong condition, wrong treatment, wrong sample_type,
# and a null-response (healthy) subject.
FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100
prj1,s0,melanoma,50,M,miraclib,yes,smp1,PBMC,7,20,20,30,40,90
prj1,s1,melanoma,60,F,miraclib,no,smp2,PBMC,0,5,20,30,40,105
prj1,s1,melanoma,60,F,miraclib,no,smp3,PBMC,7,15,20,30,40,95
prj1,s2,melanoma,55,M,miraclib,yes,smp4,WB,0,10,20,30,40,100
prj1,s3,carcinoma,55,M,miraclib,yes,smp5,PBMC,0,10,20,30,40,100
prj1,s4,melanoma,55,M,phauximab,no,smp6,PBMC,0,10,20,30,40,100
prj2,s5,healthy,40,F,none,,smp7,PBMC,0,5,5,5,5,5
"""


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_cohort_filters_to_melanoma_miraclib_pbmc_with_response(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    # Only s0 and s1 qualify: 2 subjects x 2 timepoints x 5 populations = 20 rows
    assert set(cohort["subject_id"]) == {"s0", "s1"}
    assert len(cohort) == 20
    assert set(cohort["response"]) == {"yes", "no"}
    # Excluded units must not appear
    assert "s2" not in set(cohort["subject_id"])  # WB
    assert "s3" not in set(cohort["subject_id"])  # carcinoma
    assert "s4" not in set(cohort["subject_id"])  # phauximab
    assert "s5" not in set(cohort["subject_id"])  # healthy / null response


def test_cohort_has_expected_columns(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    for col in ["sample", "total_count", "population", "count", "percentage",
                "subject_id", "time", "response", "sex"]:
        assert col in cohort.columns


def test_timepoint_filter(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path), timepoint=0)
    assert set(cohort["time"]) == {0}
    assert len(cohort) == 10  # 2 subjects x 1 timepoint x 5 populations


def test_subject_means_collapses_to_one_row_per_subject_population(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    sm = subject_means(cohort)
    assert list(sm.columns) == ["subject_id", "response", "population", "percentage"]
    assert len(sm) == 10  # 2 subjects x 5 populations
    # s0 b_cell: sample smp0 total=200 ->5%, smp1 total=200 ->10%; mean=7.5
    s0_b = sm[(sm.subject_id == "s0") & (sm.population == "b_cell")].iloc[0]
    assert abs(s0_b["percentage"] - 7.5) < 1e-9
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `cd /Users/zoyagoel/teiko && python -m pytest tests/test_comparison_cohort.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'analysis.comparison'`.

- [ ] **Step 4: Create `analysis/comparison.py` with cohort functions**

```python
"""Part 3: compare cell-population frequencies between miraclib responders and
non-responders (melanoma, PBMC). Builds on Part 2's frequency summary table."""
import pandas as pd

from load_data import POPULATION_COLUMNS
from analysis.frequencies import sample_frequencies

POPULATIONS = list(POPULATION_COLUMNS)

_METADATA_QUERY = """
SELECT sa.sample_id                         AS sample,
       sa.subject_id                        AS subject_id,
       sa.sample_type                       AS sample_type,
       sa.time_from_treatment_start         AS time,
       su.condition                         AS condition,
       su.treatment                         AS treatment,
       su.response                          AS response,
       su.sex                               AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
"""


def cohort_frequencies(conn, condition="melanoma", treatment="miraclib",
                       sample_type="PBMC", timepoint=None):
    """Per-sample relative frequencies for the response-comparison cohort.

    Reuses Part 2's sample_frequencies (single source of truth) and joins the
    sample/subject metadata, then filters to the cohort and to responders /
    non-responders only (response in {'yes','no'}).
    """
    freq = sample_frequencies(conn)
    meta = pd.read_sql_query(_METADATA_QUERY, conn)
    df = freq.merge(meta, on="sample", how="inner")
    mask = (
        (df["condition"] == condition)
        & (df["treatment"] == treatment)
        & (df["sample_type"] == sample_type)
        & (df["response"].isin(["yes", "no"]))
    )
    df = df[mask]
    if timepoint is not None:
        df = df[df["time"] == timepoint]
    keep = ["sample", "total_count", "population", "count", "percentage",
            "subject_id", "time", "response", "sex"]
    return df[keep].reset_index(drop=True)


def subject_means(cohort):
    """Collapse to one value per (subject, population): mean % across timepoints.

    Restores independence (one row per subject) for the non-parametric and
    bootstrap steps.
    """
    return (
        cohort.groupby(["subject_id", "response", "population"])["percentage"]
        .mean()
        .reset_index()
    )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_comparison_cohort.py -v`
Expected: PASS (4 tests), pristine.

- [ ] **Step 6: Commit**

```bash
git add requirements.txt analysis/comparison.py tests/test_comparison_cohort.py
git commit -m "feat(part3): cohort extraction and per-subject collapse"
```

---

### Task 2: Statistical primitives

**Files:**
- Modify: `analysis/comparison.py` (append functions)
- Test: `tests/test_comparison_stats.py`

**Interfaces:**
- Consumes: `POPULATIONS`, a `cohort` DataFrame (columns incl. `subject_id, response, time, population, percentage`), and `subject_means`.
- Produces (all per-population, one row each, `population` column):
  - `benjamini_hochberg(pvalues) -> numpy.ndarray` of BH-adjusted q-values.
  - `mannwhitney_by_population(cohort, populations=POPULATIONS) -> DataFrame` cols: `population, mw_U, mw_p, rank_biserial, n_responders, n_nonresponders, median_resp, median_nonresp`.
  - `mixedlm_by_population(cohort, populations=POPULATIONS) -> DataFrame` cols: `population, coef, ci_low, ci_high, mixed_p, group_var`.
  - `bootstrap_diff_ci(cohort, populations=POPULATIONS, n_boot=2000, seed=0) -> DataFrame` cols: `population, mean_diff, boot_ci_low, boot_ci_high`.

- [ ] **Step 1: Write the failing test `tests/test_comparison_stats.py`**

```python
import numpy as np
import pandas as pd

from analysis.comparison import (
    benjamini_hochberg,
    mannwhitney_by_population,
    mixedlm_by_population,
    bootstrap_diff_ci,
)


def _synthetic_cohort(effect=5.0, seed=0, n_per=40):
    """Two populations: 'A' has a planted responder effect, 'B' has none.
    30/40 subjects per group, 3 timepoints each, subject-level random intercept."""
    rng = np.random.default_rng(seed)
    rows = []
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(n_per):
            sid = f"{resp}_{s}"
            intercept = rng.normal(0, 1)  # subject random effect
            for t in (0, 7, 14):
                rows.append(dict(subject_id=sid, response=resp, time=t,
                                 population="A",
                                 percentage=20 + effect * is_r + intercept + rng.normal(0, 1)))
                rows.append(dict(subject_id=sid, response=resp, time=t,
                                 population="B",
                                 percentage=20 + intercept + rng.normal(0, 1)))
    return pd.DataFrame(rows)


def test_benjamini_hochberg_known_values():
    q = benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
    # BH: q_i = p_i * n / rank, then cumulative-min from the top.
    # All equal to 0.05 here since p_i * 5 / i == 0.05 for every i.
    assert np.allclose(q, [0.05, 0.05, 0.05, 0.05, 0.05])
    # q-values are clipped to [0,1] and monotone non-decreasing in p
    q2 = benjamini_hochberg([0.9, 0.001])
    assert np.all(q2 <= 1.0)
    assert q2[1] <= q2[0]


def test_mannwhitney_detects_planted_effect():
    res = mannwhitney_by_population(_synthetic_cohort(effect=5.0), ["A", "B"])
    a = res[res.population == "A"].iloc[0]
    b = res[res.population == "B"].iloc[0]
    assert a["mw_p"] < 0.05          # A differs
    assert a["rank_biserial"] > 0    # responders higher
    assert b["mw_p"] > 0.05          # B does not
    assert a["n_responders"] == 40 and a["n_nonresponders"] == 40


def test_mixedlm_recovers_effect_and_reports_group_var():
    res = mixedlm_by_population(_synthetic_cohort(effect=5.0), ["A", "B"])
    a = res[res.population == "A"].iloc[0]
    assert abs(a["coef"] - 5.0) < 1.5     # coefficient near planted effect
    assert a["mixed_p"] < 0.05
    assert a["ci_low"] < a["coef"] < a["ci_high"]
    assert "group_var" in res.columns and a["group_var"] >= 0


def test_bootstrap_ci_brackets_diff_and_is_reproducible():
    coh = _synthetic_cohort(effect=5.0)
    r1 = bootstrap_diff_ci(coh, ["A", "B"], n_boot=500, seed=0)
    r2 = bootstrap_diff_ci(coh, ["A", "B"], n_boot=500, seed=0)
    a = r1[r1.population == "A"].iloc[0]
    assert a["boot_ci_low"] < a["mean_diff"] < a["boot_ci_high"]
    assert a["boot_ci_low"] > 0            # whole CI above 0 for a strong effect
    pd.testing.assert_frame_equal(r1, r2)  # seeded => identical
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_comparison_stats.py -v`
Expected: FAIL — ImportError (`benjamini_hochberg` etc. not defined).

- [ ] **Step 3: Append the statistical primitives to `analysis/comparison.py`**

```python
import warnings

import numpy as np
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR-adjusted q-values (order-preserving, clipped [0,1])."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downward
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def mannwhitney_by_population(cohort, populations=POPULATIONS):
    """Two-sided Mann-Whitney U on per-subject means, responder vs non-responder."""
    sm = subject_means(cohort)
    rows = []
    for pop in populations:
        d = sm[sm["population"] == pop]
        r = d[d["response"] == "yes"]["percentage"].to_numpy()
        nr = d[d["response"] == "no"]["percentage"].to_numpy()
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        n1, n2 = len(r), len(nr)
        # rank-biserial correlation: +1 when responders stochastically dominate
        rank_biserial = 2.0 * u / (n1 * n2) - 1.0
        rows.append(dict(population=pop, mw_U=float(u), mw_p=float(p),
                         rank_biserial=rank_biserial,
                         n_responders=n1, n_nonresponders=n2,
                         median_resp=float(np.median(r)),
                         median_nonresp=float(np.median(nr))))
    return pd.DataFrame(rows)


def mixedlm_by_population(cohort, populations=POPULATIONS):
    """Linear mixed model per population: percentage ~ response (+ time),
    random intercept per subject. Suppresses the two benign diagnostics that
    fire when the random-effect variance collapses to ~0 (ConvergenceWarning and
    the "Random effects covariance is singular" UserWarning) and records the
    random-effect variance (group_var)."""
    rows = []
    for pop in populations:
        d = cohort[cohort["population"] == pop]
        formula = "percentage ~ response"
        if d["time"].nunique() > 1:
            formula += " + C(time)"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.filterwarnings(
                "ignore", message="Random effects covariance is singular")
            model = smf.mixedlm(formula, d, groups=d["subject_id"])
            fit = model.fit(method="lbfgs")
        term = "response[T.yes]"
        ci = fit.conf_int()
        rows.append(dict(population=pop,
                         coef=float(fit.params[term]),
                         ci_low=float(ci.loc[term, 0]),
                         ci_high=float(ci.loc[term, 1]),
                         mixed_p=float(fit.pvalues[term]),
                         group_var=float(fit.cov_re.iloc[0, 0])))
    return pd.DataFrame(rows)


def bootstrap_diff_ci(cohort, populations=POPULATIONS, n_boot=2000, seed=0):
    """Bootstrap 95% CI for the responder - non-responder difference in
    per-subject mean %, per population. Seeded for reproducibility."""
    rng = np.random.default_rng(seed)
    sm = subject_means(cohort)
    rows = []
    for pop in populations:
        d = sm[sm["population"] == pop]
        r = d[d["response"] == "yes"]["percentage"].to_numpy()
        nr = d[d["response"] == "no"]["percentage"].to_numpy()
        mean_diff = r.mean() - nr.mean()
        boots = np.empty(n_boot)
        for b in range(n_boot):
            rb = rng.choice(r, len(r), replace=True).mean()
            nrb = rng.choice(nr, len(nr), replace=True).mean()
            boots[b] = rb - nrb
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append(dict(population=pop, mean_diff=float(mean_diff),
                         boot_ci_low=float(lo), boot_ci_high=float(hi)))
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_comparison_stats.py -v`
Expected: PASS (4 tests), pristine (no ConvergenceWarning leaking to output).

- [ ] **Step 5: Commit**

```bash
git add analysis/comparison.py tests/test_comparison_stats.py
git commit -m "feat(part3): mixed model, Mann-Whitney, FDR, bootstrap primitives"
```

---

### Task 3: `compare_responders` orchestration

**Files:**
- Modify: `analysis/comparison.py` (append)
- Test: `tests/test_comparison.py`

**Interfaces:**
- Consumes: `cohort_frequencies`, `mixedlm_by_population`, `mannwhitney_by_population`, `bootstrap_diff_ci`, `benjamini_hochberg`.
- Produces: `analysis.comparison.compare_responders(conn, condition='melanoma', treatment='miraclib', sample_type='PBMC', alpha=0.05, seed=0) -> DataFrame`, one row per population (canonical order), with columns: `population, n_responders, n_nonresponders, mean_resp, mean_nonresp, mean_diff, boot_ci_low, boot_ci_high, coef, ci_low, ci_high, mixed_p, q_mixed, group_var, mw_U, mw_p, q_mw, rank_biserial, median_resp, median_nonresp, significant_primary, concordant`.
  - `q_mixed` = BH over `mixed_p`; `q_mw` = BH over `mw_p`.
  - `significant_primary` = `q_mixed < alpha`; `concordant` = `(q_mixed < alpha) & (q_mw < alpha)`.

- [ ] **Step 1: Write the failing test `tests/test_comparison.py`**

```python
import numpy as np
import pandas as pd

from db.connection import get_connection
from load_data import build_database
from analysis.comparison import compare_responders, POPULATIONS


def _fixture_csv(seed=0):
    """~24 melanoma+miraclib+PBMC subjects (12 responders/12 non), 3 timepoints,
    with a planted cd4_t_cell elevation in responders; other populations flat."""
    rng = np.random.default_rng(seed)
    header = ("project,subject,condition,age,sex,treatment,response,sample,"
              "sample_type,time_from_treatment_start,b_cell,cd8_t_cell,"
              "cd4_t_cell,nk_cell,monocyte")
    lines = [header]
    sample_no = 0
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(12):
            sid = f"{resp}{s:02d}"
            for t in (0, 7, 14):
                cd4 = 300 + 40 * is_r + rng.integers(-10, 10)
                b = 100 + rng.integers(-10, 10)
                cd8 = 250 + rng.integers(-10, 10)
                nk = 150 + rng.integers(-10, 10)
                mono = 200 + rng.integers(-10, 10)
                lines.append(
                    f"prj1,{sid},melanoma,50,M,miraclib,{resp},smp{sample_no:04d},"
                    f"PBMC,{t},{b},{cd8},{cd4},{nk},{mono}")
                sample_no += 1
    return "\n".join(lines) + "\n"


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(_fixture_csv())
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_result_has_one_row_per_population_in_order(tmp_path):
    res = compare_responders(_conn(tmp_path))
    assert list(res["population"]) == POPULATIONS
    assert len(res) == 5


def test_result_has_all_expected_columns(tmp_path):
    res = compare_responders(_conn(tmp_path))
    for col in ["population", "n_responders", "n_nonresponders", "mean_resp",
                "mean_nonresp", "mean_diff", "boot_ci_low", "boot_ci_high",
                "coef", "ci_low", "ci_high", "mixed_p", "q_mixed", "group_var",
                "mw_U", "mw_p", "q_mw", "rank_biserial", "median_resp",
                "median_nonresp", "significant_primary", "concordant"]:
        assert col in res.columns


def test_fdr_and_flags_are_well_formed(tmp_path):
    res = compare_responders(_conn(tmp_path))
    # q-values are valid probabilities and never smaller than their raw p
    assert (res["q_mixed"] >= res["mixed_p"] - 1e-9).all()
    assert (res["q_mixed"] <= 1.0).all() and (res["q_mixed"] >= 0.0).all()
    # flags are boolean and concordant implies primary
    assert res["significant_primary"].dtype == bool
    assert res["concordant"].dtype == bool
    assert (~res["concordant"] | res["significant_primary"]).all()
    assert res["n_responders"].iloc[0] == 12
    assert res["n_nonresponders"].iloc[0] == 12


def test_planted_cd4_effect_is_detected(tmp_path):
    res = compare_responders(_conn(tmp_path)).set_index("population")
    # cd4 was planted higher in responders
    assert res.loc["cd4_t_cell", "mean_diff"] > 0
    assert res.loc["cd4_t_cell", "significant_primary"]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_comparison.py -v`
Expected: FAIL — `ImportError: cannot import name 'compare_responders'`.

- [ ] **Step 3: Append `compare_responders` to `analysis/comparison.py`**

```python
def compare_responders(conn, condition="melanoma", treatment="miraclib",
                       sample_type="PBMC", alpha=0.05, seed=0):
    """Full responder vs non-responder comparison, one row per population.

    Combines the mixed-effects model (primary), Mann-Whitney U on per-subject
    means (non-parametric cross-check), bootstrap effect-size CIs, and
    Benjamini-Hochberg FDR across the five populations.
    """
    cohort = cohort_frequencies(conn, condition, treatment, sample_type)

    mm = mixedlm_by_population(cohort, POPULATIONS)
    mw = mannwhitney_by_population(cohort, POPULATIONS)
    boot = bootstrap_diff_ci(cohort, POPULATIONS, seed=seed)

    sm = subject_means(cohort)
    group_means = (
        sm.pivot_table(index="population", columns="response",
                       values="percentage", aggfunc="mean")
        .rename(columns={"yes": "mean_resp", "no": "mean_nonresp"})
        .reset_index()
    )

    res = (mm.merge(mw, on="population")
              .merge(boot, on="population")
              .merge(group_means, on="population"))

    res["q_mixed"] = benjamini_hochberg(res["mixed_p"].to_numpy())
    res["q_mw"] = benjamini_hochberg(res["mw_p"].to_numpy())
    res["significant_primary"] = res["q_mixed"] < alpha
    res["concordant"] = (res["q_mixed"] < alpha) & (res["q_mw"] < alpha)

    # canonical population order + tidy column order
    res["population"] = pd.Categorical(res["population"], categories=POPULATIONS,
                                       ordered=True)
    res = res.sort_values("population").reset_index(drop=True)
    res["population"] = res["population"].astype(str)
    columns = ["population", "n_responders", "n_nonresponders", "mean_resp",
               "mean_nonresp", "mean_diff", "boot_ci_low", "boot_ci_high",
               "coef", "ci_low", "ci_high", "mixed_p", "q_mixed", "group_var",
               "mw_U", "mw_p", "q_mw", "rank_biserial", "median_resp",
               "median_nonresp", "significant_primary", "concordant"]
    return res[columns]
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_comparison.py -v`
Expected: PASS (4 tests), pristine.

- [ ] **Step 5: Run the full suite**

Run: `python -m pytest -v`
Expected: PASS (all Parts 1–3 tests).

- [ ] **Step 6: Verify against the real cohort and record the finding**

Ensure the DB exists (`python load_data.py`), then run:

```bash
python -c "
import pandas as pd
from db.connection import get_connection
from analysis.comparison import compare_responders
pd.set_option('display.width', 240)
res = compare_responders(get_connection())
cols = ['population','mean_resp','mean_nonresp','mean_diff','coef','mixed_p','q_mixed','mw_p','q_mw','group_var','significant_primary','concordant']
print(res[cols].round(4).to_string(index=False))
print()
print('primary hits (q_mixed<0.05):', list(res.loc[res.significant_primary,'population']) or 'NONE')
print('concordant (both tests):    ', list(res.loc[res.concordant,'population']) or 'NONE')
"
```

Expected on the real data: `cd4_t_cell` is the only **primary hit** (`q_mixed`≈0.02–0.03, `mean_diff`≈+0.6), it is **not concordant** (`q_mw`≈0.06), all other populations non-significant, and `group_var` is ≈0 for every population. If the qualitative result differs (e.g. a different or additional primary hit, or cd4 not flagged), stop and investigate before committing.

- [ ] **Step 7: Commit**

```bash
git add analysis/comparison.py tests/test_comparison.py
git commit -m "feat(part3): compare_responders orchestration with FDR and flags"
```

---

### Task 4: Boxplot visualization

**Files:**
- Modify: `analysis/comparison.py` (append)
- Test: `tests/test_comparison_viz.py`

**Interfaces:**
- Consumes: `subject_means`, `POPULATIONS`, a `cohort` DataFrame.
- Produces: `analysis.comparison.responder_boxplots(cohort, populations=POPULATIONS) -> plotly.graph_objects.Figure` — one boxplot panel per population, responders vs non-responders, built from per-subject means with points overlaid.

- [ ] **Step 1: Write the failing test `tests/test_comparison_viz.py`**

```python
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from analysis.comparison import responder_boxplots


def _cohort():
    rng = np.random.default_rng(0)
    rows = []
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(10):
            for t in (0, 7, 14):
                for pop in ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]:
                    rows.append(dict(subject_id=f"{resp}{s}", response=resp, time=t,
                                     population=pop,
                                     percentage=20 + (5 * is_r if pop == "cd4_t_cell" else 0)
                                     + rng.normal(0, 1)))
    return pd.DataFrame(rows)


def test_returns_plotly_figure_with_data():
    fig = responder_boxplots(_cohort())
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0


def test_covers_all_five_populations():
    fig = responder_boxplots(_cohort())
    text = str(fig.to_dict())
    for pop in ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]:
        assert pop in text


def test_both_response_groups_present():
    fig = responder_boxplots(_cohort())
    text = str(fig.to_dict())
    assert "yes" in text and "no" in text
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_comparison_viz.py -v`
Expected: FAIL — `ImportError: cannot import name 'responder_boxplots'`.

- [ ] **Step 3: Append `responder_boxplots` to `analysis/comparison.py`**

```python
import plotly.express as px


def responder_boxplots(cohort, populations=POPULATIONS):
    """One boxplot panel per population comparing responders vs non-responders,
    built from per-subject mean % (individual points overlaid)."""
    sm = subject_means(cohort)
    fig = px.box(
        sm, x="response", y="percentage", color="response",
        facet_col="population", facet_col_wrap=len(populations),
        points="all", boxmode="overlay",
        category_orders={"population": list(populations),
                         "response": ["no", "yes"]},
        labels={"percentage": "Relative frequency (%)", "response": "Response"},
    )
    # strip the "population=" prefix plotly adds to each facet title
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(
        showlegend=False,
        title="Cell-population relative frequency: responders vs non-responders",
    )
    return fig
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_comparison_viz.py -v`
Expected: PASS (3 tests), pristine.

- [ ] **Step 5: Run the full suite**

Run: `python -m pytest -v`
Expected: PASS (all Parts 1–3).

- [ ] **Step 6: Commit**

```bash
git add analysis/comparison.py tests/test_comparison_viz.py
git commit -m "feat(part3): responder vs non-responder boxplots"
```

---

## Self-Review

**Spec coverage (against instructions Part 3 and overarching design §5):**
- Compare frequencies of melanoma/miraclib/PBMC responders vs non-responders → `cohort_frequencies` (Task 1) + `compare_responders` (Task 3). ✓
- Use only PBMC samples, response yes/no → cohort filter (Task 1), tested by exclusion cases. ✓
- "Using the data reported in the summary table" → `cohort_frequencies` reuses `sample_frequencies` (Task 1). ✓
- Boxplot per immune cell population, responders vs non-responders → `responder_boxplots` (Task 4). ✓
- Report which populations significantly differ, with supporting statistics → mixed model + MWU + FDR + effect sizes + flags in `compare_responders` (Tasks 2–3); real-data finding recorded (Task 3 Step 6). ✓
- Robust justification for each test → Methodology section above. ✓

**Placeholder scan:** No TBD/TODO; every code and test step is complete. ✓

**Type consistency:** `POPULATIONS`, `cohort_frequencies`, `subject_means`,
`benjamini_hochberg`, `mannwhitney_by_population`, `mixedlm_by_population`,
`bootstrap_diff_ci`, `compare_responders`, `responder_boxplots` are named
identically across interface blocks, code, and tests. Column names produced by
each primitive match those merged in `compare_responders` and asserted in the
tests (`mixed_p`, `mw_p`, `coef`, `group_var`, `mean_diff`, `boot_ci_low/high`,
`rank_biserial`, `q_mixed`, `q_mw`, `significant_primary`, `concordant`). ✓

**Determinism / test-hygiene notes:** the mixed model suppresses only
`ConvergenceWarning`; synthetic fixtures and the bootstrap are seeded; the
real-data scientific result is *verified and reported* (Task 3 Step 6) but not
asserted in CI, so tests don't become brittle against the dataset.
