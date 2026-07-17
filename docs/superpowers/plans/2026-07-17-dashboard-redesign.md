# Dashboard Redesign (Report/Paper Style) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the dashboard from an interactive explorer to a report/paper style — three tabs (one per instruction part), Part 3 a fixed statistical "paper" with results front and center, and Part 4 surfacing degenerate subsets.

**Architecture:** Keep the thin render layer, `data.py` loaders, theme, and self-building DB. Add two pure display helpers to `components.py`. Rewrite the three section modules as fixed reports and simplify `app.py` (minimal sidebar, part-named tabs).

**Tech Stack:** Python 3.11+, Streamlit, pandas, plotly.

## Global Constraints

- Render layer only; no changes to `analysis/`, the pipeline, or the DB.
- Part 3 is FIXED on condition=melanoma, treatment=miraclib, sample_type=PBMC — no filter selectboxes.
- All methods AND results live in the tabs; sidebar is title + one line only.
- Part 4 surfaces degenerate subsets (healthy conditions and `none` treatment → 0 responder-labelled samples).
- Section files renamed to their parts: `part2.py`, `part3.py`, `part4.py`; old `overview.py`/`response.py`/`subset.py` removed.
- Display helpers are pure and unit-tested.

## File Structure

- `dashboard/components.py` — add `format_stats_table`, `style_stats_table`, `subset_coverage`; remove `style_comparison` (unused after redesign).
- `dashboard/data.py` — add cached `load_subset_coverage`.
- `dashboard/sections/part2.py`, `part3.py`, `part4.py` — the three report tabs (replace overview/response/subset).
- `dashboard/app.py` — minimal sidebar, part-named tabs.
- `tests/test_dashboard.py` — add unit tests for the two helpers; update the AppTest smoke.

---

### Task 1: Display helpers + coverage loader

**Files:**
- Modify: `dashboard/components.py` (add helpers), `dashboard/data.py` (add loader)
- Test: `tests/test_dashboard.py` (add helper unit tests)

**Interfaces:**
- Produces:
  - `components.format_stats_table(comparison) -> DataFrame` — paper-quality display table.
  - `components.style_stats_table(display_df) -> Styler` — highlights significant rows.
  - `components.subset_coverage(metadata) -> DataFrame` cols `condition, treatment, samples, responder_labelled, degenerate`.
  - `data.load_subset_coverage()` (cached) = `subset_coverage(load_sample_metadata())`.

- [ ] **Step 1: Add the helpers to `dashboard/components.py`**

Add these functions (keep `PALETTE`, `configure_plotly`, `kpi_row`, `finding_summary`; you MAY leave `style_comparison` for now — it is removed in Task 3). Add `import pandas as pd` at the top if not present.

```python
def format_stats_table(comparison):
    """Paper-quality per-population display table from a compare_responders result."""
    df = comparison
    return pd.DataFrame({
        "Population": df["population"].values,
        "Δ resp−non (pp)": df["mean_diff"].round(2).values,
        "95% CI (pp)": [f"[{lo:.2f}, {hi:.2f}]"
                        for lo, hi in zip(df["boot_ci_low"], df["boot_ci_high"])],
        "Mixed p": df["mixed_p"].round(4).values,
        "FDR q": df["q_mixed"].round(4).values,
        "MWU q": df["q_mw"].round(4).values,
        "Effect size": df["rank_biserial"].round(3).values,
        "Significant": df["significant_primary"].map({True: "yes", False: "no"}).values,
        "Concordant": df["concordant"].map({True: "yes", False: "no"}).values,
    })


def style_stats_table(display_df):
    """Highlight rows that are significant (Significant == 'yes')."""
    def _highlight(row):
        color = "background-color: #FFF3CD" if row.get("Significant") == "yes" else ""
        return [color] * len(row)
    return display_df.style.apply(_highlight, axis=1)


def subset_coverage(metadata):
    """Baseline (t=0) PBMC subset sizes per condition x treatment, flagging the
    degenerate subsets (no responder/non-responder label)."""
    base = metadata[(metadata["sample_type"] == "PBMC") & (metadata["time"] == 0)]
    rows = []
    for (condition, treatment), grp in base.groupby(["condition", "treatment"]):
        labelled = int(grp["response"].isin(["yes", "no"]).sum())
        rows.append({
            "condition": condition,
            "treatment": treatment,
            "samples": len(grp),
            "responder_labelled": labelled,
            "degenerate": labelled == 0,
        })
    return (pd.DataFrame(rows)
            .sort_values(["condition", "treatment"])
            .reset_index(drop=True))
```

- [ ] **Step 2: Add the cached loader to `dashboard/data.py`**

Add the import and loader (place the import with the other component imports if any, else near the top):

```python
from dashboard.components import subset_coverage


@st.cache_data
def load_subset_coverage():
    return subset_coverage(load_sample_metadata())
```

- [ ] **Step 3: Write the failing helper tests (append to `tests/test_dashboard.py`)**

```python
def _comparison_like():
    import pandas as pd
    return pd.DataFrame([
        {"population": "cd4_t_cell", "mean_diff": 0.6355, "boot_ci_low": 0.21,
         "boot_ci_high": 1.07, "mixed_p": 0.005, "q_mixed": 0.0248, "q_mw": 0.0621,
         "rank_biserial": 0.08, "significant_primary": True, "concordant": False},
        {"population": "b_cell", "mean_diff": -0.199, "boot_ci_low": -0.5,
         "boot_ci_high": 0.1, "mixed_p": 0.17, "q_mixed": 0.33, "q_mw": 0.43,
         "rank_biserial": -0.03, "significant_primary": False, "concordant": False},
    ])


def test_format_stats_table_columns_and_values():
    from dashboard.components import format_stats_table
    out = format_stats_table(_comparison_like())
    assert list(out.columns) == [
        "Population", "Δ resp−non (pp)", "95% CI (pp)", "Mixed p", "FDR q",
        "MWU q", "Effect size", "Significant", "Concordant"]
    cd4 = out[out["Population"] == "cd4_t_cell"].iloc[0]
    assert cd4["Δ resp−non (pp)"] == 0.64
    assert cd4["95% CI (pp)"] == "[0.21, 1.07]"
    assert cd4["Significant"] == "yes"
    assert cd4["Concordant"] == "no"


def test_subset_coverage_flags_degenerate():
    import pandas as pd
    from dashboard.components import subset_coverage
    meta = pd.DataFrame([
        dict(sample="a", subject_id="s1", project="prj1", condition="melanoma",
             treatment="miraclib", sample_type="PBMC", time=0, response="yes", sex="M"),
        dict(sample="b", subject_id="s2", project="prj1", condition="melanoma",
             treatment="miraclib", sample_type="PBMC", time=0, response="no", sex="F"),
        dict(sample="c", subject_id="s3", project="prj2", condition="healthy",
             treatment="none", sample_type="PBMC", time=0, response=None, sex="M"),
        # non-baseline row (t=7) must be excluded from coverage
        dict(sample="d", subject_id="s1", project="prj1", condition="melanoma",
             treatment="miraclib", sample_type="PBMC", time=7, response="yes", sex="M"),
    ])
    cov = subset_coverage(meta).set_index(["condition", "treatment"])
    assert bool(cov.loc[("melanoma", "miraclib"), "degenerate"]) is False
    assert int(cov.loc[("melanoma", "miraclib"), "responder_labelled"]) == 2
    assert int(cov.loc[("melanoma", "miraclib"), "samples"]) == 2  # t=7 excluded
    assert bool(cov.loc[("healthy", "none"), "degenerate"]) is True
    assert int(cov.loc[("healthy", "none"), "responder_labelled"]) == 0
```

- [ ] **Step 4: Run the tests**

Run: `cd /Users/zoyagoel/teiko && python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS (existing tests + 2 new), pristine.

- [ ] **Step 5: Full suite**

Run: `python3 -m pytest -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add dashboard/components.py dashboard/data.py tests/test_dashboard.py
git commit -m "feat(dashboard): stats-table + subset-coverage display helpers"
```

---

### Task 2: App shell + Part 2 + Part 4 tabs (report style)

**Files:**
- Create: `dashboard/sections/part2.py`, `dashboard/sections/part4.py`, `dashboard/sections/part3.py` (stub)
- Modify: `dashboard/app.py`
- Delete: `dashboard/sections/overview.py`, `dashboard/sections/response.py`, `dashboard/sections/subset.py`
- Modify: `tests/test_dashboard.py` (update the AppTest smoke)

**Interfaces:**
- Produces: `part2.render()`, `part4.render()`, `part3.render()` (stub), a rewritten `app.py`.

- [ ] **Step 1: Create `dashboard/sections/part2.py`**

```python
"""Part 2 · Cell-Type Frequencies: the per-sample relative-frequency summary."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_dataset_summary, load_frequencies
from dashboard.components import kpi_row


def render():
    st.subheader("Part 2 · Relative cell-population frequencies")
    st.write("For each sample, the relative frequency of each immune cell "
             "population as a percentage of that sample's total cell count.")

    summary = load_dataset_summary()
    kpi_row([("Samples", summary["samples"]),
             ("Subjects", summary["subjects"]),
             ("Projects", summary["projects"])])

    freq = load_frequencies()
    mean_comp = (freq.groupby("population", as_index=False)["percentage"].mean()
                 .sort_values("percentage", ascending=False))
    fig = px.bar(mean_comp, x="population", y="percentage",
                 title="Mean relative frequency by population (all samples)",
                 labels={"percentage": "Mean %", "population": "Population"})
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("**Per-sample frequency table**")
    st.dataframe(freq, use_container_width=True, hide_index=True)
```

- [ ] **Step 2: Create `dashboard/sections/part4.py`**

```python
"""Part 4 · Subset Analysis: baseline breakdowns, the posed answer, and a
degenerate-subset coverage callout."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_baseline, load_final_answer, load_subset_coverage
from dashboard.components import kpi_row
from analysis.subsets import (
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)


def _bar(df, x, y, title):
    st.plotly_chart(px.bar(df, x=x, y=y, title=title), use_container_width=True)


def render():
    st.subheader("Part 4 · Baseline subset analysis")
    st.write("Melanoma PBMC samples at baseline (time = 0) from miraclib-treated "
             "patients, broken down by project, response, and sex.")

    subset = load_baseline("melanoma", "miraclib", "PBMC")
    kpi_row([("Baseline samples", len(subset)),
             ("Subjects", subset["subject_id"].nunique())])

    c1, c2, c3 = st.columns(3)
    with c1:
        _bar(samples_per_project(subset), "project", "n_samples",
             "Samples per project")
    with c2:
        _bar(subjects_by_response(subset), "response", "n_subjects",
             "Subjects by response")
    with c3:
        _bar(subjects_by_sex(subset), "sex", "n_subjects", "Subjects by sex")

    st.divider()
    st.metric(
        "Avg B cells — melanoma male responders at baseline "
        "(all sample & treatment types)",
        f"{load_final_answer():.2f}",
    )

    st.divider()
    st.markdown("#### Data coverage — which subsets are degenerate")
    st.write("Baseline (time = 0) PBMC subset sizes by condition and treatment. "
             "A subset is **degenerate** when it has no responder/non-responder "
             "label — healthy patients and untreated (`none`) arms carry no "
             "response, so they cannot be compared.")
    st.dataframe(load_subset_coverage(), use_container_width=True, hide_index=True)
```

- [ ] **Step 3: Create `dashboard/sections/part3.py` (stub)**

```python
"""Part 3 · Statistical Analysis. Full report implemented in the next task."""
import streamlit as st


def render():
    st.subheader("Part 3 · Statistical analysis")
    st.info("Statistical analysis report — implemented next.")
```

- [ ] **Step 4: Rewrite `dashboard/app.py`**

```python
"""Loblaw Bio immune-cell analysis dashboard."""
import os
import sys

# streamlit run dashboard/app.py puts dashboard/ on sys.path, not the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st

st.set_page_config(page_title="Loblaw Bio — Immune Cell Analysis",
                   layout="wide")

from dashboard.components import configure_plotly
from dashboard.sections import part2, part3, part4

configure_plotly()

st.sidebar.title("Loblaw Bio")
st.sidebar.caption("Immune cell population analysis — miraclib trial "
                   "(source: cell-count.csv)")

st.title("Immune Cell Population Analysis")

tab2, tab3, tab4 = st.tabs([
    "Part 2 · Cell-Type Frequencies",
    "Part 3 · Statistical Analysis",
    "Part 4 · Subset Analysis",
])
with tab2:
    part2.render()
with tab3:
    part3.render()
with tab4:
    part4.render()
```

- [ ] **Step 5: Delete the old section files**

```bash
rm dashboard/sections/overview.py dashboard/sections/response.py dashboard/sections/subset.py
```

- [ ] **Step 6: Update the AppTest smoke test in `tests/test_dashboard.py`**

Replace the existing `test_app_runs_and_renders_response_finding` (or `test_app_runs_without_exception`) with:

```python
def test_app_runs_without_exception():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
```

- [ ] **Step 7: Run the tests**

Run: `python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS. The app renders Part 2 (table), the Part 3 stub, and Part 4 (breakdowns + coverage) with no exception.

- [ ] **Step 8: Full suite**

Run: `python3 -m pytest -v`
Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add dashboard/app.py dashboard/sections tests/test_dashboard.py
git commit -m "feat(dashboard): report-style app shell, Part 2 and Part 4 tabs"
```

---

### Task 3: Part 3 statistical-analysis paper

**Files:**
- Modify: `dashboard/sections/part3.py` (replace stub with the full paper)
- Modify: `dashboard/components.py` (remove now-unused `style_comparison`)
- Modify: `tests/test_dashboard.py` (strengthen the smoke test)

**Interfaces:**
- Consumes: `data.load_cohort, load_comparison, load_boxplot_figure`; `components.finding_summary, format_stats_table, style_stats_table`.
- Produces: full `part3.render()`.

- [ ] **Step 1: Replace `dashboard/sections/part3.py` with the full report**

```python
"""Part 3 · Statistical Analysis: a fixed responder vs non-responder report on
the required cohort (melanoma + miraclib + PBMC), presented as a paper."""
import streamlit as st

from dashboard.data import load_cohort, load_comparison, load_boxplot_figure
from dashboard.components import finding_summary, format_stats_table, style_stats_table

_CONDITION, _TREATMENT, _SAMPLE_TYPE = "melanoma", "miraclib", "PBMC"


def render():
    st.subheader("Part 3 · Statistical analysis: responders vs non-responders")

    st.markdown("#### Question")
    st.write("Do immune cell population relative frequencies differ between "
             "miraclib **responders** and **non-responders** in melanoma PBMC "
             "samples — and could any population predict response to miraclib?")

    cohort = load_cohort(_CONDITION, _TREATMENT, _SAMPLE_TYPE)
    n_subjects = cohort["subject_id"].nunique()
    counts = cohort.drop_duplicates("subject_id")["response"].value_counts()
    st.markdown("#### Cohort")
    st.write(
        f"Melanoma + miraclib + PBMC: **{n_subjects} subjects** "
        f"({int(counts.get('yes', 0))} responders, "
        f"{int(counts.get('no', 0))} non-responders), each sampled at "
        "timepoints 0, 7, and 14.")

    st.markdown("#### Methods")
    st.markdown(
        "- **Primary — linear mixed-effects model** "
        "(`percentage ~ response + C(time)`, random intercept per subject): "
        "uses all timepoints while respecting the repeated measures per subject "
        "(avoids pseudoreplication).\n"
        "- **Cross-check — Mann-Whitney U** on per-subject mean %: "
        "distribution-free, one independent observation per subject.\n"
        "- **Multiple comparisons — Benjamini-Hochberg FDR** across the five "
        "populations (reported as q).\n"
        "- **Effect size — bootstrap 95% CI** of the responder − non-responder "
        "difference in mean %.\n"
        "- A population is a **candidate** if significant under the primary "
        "model, and **confirmed** only if also significant under the "
        "non-parametric cross-check."
    )

    result = load_comparison(_CONDITION, _TREATMENT, _SAMPLE_TYPE)

    st.markdown("#### Results")
    st.info(finding_summary(result))
    st.plotly_chart(
        load_boxplot_figure(_CONDITION, _TREATMENT, _SAMPLE_TYPE),
        use_container_width=True)
    st.markdown("**Per-population statistics**")
    st.caption("Highlighted rows are significant under the primary mixed-effects "
               "model (FDR q < 0.05).")
    st.dataframe(style_stats_table(format_stats_table(result)),
                 use_container_width=True, hide_index=True)

    st.markdown("#### Interpretation")
    confirmed = result.loc[result["concordant"], "population"].tolist()
    candidates = result.loc[result["significant_primary"]
                            & ~result["concordant"], "population"].tolist()
    group_var = float(result["group_var"].abs().max())
    lines = []
    if confirmed:
        lines.append(f"**Confirmed:** {', '.join(confirmed)} differ between "
                     "responders and non-responders under both the mixed model "
                     "and the non-parametric cross-check.")
    if candidates:
        lines.append(f"**Candidate(s):** {', '.join(candidates)} reach "
                     "significance under the primary model but not the "
                     "conservative per-subject test — worth validating, not an "
                     "established biomarker.")
    if not confirmed and not candidates:
        lines.append("No population shows a significant difference between "
                     "responders and non-responders.")
    lines.append(f"The between-subject (random-effect) variance is "
                 f"~{group_var:.2g} across populations, so the repeated "
                 "timepoints add little within-subject correlation here.")
    st.markdown("\n\n".join(lines))
```

- [ ] **Step 2: Remove the now-unused `style_comparison` from `dashboard/components.py`**

Delete the `style_comparison` function (Task 3's `style_stats_table` replaces it; no code references `style_comparison` after the section rewrites).

- [ ] **Step 3: Strengthen the smoke test in `tests/test_dashboard.py`**

Replace `test_app_runs_without_exception` with:

```python
def test_app_runs_and_renders_part3_finding():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
    infos = " ".join(el.value for el in at.info)
    assert "responder" in infos.lower() or "cell population" in infos.lower()
```

- [ ] **Step 4: Run the tests**

Run: `python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS, pristine. The smoke test runs the real fixed-cohort analysis (compare_responders) and confirms the Part 3 finding renders.

- [ ] **Step 5: Full suite**

Run: `python3 -m pytest -v`
Expected: PASS.

- [ ] **Step 6: Launch and eyeball**

Run: `python3 -m streamlit run dashboard/app.py --server.headless true` (or `make dashboard`).
Expected: three part-named tabs; Part 3 reads as a report (Question / Cohort / Methods / Results with finding + boxplots + full stats table + Interpretation) with no filter controls; Part 4 shows breakdowns, the 10206.15 metric, and the degenerate-coverage table. Stop the server after confirming.

- [ ] **Step 7: Commit**

```bash
git add dashboard/sections/part3.py dashboard/components.py tests/test_dashboard.py
git commit -m "feat(dashboard): Part 3 statistical-analysis paper"
```

---

## Self-Review

**Spec coverage (against the redesign design doc):**
- 3 tabs named per instruction part → `app.py` (Task 2). ✓
- Part 2 plain summary table → `part2.py` (Task 2). ✓
- Part 3 fixed paper with methods AND results (finding, boxplots, full stats table always visible, interpretation) → `part3.py` (Task 3). ✓
- Part 4 breakdowns + answer + degenerate-subset callout → `part4.py` + `subset_coverage` (Tasks 1–2). ✓
- Minimal sidebar; results in tabs → `app.py` (Task 2). ✓
- No dead-combo filters (Part 3 fixed) → `part3.py` (Task 3). ✓

**Placeholder scan:** No TBD/TODO in shipped code; the Task 2 `part3.py` stub is explicitly replaced in Task 3. ✓

**Type consistency:** `format_stats_table` output columns match `style_stats_table`'s `"Significant"` check and the test's expected list; `subset_coverage` columns (`condition, treatment, samples, responder_labelled, degenerate`) match the loader and test; `part3.render()` uses `finding_summary`/`format_stats_table`/`style_stats_table` with the `compare_responders` column names (`concordant`, `significant_primary`, `group_var`, `population`). ✓

**Note:** `part2.py` renders the full 52,500-row frequency table via `st.dataframe`, which is virtualized/searchable in Streamlit — acceptable and matches the "show the summary table" requirement.
