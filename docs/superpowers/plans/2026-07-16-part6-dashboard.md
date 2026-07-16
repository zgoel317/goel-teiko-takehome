# Part 6: Streamlit Dashboard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An interactive Streamlit dashboard (3 tabs) surfacing Parts 2–4 plus the layered statistical analysis, reading the existing `analysis/` functions live, deployable to Streamlit Community Cloud.

**Architecture:** A thin render layer. `dashboard/data.py` provides a cached, self-building DB connection and cached analysis wrappers; `dashboard/components.py` holds the palette, plot config, KPI/dataframe helpers, and the pure `finding_summary`; `dashboard/sections/*` render the three tabs; `dashboard/app.py` wires theme, sidebar, and tabs. No analysis logic is re-implemented.

**Tech Stack:** Python 3.11+, Streamlit, pandas, plotly.

## Global Constraints

- Render layer only: call existing `analysis/` functions; never re-implement analysis math.
- **Self-building DB:** `ensure_database` builds `cell_counts.db` from `cell-count.csv` if missing (so the deployed app works with no `make pipeline`).
- Analysis calls cached with `@st.cache_data` (keyed on filter strings); connection cached with `@st.cache_resource`.
- Consistent Plotly palette across every chart (`PALETTE` in `components.py`), applied via `px.defaults`.
- Every tab guards empty/degenerate cohorts with a friendly message, never a traceback.
- Response tab defaults to condition=melanoma, treatment=miraclib, sample_type=PBMC (the required view).
- Files stay focused per the decomposition in the design doc.

## File Structure

- `.streamlit/config.toml` — light theme.
- `dashboard/__init__.py`, `dashboard/sections/__init__.py` — empty package markers.
- `dashboard/data.py` — `ensure_database`, cached `get_connection`, cached loaders.
- `dashboard/components.py` — `PALETTE`, `configure_plotly`, `kpi_row`, `style_comparison`, `finding_summary`.
- `dashboard/sections/overview.py`, `response.py`, `subset.py` — tab renderers (`render()`).
- `dashboard/app.py` — entry point.
- `tests/test_dashboard.py` — `ensure_database` + `finding_summary` units; `AppTest` smoke test.

---

### Task 1: Theme + data layer + components

**Files:**
- Create: `.streamlit/config.toml`, `dashboard/__init__.py`, `dashboard/sections/__init__.py`, `dashboard/data.py`, `dashboard/components.py`
- Test: `tests/test_dashboard.py`

**Interfaces:**
- Consumes: `load_data.build_database, CSV_PATH, DB_PATH`; `db.connection.get_connection`; the `analysis/` functions.
- Produces:
  - `dashboard.data.ensure_database(db_path=DB_PATH, csv_path=CSV_PATH) -> Path` (builds only if missing).
  - `dashboard.data.get_connection()` (cached), and cached loaders: `load_sample_metadata`, `load_frequencies`, `load_frequencies_annotated`, `load_dataset_summary`, `load_comparison`, `load_cohort`, `load_baseline`, `load_boxplot_figure`, `load_final_answer`.
  - `dashboard.components.PALETTE`, `configure_plotly()`, `kpi_row(items)`, `style_comparison(df)`, `finding_summary(comparison, alpha=0.05) -> str`.

- [ ] **Step 1: Create package markers and theme**

```bash
: > dashboard/__init__.py
: > dashboard/sections/__init__.py
mkdir -p .streamlit
```

`.streamlit/config.toml`:

```toml
[theme]
primaryColor = "#0072B2"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F5F7FA"
textColor = "#1A1A1A"
font = "sans serif"

[server]
headless = true
```

- [ ] **Step 2: Write the failing test `tests/test_dashboard.py`**

```python
import pandas as pd

from dashboard.data import ensure_database
from dashboard.components import finding_summary


def test_ensure_database_builds_when_missing(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(
        "project,subject,condition,age,sex,treatment,response,sample,sample_type,"
        "time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte\n"
        "prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100\n"
    )
    db = tmp_path / "cell_counts.db"
    assert not db.exists()
    ensure_database(db_path=db, csv_path=csv)
    assert db.exists()


def test_ensure_database_reuses_existing(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(
        "project,subject,condition,age,sex,treatment,response,sample,sample_type,"
        "time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte\n"
        "prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100\n"
    )
    db = tmp_path / "cell_counts.db"
    ensure_database(db_path=db, csv_path=csv)
    mtime = db.stat().st_mtime_ns
    ensure_database(db_path=db, csv_path=csv)  # must NOT rebuild
    assert db.stat().st_mtime_ns == mtime


def _comparison(significant, concordant, mean_diff=0.64, q_mixed=0.02, q_mw=0.06):
    return pd.DataFrame([
        {"population": "cd4_t_cell", "mean_diff": mean_diff, "q_mixed": q_mixed,
         "q_mw": q_mw, "significant_primary": significant, "concordant": concordant},
        {"population": "b_cell", "mean_diff": -0.2, "q_mixed": 0.4, "q_mw": 0.5,
         "significant_primary": False, "concordant": False},
    ])


def test_finding_summary_none_significant():
    text = finding_summary(_comparison(significant=False, concordant=False))
    assert "No cell population" in text


def test_finding_summary_candidate_non_concordant():
    text = finding_summary(_comparison(significant=True, concordant=False))
    assert "cd4_t_cell" in text
    assert "candidate" in text.lower()
    assert "higher in responders" in text.lower()


def test_finding_summary_concordant_confirmed():
    text = finding_summary(_comparison(significant=True, concordant=True))
    assert "cd4_t_cell" in text
    assert "confirmed" in text.lower()
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `cd /Users/zoyagoel/teiko && python3 -m pytest tests/test_dashboard.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'dashboard.data'`.

- [ ] **Step 4: Create `dashboard/components.py`**

```python
"""Shared dashboard presentation helpers: palette, plot config, KPI tiles,
dataframe styling, and the plain-language findings summary."""
import plotly.express as px
import streamlit as st

# Okabe-Ito colour-blind-safe categorical palette (5 populations / 2 groups).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]


def configure_plotly():
    """Apply the shared palette + a clean template to every Plotly Express chart."""
    px.defaults.color_discrete_sequence = PALETTE
    px.defaults.template = "plotly_white"


def kpi_row(items):
    """Render a row of st.metric tiles. items: list of (label, value)."""
    cols = st.columns(len(items))
    for col, (label, value) in zip(cols, items):
        col.metric(label, value)


def style_comparison(df):
    """Style the Part 3 statistics table: highlight significant_primary rows."""
    def _highlight(row):
        color = "background-color: #FFF3CD" if row.get("significant_primary") else ""
        return [color] * len(row)
    return df.style.apply(_highlight, axis=1).format(precision=3)


def finding_summary(comparison, alpha=0.05):
    """Plain-language markdown summary of the responder comparison result."""
    primary = comparison[comparison["significant_primary"]]
    if primary.empty:
        return ("**No cell population** shows a statistically significant "
                "difference in relative frequency between responders and "
                f"non-responders (all FDR-adjusted q ≥ {alpha:.2f}).")
    lines = []
    for _, r in primary.iterrows():
        direction = "higher" if r["mean_diff"] > 0 else "lower"
        line = (f"**{r['population']}** is {direction} in responders "
                f"({r['mean_diff']:+.2f} pp; mixed-model q={r['q_mixed']:.3f})")
        if r["concordant"]:
            line += (" — **confirmed** by both the mixed model and the "
                     "non-parametric test.")
        else:
            line += (f" — but **not concordant** with the conservative "
                     f"per-subject test (q={r['q_mw']:.3f}); treat as a "
                     f"**candidate**, not confirmed.")
        lines.append(line)
    return "  \n".join(lines)
```

- [ ] **Step 5: Create `dashboard/data.py`**

```python
"""Cached data access for the dashboard: a self-building DB connection and
cached wrappers over the analysis functions."""
from pathlib import Path

import pandas as pd
import streamlit as st

from load_data import build_database, CSV_PATH, DB_PATH
from db.connection import get_connection as _open_connection
from analysis.frequencies import sample_frequencies
from analysis.comparison import (
    compare_responders,
    cohort_frequencies,
    responder_boxplots,
)
from analysis.subsets import baseline_subset
from analysis.questions import avg_bcells_melanoma_male_responders_at_baseline

_META_QUERY = """
SELECT sa.sample_id                 AS sample,
       sa.subject_id                AS subject_id,
       su.project_id                AS project,
       su.condition                 AS condition,
       su.treatment                 AS treatment,
       sa.sample_type               AS sample_type,
       sa.time_from_treatment_start AS time,
       su.response                  AS response,
       su.sex                       AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
"""


def ensure_database(db_path=DB_PATH, csv_path=CSV_PATH):
    """Build the SQLite DB from the CSV only if it does not already exist."""
    db_path = Path(db_path)
    if not db_path.exists():
        build_database(csv_path=csv_path, db_path=db_path)
    return db_path


@st.cache_resource
def get_connection():
    ensure_database()
    return _open_connection()


@st.cache_data
def load_sample_metadata():
    return pd.read_sql_query(_META_QUERY, get_connection())


@st.cache_data
def load_frequencies():
    return sample_frequencies(get_connection())


@st.cache_data
def load_frequencies_annotated():
    return load_frequencies().merge(load_sample_metadata(), on="sample", how="left")


@st.cache_data
def load_dataset_summary():
    m = load_sample_metadata()
    return {
        "samples": int(m["sample"].nunique()),
        "subjects": int(m["subject_id"].nunique()),
        "projects": int(m["project"].nunique()),
    }


@st.cache_data
def load_comparison(condition, treatment, sample_type):
    return compare_responders(get_connection(), condition=condition,
                              treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_cohort(condition, treatment, sample_type):
    return cohort_frequencies(get_connection(), condition=condition,
                              treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_boxplot_figure(condition, treatment, sample_type):
    return responder_boxplots(load_cohort(condition, treatment, sample_type))


@st.cache_data
def load_baseline(condition, treatment, sample_type):
    return baseline_subset(get_connection(), condition=condition,
                           treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_final_answer():
    return avg_bcells_melanoma_male_responders_at_baseline(get_connection())
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS (5 tests), pristine.

- [ ] **Step 7: Commit**

```bash
git add .streamlit dashboard tests/test_dashboard.py
git commit -m "feat(part6): dashboard data layer, components, theme"
```

---

### Task 2: App shell + Overview + Subset tabs

**Files:**
- Create: `dashboard/sections/overview.py`, `dashboard/sections/subset.py`, `dashboard/sections/response.py` (stub), `dashboard/app.py`
- Modify: `tests/test_dashboard.py` (add AppTest smoke test)

**Interfaces:**
- Consumes: `dashboard.data` loaders; `dashboard.components`.
- Produces: `overview.render()`, `subset.render()`, `response.render()` (stub for now), and a runnable `dashboard/app.py`.

- [ ] **Step 1: Create `dashboard/sections/overview.py`**

```python
"""Overview tab (Part 2): dataset KPIs + relative-frequency table + composition."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_dataset_summary, load_frequencies_annotated
from dashboard.components import kpi_row

_ALL = "All"


def _select(label, values):
    return st.selectbox(label, [_ALL] + sorted(values))


def render():
    st.subheader("Data overview — relative cell-population frequencies")
    summary = load_dataset_summary()
    kpi_row([("Samples", summary["samples"]),
             ("Subjects", summary["subjects"]),
             ("Projects", summary["projects"])])

    data = load_frequencies_annotated()
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        condition = _select("Condition", data["condition"].dropna().unique())
    with c2:
        treatment = _select("Treatment", data["treatment"].dropna().unique())
    with c3:
        sample_type = _select("Sample type", data["sample_type"].dropna().unique())
    with c4:
        populations = st.multiselect("Populations",
                                     sorted(data["population"].unique()))

    view = data
    if condition != _ALL:
        view = view[view["condition"] == condition]
    if treatment != _ALL:
        view = view[view["treatment"] == treatment]
    if sample_type != _ALL:
        view = view[view["sample_type"] == sample_type]
    if populations:
        view = view[view["population"].isin(populations)]

    if view.empty:
        st.info("No samples match these filters.")
        return

    mean_comp = (view.groupby("population", as_index=False)["percentage"].mean()
                 .sort_values("percentage", ascending=False))
    fig = px.bar(mean_comp, x="population", y="percentage",
                 title="Mean relative frequency by population (%)",
                 labels={"percentage": "Mean %", "population": "Population"})
    st.plotly_chart(fig, use_container_width=True)

    st.dataframe(
        view[["sample", "total_count", "population", "count", "percentage"]]
        .reset_index(drop=True),
        use_container_width=True, hide_index=True,
    )
```

- [ ] **Step 2: Create `dashboard/sections/subset.py`**

```python
"""Subset Explorer tab (Part 4): baseline subset breakdowns + posed answer."""
import plotly.express as px
import streamlit as st

from dashboard.data import (
    load_baseline,
    load_sample_metadata,
    load_final_answer,
)
from dashboard.components import kpi_row
from analysis.subsets import (
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)


def render():
    st.subheader("Baseline subset explorer (time = 0)")
    meta = load_sample_metadata()
    c1, c2 = st.columns(2)
    with c1:
        condition = st.selectbox("Condition",
                                 sorted(meta["condition"].dropna().unique()),
                                 index=_default_index(meta["condition"], "melanoma"))
    with c2:
        treatment = st.selectbox("Treatment",
                                 sorted(meta["treatment"].dropna().unique()),
                                 index=_default_index(meta["treatment"], "miraclib"))

    subset = load_baseline(condition, treatment, "PBMC")
    if subset.empty:
        st.info("No baseline PBMC samples match this condition/treatment.")
    else:
        kpi_row([("Samples", len(subset)),
                 ("Subjects", subset["subject_id"].nunique())])
        col1, col2, col3 = st.columns(3)
        with col1:
            _bar(samples_per_project(subset), "project", "n_samples",
                 "Samples per project")
        with col2:
            _bar(subjects_by_response(subset), "response", "n_subjects",
                 "Subjects by response")
        with col3:
            _bar(subjects_by_sex(subset), "sex", "n_subjects", "Subjects by sex")

    st.divider()
    answer = load_final_answer()
    st.metric(
        "Avg B cells — melanoma male responders at baseline "
        "(all sample & treatment types)",
        f"{answer:.2f}",
    )


def _default_index(series, value):
    options = sorted(series.dropna().unique())
    return options.index(value) if value in options else 0


def _bar(df, x, y, title):
    fig = px.bar(df, x=x, y=y, title=title)
    st.plotly_chart(fig, use_container_width=True)
```

- [ ] **Step 3: Create `dashboard/sections/response.py` (stub)**

```python
"""Response Analysis tab (Part 3). Full implementation in the next task."""
import streamlit as st


def render():
    st.subheader("Response analysis")
    st.info("Response analysis is being wired up.")
```

- [ ] **Step 4: Create `dashboard/app.py`**

```python
"""Loblaw Bio immune-cell analysis dashboard."""
import os
import sys

# `streamlit run dashboard/app.py` puts dashboard/ on sys.path, not the repo
# root, so project imports (load_data, analysis.*, dashboard.*) would fail both
# locally (`make dashboard`) and on Streamlit Community Cloud. Put the repo root
# (this file's grandparent) on the path before any project import.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st

st.set_page_config(page_title="Loblaw Bio — Immune Cell Analysis",
                   layout="wide")

from dashboard.components import configure_plotly
from dashboard.sections import overview, response, subset

configure_plotly()

st.sidebar.title("Loblaw Bio")
st.sidebar.caption("Immune cell population analysis — miraclib trial")
with st.sidebar.expander("About & methodology"):
    st.markdown(
        "Relative cell-population frequencies from `cell-count.csv`. "
        "Response analysis compares miraclib responders vs non-responders "
        "(melanoma, PBMC) using a linear mixed-effects model (primary), a "
        "Mann-Whitney U cross-check on per-subject means, Benjamini-Hochberg "
        "FDR correction, and bootstrap effect-size CIs."
    )

st.title("Immune Cell Population Analysis")

tab_overview, tab_response, tab_subset = st.tabs(
    ["Overview", "Response Analysis", "Subset Explorer"])
with tab_overview:
    overview.render()
with tab_response:
    response.render()
with tab_subset:
    subset.render()
```

- [ ] **Step 5: Add the AppTest smoke test to `tests/test_dashboard.py`**

Append:

```python
def test_app_runs_without_exception():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
```

- [ ] **Step 6: Run the tests**

Run: `python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS (6 tests). The smoke test builds/uses the real DB and renders all tabs (response stub) with no exception.

- [ ] **Step 7: Run the full suite**

Run: `python3 -m pytest -v`
Expected: PASS (all parts).

- [ ] **Step 8: Commit**

```bash
git add dashboard tests/test_dashboard.py
git commit -m "feat(part6): app shell with Overview and Subset tabs"
```

---

### Task 3: Response Analysis tab (layered statistics)

**Files:**
- Modify: `dashboard/sections/response.py` (replace stub with full implementation)
- Modify: `tests/test_dashboard.py` (strengthen smoke test)

**Interfaces:**
- Consumes: `dashboard.data.load_sample_metadata, load_cohort, load_comparison, load_boxplot_figure`; `dashboard.components.finding_summary, style_comparison`.
- Produces: full `response.render()`.

- [ ] **Step 1: Replace `dashboard/sections/response.py` with the full implementation**

```python
"""Response Analysis tab (Part 3): layered responder vs non-responder stats."""
import streamlit as st

from dashboard.data import (
    load_sample_metadata,
    load_cohort,
    load_comparison,
    load_boxplot_figure,
)
from dashboard.components import finding_summary, style_comparison

_MIN_PER_GROUP = 3


def _default_index(series, value):
    options = sorted(series.dropna().unique())
    return options.index(value) if value in options else 0


def render():
    st.subheader("Responders vs non-responders")
    meta = load_sample_metadata()
    c1, c2, c3 = st.columns(3)
    with c1:
        condition = st.selectbox("Condition",
                                 sorted(meta["condition"].dropna().unique()),
                                 index=_default_index(meta["condition"], "melanoma"))
    with c2:
        treatment = st.selectbox("Treatment",
                                 sorted(meta["treatment"].dropna().unique()),
                                 index=_default_index(meta["treatment"], "miraclib"))
    with c3:
        sample_type = st.selectbox("Sample type",
                                   sorted(meta["sample_type"].dropna().unique()),
                                   index=_default_index(meta["sample_type"], "PBMC"))

    cohort = load_cohort(condition, treatment, sample_type)
    if cohort.empty:
        st.info("No responders/non-responders match this cohort.")
        return
    per_group = cohort.drop_duplicates("subject_id")["response"].value_counts()
    if set(per_group.index) != {"yes", "no"} or per_group.min() < _MIN_PER_GROUP:
        st.warning(
            f"Not enough subjects to compare: need at least {_MIN_PER_GROUP} "
            "responders and non-responders in this cohort.")
        return

    result = load_comparison(condition, treatment, sample_type)

    st.markdown("#### Finding")
    st.info(finding_summary(result))

    st.plotly_chart(load_boxplot_figure(condition, treatment, sample_type),
                    use_container_width=True)

    with st.expander("Full statistics table"):
        st.caption("Highlighted rows are significant under the primary "
                   "mixed-effects model (FDR q < 0.05).")
        st.dataframe(style_comparison(result), use_container_width=True)

    with st.expander("Methodology & assumptions"):
        st.markdown(
            "- **Primary — linear mixed-effects model** "
            "(`frequency ~ response + time`, random intercept per subject): "
            "uses all timepoints while respecting repeated measures.\n"
            "- **Cross-check — Mann-Whitney U** on per-subject means: "
            "distribution-free, independent units.\n"
            "- **Multiple comparisons — Benjamini-Hochberg FDR** across the "
            "five populations.\n"
            "- **Effect size — bootstrap 95% CI** of the responder − "
            "non-responder difference.\n"
            "- A population is a **candidate** if significant under the mixed "
            "model, and **confirmed** only if also significant under the "
            "non-parametric test."
        )
```

- [ ] **Step 2: Strengthen the smoke test in `tests/test_dashboard.py`**

Replace `test_app_runs_without_exception` with:

```python
def test_app_runs_and_renders_response_finding():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
    # The Response tab computes and renders a plain-language finding via st.info.
    infos = " ".join(el.value for el in at.info)
    assert "responder" in infos.lower() or "cell population" in infos.lower()
```

- [ ] **Step 3: Run the tests**

Run: `python3 -m pytest tests/test_dashboard.py -v`
Expected: PASS (6 tests), pristine. The smoke test now exercises the full Response tab (real `compare_responders`) with no exception and confirms a finding rendered.

- [ ] **Step 4: Run the full suite**

Run: `python3 -m pytest -v`
Expected: PASS (all parts).

- [ ] **Step 5: Launch the app and verify manually**

Run: `python3 -m streamlit run dashboard/app.py --server.headless true` (or `make dashboard`).
Expected: app starts; Overview shows KPIs + frequency table; Response shows the cd4_t_cell candidate finding + boxplots + expandable stats; Subset shows breakdowns + the 10206.15 metric. Stop the server after confirming (Ctrl-C).

- [ ] **Step 6: Commit**

```bash
git add dashboard/sections/response.py tests/test_dashboard.py
git commit -m "feat(part6): Response Analysis tab with layered statistics"
```

---

## Self-Review

**Spec coverage (against dashboard design + submission):**
- Part 2 table + interactivity → Overview tab (Task 2). ✓
- Part 3 boxplots + statistics, layered → Response tab (Task 3), `finding_summary` (Task 1). ✓
- Part 4 subset breakdowns + final answer → Subset tab (Task 2). ✓
- Self-building DB for deployment → `ensure_database` (Task 1), tested. ✓
- Consistent palette / theme → `PALETTE`, `configure_plotly`, `config.toml` (Task 1). ✓
- Empty/degenerate cohort guards → each `render()` (Tasks 2–3). ✓
- Testable logic + smoke test → `ensure_database`/`finding_summary` units + `AppTest` (Tasks 1–3). ✓

**Placeholder scan:** No TBD/TODO in shipped code. The Task 2 `response.py` stub is explicitly replaced in Task 3. ✓

**Type consistency:** Loader names in `data.py` match imports in the sections; `finding_summary`/`style_comparison`/`kpi_row`/`configure_plotly` signatures match their call sites; `compare_responders` output columns (`significant_primary`, `concordant`, `mean_diff`, `q_mixed`, `q_mw`, `population`) match `finding_summary` and `style_comparison` usage. ✓

**Note:** `load_sample_metadata` uses a small metadata query local to `data.py` (annotation for display) — a deliberate, minor duplication of a SELECT rather than reaching into `analysis.comparison`'s private query, keeping the tested analysis module untouched.
