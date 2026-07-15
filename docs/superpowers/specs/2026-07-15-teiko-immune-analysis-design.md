# Loblaw Bio — Immune Cell Population Analysis: Overarching Design

**Date:** 2026-07-15
**Status:** Approved design (overarching). Per-part specs to follow.

---

## 1. Context

Bob Loblaw (Loblaw Bio) is running a clinical trial for the drug candidate
**miraclib** and needs to understand how it affects immune cell populations.
We are given `cell-count.csv` — cell counts for five immune populations per
biological sample, plus sample/subject metadata — and must deliver:

1. A SQLite relational database + `load_data.py` loader (Part 1).
2. A per-sample relative-frequency summary table (Part 2).
3. A robust responder-vs-non-responder statistical comparison (Part 3).
4. A baseline-subset exploratory query with breakdown counts (Part 4).
5. A specific computed answer (avg B cells, melanoma males, responders, t=0).
6. An interactive dashboard (public link) surfacing all of the above.
7. A `Makefile` with `setup`, `pipeline`, `dashboard` targets, and a README.

### Data profile (verified from the CSV)

- **10,500 samples**, **3,500 subjects**, **3 projects** (prj1–prj3).
- Conditions: melanoma (5,175), carcinoma (3,903), healthy (1,422).
- Treatments: miraclib (4,695), phauximab (4,383), none (1,422).
- Response: yes (4,611), no (4,467), blank (1,422 — all healthy/`none`).
- Sample types: PBMC (7,500), WB (3,000). Timepoints: 0 / 7 / 14 (3,500 each).
- Sex: F (5,070), M (5,430).
- Each subject contributes 3 samples (one per timepoint).
- **`condition`, `age`, `sex`, `treatment`, `response` are constant per
  subject**; `sample_type` and `time_from_treatment_start` vary per sample.

### CSV → domain column mapping

The instructions' prose names differ from the CSV headers:

| Instructions term | CSV column     |
|-------------------|----------------|
| sample_id         | `sample`       |
| indication        | `condition`    |
| gender            | `sex`          |
| (n/a)             | `project`, `subject`, `age`, `sample_type` |

---

## 2. Architecture

```
CSV ──► load_data.py ──► SQLite (normalized, indexed)
                              │
                              ▼
                    analysis/  (pure functions: conn + params → DataFrame/figure)
                       ┌──────────┴──────────┐
                       ▼                      ▼
                 pipeline.py             dashboard/app.py
           (calls fns, WRITES         (calls SAME fns live,
            outputs/ artifacts)        cached, RENDERS them)
```

Three isolated responsibilities: **compute** (`analysis/`), **persist**
(`pipeline.py`), **render** (`dashboard/app.py`). Aggregation is pushed into
SQL; Python/dashboard consume small result sets only.

---

## 3. Database schema

Normalized, **long-format** cell counts. Five tables:

| Table         | Columns                                                                                          | Grain                         |
|---------------|--------------------------------------------------------------------------------------------------|-------------------------------|
| `projects`    | `project_id` (PK)                                                                                | one row per project           |
| `subjects`    | `subject_id` (PK), `project_id` (FK), `condition`, `age`, `sex`, `treatment`, `response`         | one row per subject           |
| `samples`     | `sample_id` (PK), `subject_id` (FK), `sample_type`, `time_from_treatment_start`                  | one row per sample            |
| `populations` | `population_id` (PK), `name` (UNIQUE)                                                             | lookup (b_cell, cd8_t_cell…)  |
| `cell_counts` | `sample_id` (FK), `population_id` (FK), `count`  — PK `(sample_id, population_id)`                | **one row per sample×pop**    |

Subject-level attributes (`condition`, `age`, `sex`, `treatment`, `response`)
live on `subjects` (no per-row duplication); sample-varying attributes live on
`samples`. `response` is nullable (healthy / `none` treatment).

**Indexes:** all FKs; plus filter columns used by analytics —
`subjects(condition)`, `subjects(treatment)`, `subjects(response)`,
`samples(sample_type)`, `samples(time_from_treatment_start)`, and a composite
`cell_counts(population_id, sample_id)`.

### Rationale: why long-format over wide

The wide alternative (five count columns on `samples`) mirrors the CSV but
fails the scaling requirement in three ways:

1. **New populations require schema migration.** Different projects/assays
   measure different panels; a new marker in wide format is an
   `ALTER TABLE ADD COLUMN` plus sparse NULL columns across all projects. In
   long format it is simply new rows — the schema never changes.
2. **Population becomes a first-class queryable dimension.** Part 2's relative
   frequency is `count / SUM(count) OVER (PARTITION BY sample_id)`; Part 3 is a
   `GROUP BY population`. Both are one-liners in long format; wide format needs
   hardcoded column names and UNPIVOT.
3. **Controlled vocabulary** via the `populations` lookup: one canonical name,
   space savings at scale, and a home for future population metadata (lineage,
   marker definitions) without touching the fact table.

**Cost:** 5× rows (52,500 today; ~1M at "hundreds of projects" scale — trivial
for SQLite) and a pivot to reconstruct a wide view. Acceptable trade for
extensibility.

### Scaling story (for README)

- Push aggregation into SQL (`GROUP BY`, indexed filters); the app consumes
  small aggregates whose size is independent of raw row count.
- The normalized schema ports 1:1 to Postgres/DuckDB; at warehouse scale,
  `cell_counts` is a natural fact table for a star schema and columnar storage.
- Materialized views / a `sample_frequencies` cache can be added if per-sample
  aggregation ever becomes hot — without changing the base schema.

---

## 4. Analysis layer (`analysis/`)

Pure functions: each takes a DB connection + explicit filter params and
**returns** data (a `DataFrame` or a plot figure object). No printing, no file
writes, no Streamlit imports.

| Module          | Responsibility                                                              |
|-----------------|-----------------------------------------------------------------------------|
| `frequencies.py`| Part 2 — per-sample relative frequency table (`sample, total_count, population, count, percentage`). |
| `comparison.py` | Part 3 — responder vs non-responder stats + boxplots.                       |
| `subsets.py`    | Part 4 — baseline melanoma/miraclib/PBMC subset + project/response/sex counts. |
| `questions.py`  | The final computed answer (avg B cells, melanoma males, responders, t=0).    |

---

## 5. Statistical methodology (Part 3 — publication-grade)

**Scope:** melanoma + miraclib + PBMC samples, responders (`response=yes`) vs
non-responders (`response=no`), per cell population (relative frequency %).

- **Primary model — mixed effects:** `frequency ~ response + (1 | subject)`
  (subject random intercept), fit per population. This correctly uses all three
  timepoints while respecting that each subject contributes repeated,
  non-independent samples (avoids pseudoreplication). Report the `response`
  coefficient, its CI, and p-value.
- **Cross-check — Mann-Whitney U** per population (non-parametric, robust to
  skew/outliers) on a per-subject unit, with **Benjamini-Hochberg FDR**
  correction across the five populations. Report raw p and adjusted q.
- **Effect & uncertainty:** medians/means per group, effect size, n per group,
  bootstrap confidence intervals.
- **Visualization:** one boxplot per population (responder vs non-responder)
  with individual points overlaid.
- **Conclusion:** flag populations significant under both the mixed model and
  FDR-adjusted MWU; note any disagreement as a sensitivity caveat.

---

## 6. Dashboard & integration

**Framework:** Streamlit; **hosting:** Streamlit Community Cloud (free public
link from the GitHub repo). `make dashboard` → `streamlit run dashboard/app.py`.

**Integration pattern (load-bearing):** the dashboard calls the `analysis/`
functions **live**, wrapped in `@st.cache_data`, while `pipeline.py` calls the
**same** functions and writes their results to `outputs/`.

- **One source of truth** — on-screen numbers are byte-identical to graded
  artifacts (same code path). No logic drift.
- **Interactivity** — the dashboard passes different filter args (condition,
  treatment, timepoint) to the same functions; the pipeline passes the
  assignment's fixed filters.
- **Testable core** — pure functions returning values, independent of Streamlit
  and the filesystem.

**Rejected alternatives:** (A) dashboard reads static files only — not
interactive, goes stale; (B) dashboard recomputes independently — logic drift.

**Dashboard sections:** (1) data overview + frequency table (Part 2);
(2) responder-vs-non-responder boxplots + stats table (Part 3); (3) baseline
subset explorer + breakdown counts (Part 4); (4) the final computed answer.

---

## 7. Repository layout

```
teiko/
├── load_data.py        # REQUIRED root script: build schema + load CSV → .db
├── pipeline.py         # runs load + all analyses, writes outputs/
├── Makefile            # setup / pipeline / dashboard
├── requirements.txt
├── README.md
├── cell-count.csv
├── cell_counts.db      # generated
├── db/
│   ├── schema.sql      # DDL
│   └── connection.py   # get_connection()
├── analysis/
│   ├── frequencies.py
│   ├── comparison.py
│   ├── subsets.py
│   └── questions.py
├── outputs/            # generated tables (.csv) + plots (.png/.html)
├── dashboard/
│   └── app.py
└── tests/
```

---

## 8. Pipeline, Makefile & deployment

- `make setup` → `pip install -r requirements.txt`.
- `make pipeline` → `python load_data.py && python pipeline.py`. Idempotent:
  rebuilds the DB from scratch and regenerates all `outputs/`.
- `make dashboard` → `streamlit run dashboard/app.py`.
- `load_data.py` runs standalone (`python load_data.py`, no args, no `-m`),
  creating a `.db` file in the repo root, per the requirements.
- Deployment: push repo to GitHub, connect to Streamlit Community Cloud, obtain
  public `https://<app>.streamlit.app` link for the README.

---

## 9. Testing

Unit tests (`tests/`) over the pure analysis functions against a small fixture
DB: frequency percentages sum to 100 per sample; subset filters return expected
counts; the final computed answer matches an independently derived value; stats
functions run and return well-formed results on fixture data.

---

## 10. Requirements traceability

| Requirement                              | Where addressed              |
|------------------------------------------|------------------------------|
| SQLite schema + `load_data.py` (Part 1)  | §3, §7                       |
| Relative-frequency table (Part 2)        | §4 `frequencies.py`          |
| Responder vs non-responder stats (Part 3)| §5 `comparison.py`           |
| Baseline subset + counts (Part 4)        | §4 `subsets.py`              |
| Final computed answer                    | §4 `questions.py`            |
| Interactive dashboard + link             | §6                           |
| Schema rationale + scaling (README)      | §3                           |
| Code-structure rationale (README)        | §2, §7                       |
| `make setup/pipeline/dashboard`          | §8                           |

---

## 11. Next steps

Decompose into per-part specs → implementation plans, in dependency order:

1. **Foundation** — schema (`db/`) + `load_data.py` + connection helper.
2. **Part 2** — frequency analysis + tests.
3. **Part 4** — subset queries + final answer (simple SQL; unblocks dashboard).
4. **Part 3** — statistical comparison (heaviest; `statsmodels`, plots).
5. **Pipeline** — `pipeline.py` + `Makefile` + `outputs/`.
6. **Dashboard** — Streamlit app consuming `analysis/`.
7. **README + deployment** — scaling write-up, run instructions, live link.
```
