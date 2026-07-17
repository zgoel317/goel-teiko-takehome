# Loblaw Bio — Immune Cell Population Analysis

Analysis of how the drug candidate **miraclib** affects immune cell populations
in a clinical trial, using per-sample cell counts for five immune populations
(`b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`).

The project loads `cell-count.csv` into a normalized SQLite database, computes
relative cell-population frequencies, runs a responder-vs-non-responder
statistical comparison, explores a baseline treatment subset, and presents the
results across a three-tab, report-style Streamlit dashboard (one tab per
instruction part).

### Instructions

We are graded via GitHub Codespaces using the `Makefile` targets below.

```bash
make setup      # install all dependencies (from requirements.txt)
make pipeline   # build DB, load data, generate all output tables + plots
make dashboard  # start the local dashboard server
```

- `make setup` — installs all Python dependencies.
- `make pipeline` — runs the full pipeline end-to-end with no manual steps:
  initializes the SQLite database, loads `cell-count.csv` (Part 1), and
  generates all required output tables and plots (Parts 2–4) into `outputs/`.
  The step is idempotent — it rebuilds the database from scratch each run.
- `make dashboard` — starts the local Streamlit dashboard server.

The loader also runs standalone, exactly as required:

```bash
python load_data.py    # creates cell_counts.db in the repository root
```

Requirements: Python 3.11+. Dependencies are pinned in `requirements.txt`; core
libraries are `pandas`, `scipy`, `statsmodels`, `plotly`, and `streamlit`.

### Schema Explanation

The data is stored in a normalized, **long-format** SQLite schema (five tables)
rather than a single wide table mirroring the CSV.

| Table         | Columns                                                                                    | Grain                      |
|---------------|--------------------------------------------------------------------------------------------|----------------------------|
| `projects`    | `project_id` (PK)                                                                          | one row per project        |
| `subjects`    | `subject_id` (PK), `project_id` (FK), `condition`, `age`, `sex`, `treatment`, `response`   | one row per subject         |
| `samples`     | `sample_id` (PK), `subject_id` (FK), `sample_type`, `time_from_treatment_start`            | one row per sample         |
| `populations` | `population_id` (PK), `name` (UNIQUE)                                                       | lookup of cell populations |
| `cell_counts` | `sample_id` (FK), `population_id` (FK), `count` — PK `(sample_id, population_id)`           | one row per sample×pop      |

Subject-level attributes (`condition`, `age`, `sex`, `treatment`, `response`)
are stored once on `subjects` rather than repeated on every sample row;
sample-varying attributes (`sample_type`, `time_from_treatment_start`) live on
`samples`. `response` is nullable (healthy subjects / `none` treatment have no
response). Indexes cover all foreign keys and the columns analytics filter on
(`condition`, `treatment`, `response`, `sample_type`,
`time_from_treatment_start`).

**Why long-format cell counts.** Counts are stored one row per sample×population
in `cell_counts`, not as five columns:

1. **New populations require no schema change** — adding a cell population is an
   `INSERT`, not an `ALTER TABLE ADD COLUMN`. Projects/assays measuring
   different immune panels coexist without sparse NULL columns.
2. **Population is a first-class queryable dimension** — relative frequency is
   `count / SUM(count) OVER (PARTITION BY sample_id)` and per-population
   comparisons are a `GROUP BY`, with no hardcoded column names or UNPIVOT.
3. **Controlled vocabulary** via the `populations` lookup — one canonical name,
   space savings, and a place to attach population metadata later.

The cost (~5× rows, and a pivot to reconstruct a wide view) is negligible —
about 52,500 count rows today.

**How this scales** to hundreds of projects, thousands of samples, and varied
analytics:

- **Aggregation is pushed into SQL** (`GROUP BY`, indexed filters), so the
  application consumes small aggregate result sets whose size is independent of
  raw row count — query and dashboard performance do not degrade as data grows.
- **It is a natural star schema** — `cell_counts` is the fact table;
  `subjects`, `samples`, `populations`, `projects` are dimensions. It ports 1:1
  to Postgres or a columnar warehouse (DuckDB/BigQuery/Snowflake) when SQLite is
  outgrown, with no application changes, because the analysis layer talks to a
  query interface rather than to raw files.
- **Hot per-sample aggregations can be materialized** (e.g. a cached
  `sample_frequencies` view) without altering the base schema.

### Code Structure

```
.
├── load_data.py        # Part 1: builds schema + loads CSV → cell_counts.db (root, no args)
├── pipeline.py         # runs load + all analyses, writes outputs/
├── Makefile            # setup / pipeline / dashboard targets
├── requirements.txt
├── conftest.py         # puts repo root on sys.path for the test suite
├── cell-count.csv      # input data
├── cell_counts.db      # generated SQLite database (git-ignored)
├── .streamlit/
│   └── config.toml     # dashboard theme
├── db/
│   ├── schema.sql      # table definitions (DDL)
│   └── connection.py   # get_connection() helper
├── analysis/           # pure analysis functions (conn + params → DataFrame/figure)
│   ├── frequencies.py  # Part 2: per-sample relative frequency table
│   ├── comparison.py   # Part 3: responder vs non-responder stats + boxplots
│   ├── subsets.py      # Part 4: baseline subset queries + breakdown counts
│   └── questions.py    # the final computed answer
├── outputs/            # generated tables (.csv) and the boxplot (.html)
├── dashboard/          # Streamlit render layer (report style)
│   ├── app.py          # entry point: theme, sidebar, three part-named tabs
│   ├── data.py         # cached, self-building DB connection + cached analysis wrappers
│   ├── components.py   # palette, KPI tiles, table styling, plain-language finding
│   └── sections/       # one report per instruction part
│       ├── part2.py    # cell-type frequencies
│       ├── part3.py    # statistical analysis (paper style)
│       └── part4.py    # subset analysis + degenerate-subset coverage
└── tests/              # unit tests + a Streamlit AppTest smoke test
```

The codebase separates three concerns so each can be developed and tested
independently:

- **Compute** (`analysis/`) — each function takes a database connection plus
  explicit filter parameters and *returns* data (a `DataFrame` or a plot
  figure). It never prints, writes files, or imports Streamlit.
- **Persist** (`pipeline.py`) — calls the analysis functions and writes their
  results to `outputs/`, the reproducible artifacts produced by `make pipeline`.
- **Render** (`dashboard/`) — calls the *same* analysis functions (cached with
  `@st.cache_data`) and displays them across three report-style tabs, one per
  instruction part.

Because the pipeline and the dashboard are two thin consumers of one analysis
library, the numbers on screen are byte-identical to the generated artifacts —
there is no second implementation of any calculation to drift out of sync. The
dashboard also self-builds the database from `cell-count.csv` on first load
(caching the connection), so it runs on a fresh clone without `make pipeline`.

**Analysis performed:**

- **Part 2 — Data overview:** per sample, the total cell count and each
  population's relative frequency (%), as a table with columns
  `sample, total_count, population, count, percentage`.
- **Part 3 — Statistical analysis:** compares population relative frequencies
  between responders and non-responders among **melanoma + miraclib + PBMC**
  samples. The primary test is a linear mixed-effects model
  (`percentage ~ response + C(time)`, random intercept per subject) that uses all
  timepoints while respecting the repeated measures; where the between-subject
  variance is ~0 it reduces to — and is computed as — OLS. It is cross-checked
  with a Mann-Whitney U test on per-subject means, with Benjamini-Hochberg FDR
  correction across the five populations and bootstrap 95% CIs for effect sizes,
  reported alongside per-population boxplots. A population is a *candidate* if
  significant under the mixed model, and *confirmed* only if also significant
  under the non-parametric cross-check.
- **Part 4 — Subset analysis:** melanoma PBMC samples at baseline
  (`time_from_treatment_start = 0`) from miraclib-treated patients, broken down
  by project, by responder/non-responder, and by sex.

**Why this statistical design.** Each choice guards a specific threat to validity:

- *Mixed-effects model (primary).* Each subject contributes three timepoints, so
  the samples are not independent; treating them as independent
  (pseudoreplication) would understate standard errors and inflate false
  positives. A per-subject random intercept models that within-subject
  correlation, so the responder effect is judged against the number of subjects,
  not samples.
- *Mann-Whitney U on per-subject means (cross-check).* Collapsing each subject to
  one value restores independence, and the rank-based test assumes no normality —
  robust to the skew of bounded [0, 100] frequencies. Agreement between the
  parametric and non-parametric routes is the basis of the candidate/confirmed
  rule.
- *Benjamini-Hochberg FDR.* Five populations are tested at once; without
  correction the chance of at least one spurious hit is ~23%. FDR controls that
  while keeping more power than Bonferroni.
- *Bootstrap 95% CIs.* A p-value says whether an effect is distinguishable from
  zero, not whether it is large enough to matter; the CI on the
  responder − non-responder difference makes the magnitude explicit.

**Answer to the posed question** — *Considering melanoma males of all sample and
treatment types, what is the average number of B cells for responders at
time = 0?*
**`10206.15`** (mean B-cell count over the n=485 qualifying samples), computed by
`analysis/questions.py::avg_bcells_melanoma_male_responders_at_baseline`.

### Dashboard

**Live dashboard:** https://goel-teiko-takehome-dauceseepmxvm8mmxqvbiv.streamlit.app/

The dashboard is deployed on Streamlit Community Cloud (main file
`dashboard/app.py`). It self-builds the SQLite database from the committed
`cell-count.csv` on first load, so no separate data step is needed. When
deploying, select **Python 3.11** in the app's advanced settings — newer Python
versions currently lack prebuilt wheels for some dependencies.
