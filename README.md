# Loblaw Bio — Immune Cell Population Analysis

Analysis of how the drug candidate **miraclib** affects immune cell populations
in a clinical trial, using per-sample cell counts for five immune populations
(`b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`).

The project loads `cell-count.csv` into a normalized SQLite database, computes
relative cell-population frequencies, runs a responder-vs-non-responder
statistical comparison, explores a baseline treatment subset, and surfaces
everything in an interactive Streamlit dashboard.

### Instructions

We are graded via GitHub Codespaces using the `Makefile` targets below.

```bash
make setup      # install all dependencies (from requirements.txt)
make pipeline   # build DB, load data, generate all output tables + plots
make dashboard  # start the local server for the interactive dashboard
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
├── cell-count.csv      # input data
├── cell_counts.db      # generated SQLite database (git-ignored)
├── db/
│   ├── schema.sql      # table definitions (DDL)
│   └── connection.py   # get_connection() helper
├── analysis/           # pure analysis functions (conn + params → DataFrame/figure)
│   ├── frequencies.py  # Part 2: per-sample relative frequency table
│   ├── comparison.py   # Part 3: responder vs non-responder stats + boxplots
│   ├── subsets.py      # Part 4: baseline subset queries + breakdown counts
│   └── questions.py    # the final computed answer
├── outputs/            # generated tables (.csv) and plots (.png/.html)
├── dashboard/
│   └── app.py          # Streamlit dashboard
└── tests/              # unit tests over the analysis functions
```

The codebase separates three concerns so each can be developed and tested
independently:

- **Compute** (`analysis/`) — each function takes a database connection plus
  explicit filter parameters and *returns* data (a `DataFrame` or a plot
  figure). It never prints, writes files, or imports Streamlit.
- **Persist** (`pipeline.py`) — calls the analysis functions and writes their
  results to `outputs/`, the reproducible artifacts produced by `make pipeline`.
- **Render** (`dashboard/app.py`) — calls the *same* analysis functions live
  (cached with `@st.cache_data`) and displays them.

Because the pipeline and the dashboard are two thin consumers of one analysis
library, the numbers on screen are byte-identical to the generated artifacts —
there is no second implementation of any calculation to drift out of sync, and
the dashboard gains interactivity for free by passing different filter arguments
to the same functions.

**Analysis performed:**

- **Part 2 — Data overview:** per sample, the total cell count and each
  population's relative frequency (%), as a table with columns
  `sample, total_count, population, count, percentage`.
- **Part 3 — Statistical analysis:** compares population relative frequencies
  between responders and non-responders among **melanoma + miraclib + PBMC**
  samples, using a mixed-effects model (`frequency ~ response + (1 | subject)`)
  that respects repeated timepoints per subject, cross-checked with Mann-Whitney
  U tests corrected for multiple comparisons (Benjamini-Hochberg FDR), reported
  alongside per-population boxplots.
- **Part 4 — Subset analysis:** melanoma PBMC samples at baseline
  (`time_from_treatment_start = 0`) from miraclib-treated patients, broken down
  by project, by responder/non-responder, and by sex.

**Answer to the posed question** — *Considering melanoma males of all sample and
treatment types, what is the average number of B cells for responders at
time = 0?*
**`10206.15`** (mean B-cell count over the n=485 qualifying samples), computed by
`analysis/questions.py::avg_bcells_melanoma_male_responders_at_baseline`.

The full overarching design and rationale live in
[`docs/superpowers/specs/2026-07-15-teiko-immune-analysis-design.md`](docs/superpowers/specs/2026-07-15-teiko-immune-analysis-design.md).

### Dashboard

<!-- TODO: Streamlit Community Cloud URL once deployed -->
`https://<app>.streamlit.app`
