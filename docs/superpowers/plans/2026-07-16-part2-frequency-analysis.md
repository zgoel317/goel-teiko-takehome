# Part 2: Frequency Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute, for each sample, the total cell count and each population's relative frequency (%), as a tidy table.

**Architecture:** A pure function `sample_frequencies(conn)` in a new `analysis/` package. Aggregation is done in SQL with a window function (`SUM(count) OVER (PARTITION BY sample_id)`) per the design's "push aggregation into SQL" principle; the function returns a pandas DataFrame. No printing, no file writes, no Streamlit — the pipeline and dashboard (later parts) consume this function.

**Tech Stack:** Python 3.11+, sqlite3 (stdlib), pandas.

## Global Constraints

- Analysis functions are pure: take a DB connection (and explicit params), **return** a DataFrame; never print, write files, or import Streamlit.
- Output table columns, in this exact order and naming: `sample`, `total_count`, `population`, `count`, `percentage`.
  - `sample` = the sample id (the `sample` column from `cell-count.csv`, stored as `samples.sample_id`).
  - `total_count` = sum of all five population counts for that sample.
  - `population` = the cell population name (`b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`).
  - `count` = that population's cell count in that sample.
  - `percentage` = `100.0 * count / total_count` (relative frequency as a percentage), full float precision (rounding is a display concern for later consumers).
- One row per (sample, population). On the real data that is 10,500 samples × 5 = **52,500 rows**.
- Percentages within a sample sum to 100.
- Depends on Part 1: `db.connection.get_connection`, the `cell_counts` / `samples` / `populations` tables, and `load_data.build_database` (used to build fixture DBs in tests).

---

## File Structure

- `analysis/__init__.py` — empty; makes `analysis` an importable package.
- `analysis/frequencies.py` — the `sample_frequencies(conn)` function and its SQL. One clear responsibility: the Part 2 frequency table.
- `tests/test_frequencies.py` — builds a small fixture DB via `build_database`, then asserts columns, shape, totals, exact percentages, and the sum-to-100 property.

---

### Task 1: `sample_frequencies` frequency table

**Files:**
- Create: `analysis/__init__.py`
- Create: `analysis/frequencies.py`
- Test: `tests/test_frequencies.py`

**Interfaces:**
- Consumes (from Part 1): `db.connection.get_connection(db_path) -> sqlite3.Connection`; `load_data.build_database(csv_path, db_path, schema_path) -> Path`; tables `cell_counts(sample_id, population_id, count)`, `samples(sample_id, ...)`, `populations(population_id, name)`.
- Produces: `analysis.frequencies.sample_frequencies(conn) -> pandas.DataFrame` with columns `[sample, total_count, population, count, percentage]`, one row per (sample, population), ordered by `sample_id` then `population_id`.

- [ ] **Step 1: Create the `analysis` package marker**

```bash
: > analysis/__init__.py
```

(Create the `analysis/` directory if it does not exist.)

- [ ] **Step 2: Write the failing test `tests/test_frequencies.py`**

```python
from db.connection import get_connection
from load_data import build_database
from analysis.frequencies import sample_frequencies

FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100
prj1,s0,melanoma,50,M,miraclib,yes,smp1,PBMC,7,1,2,3,4,90
prj2,s1,healthy,40,F,none,,smp2,WB,0,5,5,5,5,5
"""


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_columns_and_row_count(tmp_path):
    df = sample_frequencies(_conn(tmp_path))
    assert list(df.columns) == [
        "sample", "total_count", "population", "count", "percentage"]
    assert len(df) == 15  # 3 samples x 5 populations


def test_total_count_and_exact_percentages(tmp_path):
    df = sample_frequencies(_conn(tmp_path))
    smp0 = df[df["sample"] == "smp0"]
    # 10 + 20 + 30 + 40 + 100 = 200
    assert set(smp0["total_count"]) == {200}
    b_cell = smp0[smp0["population"] == "b_cell"].iloc[0]
    assert b_cell["count"] == 10
    assert b_cell["percentage"] == 5.0        # 10/200*100
    monocyte = smp0[smp0["population"] == "monocyte"].iloc[0]
    assert monocyte["percentage"] == 50.0     # 100/200*100


def test_percentages_sum_to_100_per_sample(tmp_path):
    df = sample_frequencies(_conn(tmp_path))
    sums = df.groupby("sample")["percentage"].sum()
    assert all(abs(v - 100.0) < 1e-9 for v in sums)


def test_one_row_per_sample_population(tmp_path):
    df = sample_frequencies(_conn(tmp_path))
    assert not df.duplicated(subset=["sample", "population"]).any()
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `cd /Users/zoyagoel/teiko && python -m pytest tests/test_frequencies.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'analysis.frequencies'` (or ImportError for `sample_frequencies`).

- [ ] **Step 4: Write `analysis/frequencies.py`**

```python
"""Part 2: per-sample relative frequency of each immune cell population."""
import pandas as pd

# Relative frequency computed in SQL: the per-sample total is a window sum over
# the sample's rows, so each population's percentage is count / total * 100.
FREQUENCY_QUERY = """
SELECT
    s.sample_id                                    AS sample,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
    p.name                                         AS population,
    cc.count                                       AS count,
    100.0 * cc.count
        / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts cc
JOIN samples s     ON s.sample_id = cc.sample_id
JOIN populations p ON p.population_id = cc.population_id
ORDER BY s.sample_id, p.population_id
"""


def sample_frequencies(conn):
    """Return the per-sample relative-frequency table.

    Columns: sample, total_count, population, count, percentage.
    One row per (sample, population); percentages within a sample sum to 100.
    """
    return pd.read_sql_query(FREQUENCY_QUERY, conn)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_frequencies.py -v`
Expected: PASS (4 tests), pristine output.

- [ ] **Step 6: Run the full suite (no regressions from Part 1)**

Run: `python -m pytest -v`
Expected: PASS (all Part 1 + Part 2 tests).

- [ ] **Step 7: Verify against the real database**

First ensure the real DB exists (idempotent): `python load_data.py`
Then run:

```bash
python -c "
from db.connection import get_connection
from analysis.frequencies import sample_frequencies
df = sample_frequencies(get_connection())
print('rows:', len(df))
sums = df.groupby('sample')['percentage'].sum()
print('max abs deviation from 100:', float((sums - 100).abs().max()))
print(df.head(5).to_string(index=False))
"
```

Expected: `rows: 52500`; max abs deviation from 100 is ~0 (< 1e-9); the head shows one sample's five populations with `total_count` equal across its rows and ascending `percentage`/`population` order. If rows != 52500 or the deviation is non-trivial, stop and investigate before committing.

- [ ] **Step 8: Commit**

```bash
git add analysis/ tests/test_frequencies.py
git commit -m "feat(part2): per-sample relative frequency table"
```

---

## Self-Review

**Spec coverage (against overarching design §4 `frequencies.py` and instructions Part 2):**
- Table with columns `sample, total_count, population, count, percentage` → Task 1 `FREQUENCY_QUERY` + `test_columns_and_row_count`. ✓
- `total_count` = sum across the five populations per sample → window `SUM` + `test_total_count_and_exact_percentages`. ✓
- `percentage` = relative frequency as a percentage → `100.0 * count / total` + exact-value test. ✓
- One row per population per sample → `test_one_row_per_sample_population`; 52,500 rows on real data → Step 7. ✓
- Pure function returning a DataFrame, SQL-side aggregation → `sample_frequencies` (no I/O, window function). ✓

**Placeholder scan:** No TBD/TODO; every code and test step is complete. ✓

**Type consistency:** `sample_frequencies(conn) -> DataFrame` and the column list `[sample, total_count, population, count, percentage]` are identical across the interface block, the SQL `AS` aliases, and the tests. Fixture math (smp0 total 200, b_cell 5.0%, monocyte 50.0%) is internally consistent. ✓

**Note:** `percentage` is full-precision float by design; any rounding (e.g. the 2-decimal display) is deferred to the pipeline/dashboard consumers, not this function.
