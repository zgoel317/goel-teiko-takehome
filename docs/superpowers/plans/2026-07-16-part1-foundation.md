# Part 1: Foundation (Schema + Data Loading) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a normalized SQLite database from `cell-count.csv` via a root-level `load_data.py` that runs with no arguments.

**Architecture:** A `db/` package holds the schema DDL (`schema.sql`) and a connection helper (`connection.py`) that enables foreign keys. `load_data.py` at the repo root reads the CSV with pandas, reshapes the five wide population columns into a long `cell_counts` table, and populates five tables in dependency order. All paths are resolved relative to the script location so it runs from anywhere.

**Tech Stack:** Python 3.11+, sqlite3 (stdlib), pandas, pytest.

## Global Constraints

- Python 3.11+.
- `load_data.py` MUST live in the repository root, run as `python load_data.py` with no CLI arguments and no `-m` module execution, and create a `.db` file in the repo root.
- Database is SQLite. Cell counts stored **long-format** (one row per sample×population), never wide.
- Subject-level attributes (`condition`, `age`, `sex`, `treatment`, `response`) live on `subjects`; sample-varying attributes (`sample_type`, `time_from_treatment_start`) live on `samples`. `response` is nullable.
- The five populations, in this canonical order, are: `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`.
- Loading is idempotent: rebuilding drops any existing DB and recreates it from scratch.
- Dependencies are pinned in `requirements.txt`.

---

## File Structure

- `requirements.txt` — pinned dependencies (this part: pandas, pytest).
- `conftest.py` — empty; ensures repo root is on `sys.path` so tests can `import load_data` / `import db.connection`.
- `db/__init__.py` — empty; makes `db` an importable package.
- `db/schema.sql` — DDL for the five tables + indexes.
- `db/connection.py` — `get_connection()` helper (foreign keys on, `Row` factory).
- `load_data.py` — root loader: `create_schema`, `read_csv`, `build_tables`, `build_database`, `main`.
- `tests/test_schema.py` — schema applies; FK constraints enforced.
- `tests/test_build_tables.py` — CSV → per-table DataFrames, long-format reshape.
- `tests/test_load_data.py` — end-to-end build on a fixture; idempotency; response NULL.

---

### Task 1: Project scaffolding + schema + connection helper

**Files:**
- Create: `requirements.txt`
- Create: `conftest.py`
- Create: `db/__init__.py`
- Create: `db/schema.sql`
- Create: `db/connection.py`
- Test: `tests/test_schema.py`

**Interfaces:**
- Consumes: nothing (first task).
- Produces:
  - `db.connection.get_connection(db_path=DEFAULT_DB_PATH) -> sqlite3.Connection` — foreign keys ON, `row_factory = sqlite3.Row`. `DEFAULT_DB_PATH` is repo-root `cell_counts.db`.
  - `db/schema.sql` — creates tables `projects`, `subjects`, `samples`, `populations`, `cell_counts` and their indexes.

- [ ] **Step 1: Create `requirements.txt`**

```text
pandas==2.2.2
pytest==8.2.2
```

- [ ] **Step 2: Create `conftest.py` (empty) and `db/__init__.py` (empty)**

Both files are empty. `conftest.py` at the repo root makes pytest treat the root as rootdir and put it on `sys.path`; `db/__init__.py` makes `db` a package.

```bash
: > conftest.py
: > db/__init__.py
```

- [ ] **Step 3: Create `db/schema.sql`**

```sql
-- Normalized, long-format schema for immune cell population data.
-- Foreign keys are enabled per-connection in db/connection.py (a PRAGMA in a
-- schema file does not persist), so this DDL relies on that being set.

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT,
    treatment  TEXT,
    response   TEXT              -- nullable: healthy / 'none' treatment have no response
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

CREATE TABLE populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id     TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES populations(population_id),
    count         INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project   ON subjects(project_id);
CREATE INDEX idx_subjects_condition ON subjects(condition);
CREATE INDEX idx_subjects_treatment ON subjects(treatment);
CREATE INDEX idx_subjects_response  ON subjects(response);
CREATE INDEX idx_samples_subject    ON samples(subject_id);
CREATE INDEX idx_samples_type       ON samples(sample_type);
CREATE INDEX idx_samples_time       ON samples(time_from_treatment_start);
CREATE INDEX idx_cell_counts_pop    ON cell_counts(population_id, sample_id);
```

- [ ] **Step 4: Create `db/connection.py`**

```python
"""SQLite connection helper."""
import sqlite3
from pathlib import Path

# Repo root is the parent of the db/ package directory.
DEFAULT_DB_PATH = Path(__file__).resolve().parent.parent / "cell_counts.db"


def get_connection(db_path=DEFAULT_DB_PATH):
    """Return a SQLite connection with foreign keys enabled and Row access."""
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
    return conn
```

- [ ] **Step 5: Write the failing test `tests/test_schema.py`**

`create_schema` is defined in `load_data.py` in Task 2, but the test file can be written now against the intended interface. To keep this task self-contained, apply the schema by reading the file directly here.

```python
import sqlite3
from pathlib import Path

import pytest

from db.connection import get_connection

SCHEMA_PATH = Path(__file__).resolve().parent.parent / "db" / "schema.sql"


def _fresh_db(tmp_path):
    conn = get_connection(tmp_path / "test.db")
    conn.executescript(SCHEMA_PATH.read_text())
    return conn


def test_schema_creates_all_tables(tmp_path):
    conn = _fresh_db(tmp_path)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"projects", "subjects", "samples", "populations",
            "cell_counts"} <= tables
    conn.close()


def test_foreign_keys_enforced(tmp_path):
    conn = _fresh_db(tmp_path)
    # A sample referencing a non-existent subject must fail.
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO samples VALUES ('smp1', 'missing', 'PBMC', 0)")
        conn.commit()
    conn.close()
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `pytest tests/test_schema.py -v`
Expected: PASS (2 tests). Both tables exist and the FK insert raises `IntegrityError`.

- [ ] **Step 7: Commit**

```bash
git add requirements.txt conftest.py db/ tests/test_schema.py
git commit -m "feat(part1): add SQLite schema and connection helper"
```

---

### Task 2: CSV read + reshape to per-table DataFrames

**Files:**
- Create: `load_data.py`
- Test: `tests/test_build_tables.py`

**Interfaces:**
- Consumes: `db.connection.get_connection` (from Task 1); `db/schema.sql`.
- Produces:
  - `load_data.POPULATION_COLUMNS: list[str]` — the five canonical names.
  - `load_data.create_schema(conn, schema_path=SCHEMA_PATH) -> None`.
  - `load_data.read_csv(csv_path=CSV_PATH) -> pd.DataFrame` — blank `response` coerced to `None`.
  - `load_data.build_tables(df) -> tuple[DataFrame, ...]` returning `(projects, subjects, samples, populations, cell_counts)` frames whose columns exactly match the schema.

- [ ] **Step 1: Create `load_data.py` with reshape logic (no `main` yet)**

```python
"""Part 1: build the SQLite schema and load cell-count.csv.

Run directly:
    python load_data.py
Creates cell_counts.db in the repository root.
"""
from pathlib import Path

import pandas as pd

from db.connection import get_connection

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "db" / "schema.sql"

POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def create_schema(conn, schema_path=SCHEMA_PATH):
    """Apply the DDL in schema_path to an open connection."""
    conn.executescript(Path(schema_path).read_text())


def read_csv(csv_path=CSV_PATH):
    """Read the raw CSV, coercing blank response cells to None (SQL NULL)."""
    df = pd.read_csv(csv_path)
    df["response"] = df["response"].where(df["response"].notna(), None)
    return df


def build_tables(df):
    """Reshape the raw dataframe into the five normalized tables."""
    projects = (
        df[["project"]]
        .drop_duplicates()
        .rename(columns={"project": "project_id"})
        .reset_index(drop=True)
    )

    subjects = (
        df[["subject", "project", "condition", "age", "sex",
            "treatment", "response"]]
        .drop_duplicates(subset="subject")
        .rename(columns={"subject": "subject_id", "project": "project_id"})
        .reset_index(drop=True)
    )

    samples = (
        df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
        .drop_duplicates(subset="sample")
        .rename(columns={"sample": "sample_id", "subject": "subject_id"})
        .reset_index(drop=True)
    )

    populations = pd.DataFrame({
        "population_id": range(1, len(POPULATION_COLUMNS) + 1),
        "name": POPULATION_COLUMNS,
    })

    name_to_id = dict(zip(populations["name"], populations["population_id"]))
    long = df.melt(
        id_vars=["sample"],
        value_vars=POPULATION_COLUMNS,
        var_name="name",
        value_name="count",
    )
    long["population_id"] = long["name"].map(name_to_id)
    cell_counts = (
        long[["sample", "population_id", "count"]]
        .rename(columns={"sample": "sample_id"})
        .reset_index(drop=True)
    )

    return projects, subjects, samples, populations, cell_counts
```

- [ ] **Step 2: Write the failing test `tests/test_build_tables.py`**

```python
import pandas as pd

from load_data import build_tables, POPULATION_COLUMNS


def _sample_df():
    # Two subjects; subject s0 has two samples, s1 has one.
    rows = [
        dict(project="prj1", subject="s0", condition="melanoma", age=50,
             sex="M", treatment="miraclib", response="yes",
             sample="smp0", sample_type="PBMC", time_from_treatment_start=0,
             b_cell=10, cd8_t_cell=20, cd4_t_cell=30, nk_cell=40, monocyte=100),
        dict(project="prj1", subject="s0", condition="melanoma", age=50,
             sex="M", treatment="miraclib", response="yes",
             sample="smp1", sample_type="PBMC", time_from_treatment_start=7,
             b_cell=1, cd8_t_cell=2, cd4_t_cell=3, nk_cell=4, monocyte=90),
        dict(project="prj2", subject="s1", condition="healthy", age=40,
             sex="F", treatment="none", response=None,
             sample="smp2", sample_type="WB", time_from_treatment_start=0,
             b_cell=5, cd8_t_cell=5, cd4_t_cell=5, nk_cell=5, monocyte=5),
    ]
    return pd.DataFrame(rows)


def test_table_row_counts():
    projects, subjects, samples, populations, cell_counts = build_tables(_sample_df())
    assert len(projects) == 2          # prj1, prj2
    assert len(subjects) == 2          # s0, s1
    assert len(samples) == 3           # smp0, smp1, smp2
    assert len(populations) == 5
    assert len(cell_counts) == 15      # 3 samples x 5 populations


def test_cell_counts_long_format_and_totals():
    _, _, _, populations, cell_counts = build_tables(_sample_df())
    # smp0 counts sum to 10+20+30+40+100 = 200
    smp0 = cell_counts[cell_counts["sample_id"] == "smp0"]
    assert len(smp0) == 5
    assert int(smp0["count"].sum()) == 200
    # population_id values are all valid ids from the lookup
    valid_ids = set(populations["population_id"])
    assert set(cell_counts["population_id"]) <= valid_ids


def test_subject_attributes_deduped():
    _, subjects, _, _, _ = build_tables(_sample_df())
    s0 = subjects[subjects["subject_id"] == "s0"]
    assert len(s0) == 1
    assert s0.iloc[0]["treatment"] == "miraclib"


def test_column_names_match_schema():
    projects, subjects, samples, populations, cell_counts = build_tables(_sample_df())
    assert list(projects.columns) == ["project_id"]
    assert list(subjects.columns) == [
        "subject_id", "project_id", "condition", "age", "sex",
        "treatment", "response"]
    assert list(samples.columns) == [
        "sample_id", "subject_id", "sample_type", "time_from_treatment_start"]
    assert list(populations.columns) == ["population_id", "name"]
    assert list(cell_counts.columns) == ["sample_id", "population_id", "count"]
```

- [ ] **Step 3: Run the tests to verify they pass**

Run: `pytest tests/test_build_tables.py -v`
Expected: PASS (4 tests).

- [ ] **Step 4: Commit**

```bash
git add load_data.py tests/test_build_tables.py
git commit -m "feat(part1): reshape CSV into normalized per-table frames"
```

---

### Task 3: End-to-end database build + `main` + real-data verification

**Files:**
- Modify: `load_data.py` (add `build_database`, `main`, `__main__` guard)
- Test: `tests/test_load_data.py`

**Interfaces:**
- Consumes: `create_schema`, `read_csv`, `build_tables`, `get_connection`.
- Produces:
  - `load_data.build_database(csv_path=CSV_PATH, db_path=DB_PATH, schema_path=SCHEMA_PATH) -> Path` — deletes any existing db, creates schema, loads all tables in FK-dependency order, commits, returns the db path.
  - `load_data.main() -> None` — builds the default DB and prints per-table row counts.

- [ ] **Step 1: Append `build_database`, `main`, and the entry point to `load_data.py`**

```python
def build_database(csv_path=CSV_PATH, db_path=DB_PATH, schema_path=SCHEMA_PATH):
    """Build the SQLite database from scratch. Idempotent."""
    db_path = Path(db_path)
    if db_path.exists():
        db_path.unlink()  # rebuild from scratch on every run

    conn = get_connection(db_path)
    try:
        create_schema(conn, schema_path)
        df = read_csv(csv_path)
        projects, subjects, samples, populations, cell_counts = build_tables(df)
        # Insert in FK-dependency order.
        projects.to_sql("projects", conn, if_exists="append", index=False)
        subjects.to_sql("subjects", conn, if_exists="append", index=False)
        samples.to_sql("samples", conn, if_exists="append", index=False)
        populations.to_sql("populations", conn, if_exists="append", index=False)
        cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)
        conn.commit()
    finally:
        conn.close()
    return db_path


def main():
    db_path = build_database()
    conn = get_connection(db_path)
    tables = ["projects", "subjects", "samples", "populations", "cell_counts"]
    counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
              for t in tables}
    conn.close()
    print(f"Built {db_path}")
    for table, n in counts.items():
        print(f"  {table}: {n} rows")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write the failing test `tests/test_load_data.py`**

```python
from db.connection import get_connection
from load_data import build_database

FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100
prj1,s0,melanoma,50,M,miraclib,yes,smp1,PBMC,7,1,2,3,4,90
prj2,s1,healthy,40,F,none,,smp2,WB,0,5,5,5,5,5
"""


def _build(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_row_counts(tmp_path):
    conn = _build(tmp_path)
    assert conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 2
    assert conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0] == 2
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 3
    assert conn.execute("SELECT COUNT(*) FROM populations").fetchone()[0] == 5
    assert conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0] == 15
    conn.close()


def test_blank_response_is_null(tmp_path):
    conn = _build(tmp_path)
    resp = conn.execute(
        "SELECT response FROM subjects WHERE subject_id = 's1'").fetchone()[0]
    assert resp is None
    conn.close()


def test_counts_join_back_to_populations(tmp_path):
    conn = _build(tmp_path)
    total = conn.execute(
        "SELECT SUM(count) FROM cell_counts WHERE sample_id = 'smp0'"
    ).fetchone()[0]
    assert total == 200
    b_cell = conn.execute(
        "SELECT count FROM cell_counts cc "
        "JOIN populations p ON p.population_id = cc.population_id "
        "WHERE cc.sample_id = 'smp0' AND p.name = 'b_cell'"
    ).fetchone()[0]
    assert b_cell == 10
    conn.close()


def test_idempotent_rebuild(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    build_database(csv_path=csv, db_path=db)  # second run must not duplicate
    conn = get_connection(db)
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 3
    conn.close()
```

- [ ] **Step 3: Run the tests to verify they pass**

Run: `pytest tests/test_load_data.py -v`
Expected: PASS (4 tests).

- [ ] **Step 4: Run the full test suite**

Run: `pytest -v`
Expected: PASS (all tests from Tasks 1–3).

- [ ] **Step 5: Build the real database and verify counts**

Run: `python load_data.py`
Expected output (a `cell_counts.db` file appears in the repo root):

```
Built /Users/zoyagoel/teiko/cell_counts.db
  projects: 3 rows
  subjects: 3500 rows
  samples: 10500 rows
  populations: 5 rows
  cell_counts: 52500 rows
```

If the counts differ, stop and investigate before committing.

- [ ] **Step 6: Commit**

```bash
git add load_data.py tests/test_load_data.py
git commit -m "feat(part1): end-to-end database build with load_data.py"
```

---

## Self-Review

**Spec coverage (against the overarching design §3, §7 and Part 1 requirements):**
- SQLite schema with five long-format tables → Task 1 (`db/schema.sql`). ✓
- `load_data.py` in root, no args, no `-m`, creates `.db` in root → Task 3 (`main` + `__main__`, `DB_PATH` at root). ✓
- Loads all rows from `cell-count.csv` → Task 3 (verified: 10,500 samples / 52,500 counts). ✓
- Subject vs sample attribute placement; nullable `response` → Task 2 (`build_tables`), Task 3 (`test_blank_response_is_null`). ✓
- Idempotent rebuild → Task 3 (`build_database` unlinks; `test_idempotent_rebuild`). ✓
- Indexes on FKs and filter columns → Task 1 (`schema.sql`). ✓
- Connection helper with FK enforcement → Task 1 (`get_connection`; `test_foreign_keys_enforced`). ✓

**Placeholder scan:** No TBD/TODO; every code and test step contains complete content. ✓

**Type consistency:** `create_schema`, `read_csv`, `build_tables`, `build_database`, `get_connection`, `POPULATION_COLUMNS`, and the `(projects, subjects, samples, populations, cell_counts)` tuple are named identically across the interface blocks and code. Column lists in `build_tables` match `schema.sql` (checked by `test_column_names_match_schema`). ✓

**Note for later parts:** `requirements.txt` here pins only `pandas` + `pytest`; the Part 3 and dashboard plans will append `scipy`, `statsmodels`, `plotly`, `streamlit`. The README already lists the full intended set.
