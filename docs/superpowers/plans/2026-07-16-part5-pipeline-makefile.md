# Part 5: Pipeline + Makefile Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One command (`make pipeline`) that initializes the database, loads the data (Part 1), and regenerates every required output table and plot (Parts 2–4) into `outputs/`; plus `make setup` and `make dashboard` targets.

**Architecture:** `pipeline.py` is a thin *persist* layer: it rebuilds the DB via Part 1's `build_database`, then calls the same pure `analysis/` functions the dashboard uses and writes their results to `outputs/`. The `Makefile` exposes the three required targets. No analysis logic lives here — only orchestration and file writing.

**Tech Stack:** Python 3.11+, pandas, plotly (HTML export), GNU Make; `streamlit` added to deps for the (next) dashboard.

## Global Constraints

- `make setup` installs all dependencies from `requirements.txt`.
- `make pipeline` runs end-to-end with no manual steps: initialize DB, load CSV (Part 1), generate all Part 2–4 output tables and plots into `outputs/`. Idempotent (safe to re-run; rebuilds from scratch).
- `make dashboard` starts the Streamlit server (`dashboard/app.py`, created in the next spec).
- The pipeline reuses the existing pure functions — no re-implementation of any analysis.
- Use `python3` in the Makefile (present both locally and in Codespaces).
- `outputs/` is a committed deliverable (the submission requires generated output files) — do NOT gitignore it.
- The plot is exported as a self-contained interactive **HTML** file (no extra native/`kaleido` dependency).

---

## Output artifacts (`outputs/`)

| File | Source | Part |
|------|--------|------|
| `frequencies.csv` | `sample_frequencies(conn)` | 2 |
| `responder_comparison.csv` | `compare_responders(conn)` | 3 |
| `responder_boxplots.html` | `responder_boxplots(cohort_frequencies(conn))` | 3 |
| `baseline_subset.csv` | `baseline_subset(conn)` | 4 |
| `baseline_samples_per_project.csv` | `samples_per_project(subset)` | 4 |
| `baseline_subjects_by_response.csv` | `subjects_by_response(subset)` | 4 |
| `baseline_subjects_by_sex.csv` | `subjects_by_sex(subset)` | 4 |
| `final_answer.txt` | `avg_bcells_melanoma_male_responders_at_baseline(conn)` | 4 |

---

## File Structure

- `pipeline.py` — root orchestration script (`run_pipeline`, `main`, `__main__`).
- `Makefile` — `setup`, `pipeline`, `dashboard` targets.
- `requirements.txt` — add `streamlit` (modify).
- `tests/test_pipeline.py` — runs `run_pipeline` against a fixture DB into a temp dir and asserts all artifacts are produced with correct content.

---

### Task 1: `pipeline.py`

**Files:**
- Create: `pipeline.py`
- Test: `tests/test_pipeline.py`

**Interfaces:**
- Consumes: `load_data.build_database`; `db.connection.get_connection`; `analysis.frequencies.sample_frequencies`; `analysis.comparison.compare_responders`, `cohort_frequencies`, `responder_boxplots`; `analysis.subsets.baseline_subset`, `samples_per_project`, `subjects_by_response`, `subjects_by_sex`; `analysis.questions.avg_bcells_melanoma_male_responders_at_baseline`.
- Produces: `pipeline.run_pipeline(csv_path=CSV_PATH, db_path=DB_PATH, outdir=OUTPUTS_DIR) -> dict` — builds the DB, writes all eight artifacts, returns a dict mapping artifact name → written path. `pipeline.main()` calls it with defaults and prints a summary.

- [ ] **Step 1: Write the failing test `tests/test_pipeline.py`**

```python
from pathlib import Path

import pandas as pd

from pipeline import run_pipeline

# Minimal but valid cohort: enough melanoma+miraclib+PBMC subjects with both
# responses and 3 timepoints for compare_responders to fit, plus a melanoma
# male responder at t=0 for the final answer.
def _fixture_csv():
    header = ("project,subject,condition,age,sex,treatment,response,sample,"
              "sample_type,time_from_treatment_start,b_cell,cd8_t_cell,"
              "cd4_t_cell,nk_cell,monocyte")
    lines = [header]
    n = 0
    for resp in ("yes", "no"):
        for s in range(6):
            sid = f"{resp}{s}"
            sex = "M" if s % 2 == 0 else "F"
            for t in (0, 7, 14):
                lines.append(
                    f"prj1,{sid},melanoma,50,{sex},miraclib,{resp},smp{n:03d},"
                    f"PBMC,{t},{100 + s},200,300,150,200")
                n += 1
    return "\n".join(lines) + "\n"


def test_pipeline_writes_all_artifacts(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(_fixture_csv())
    db = tmp_path / "cell_counts.db"
    outdir = tmp_path / "outputs"

    written = run_pipeline(csv_path=csv, db_path=db, outdir=outdir)

    expected = {
        "frequencies", "responder_comparison", "responder_boxplots",
        "baseline_subset", "baseline_samples_per_project",
        "baseline_subjects_by_response", "baseline_subjects_by_sex",
        "final_answer",
    }
    assert set(written) == expected
    for name, path in written.items():
        assert Path(path).exists(), name
        assert Path(path).stat().st_size > 0, name


def test_pipeline_outputs_have_expected_content(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(_fixture_csv())
    db = tmp_path / "cell_counts.db"
    outdir = tmp_path / "outputs"
    run_pipeline(csv_path=csv, db_path=db, outdir=outdir)

    freq = pd.read_csv(outdir / "frequencies.csv")
    assert list(freq.columns) == [
        "sample", "total_count", "population", "count", "percentage"]

    comp = pd.read_csv(outdir / "responder_comparison.csv")
    assert len(comp) == 5  # one row per population
    assert "significant_primary" in comp.columns

    answer = (outdir / "final_answer.txt").read_text().strip()
    assert float(answer) >= 0  # parses as a number

    html = (outdir / "responder_boxplots.html").read_text()
    assert "plotly" in html.lower()  # a real plotly figure was written
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/zoyagoel/teiko && python3 -m pytest tests/test_pipeline.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pipeline'`.

- [ ] **Step 3: Create `pipeline.py`**

```python
"""Runs the full data pipeline: build the database (Part 1) and generate every
required output table and plot (Parts 2-4) into outputs/.

Run directly:
    python3 pipeline.py
"""
from pathlib import Path

from load_data import build_database, CSV_PATH, DB_PATH
from db.connection import get_connection
from analysis.frequencies import sample_frequencies
from analysis.comparison import (
    compare_responders,
    cohort_frequencies,
    responder_boxplots,
)
from analysis.subsets import (
    baseline_subset,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)
from analysis.questions import avg_bcells_melanoma_male_responders_at_baseline

ROOT = Path(__file__).resolve().parent
OUTPUTS_DIR = ROOT / "outputs"


def run_pipeline(csv_path=CSV_PATH, db_path=DB_PATH, outdir=OUTPUTS_DIR):
    """Build the DB and write all Part 2-4 artifacts. Returns {name: path}."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    build_database(csv_path=csv_path, db_path=db_path)  # Part 1
    conn = get_connection(db_path)
    written = {}
    try:
        # Part 2
        freq_path = outdir / "frequencies.csv"
        sample_frequencies(conn).to_csv(freq_path, index=False)
        written["frequencies"] = freq_path

        # Part 3
        comp_path = outdir / "responder_comparison.csv"
        compare_responders(conn).to_csv(comp_path, index=False)
        written["responder_comparison"] = comp_path

        box_path = outdir / "responder_boxplots.html"
        responder_boxplots(cohort_frequencies(conn)).write_html(box_path)
        written["responder_boxplots"] = box_path

        # Part 4
        subset = baseline_subset(conn)
        subset_path = outdir / "baseline_subset.csv"
        subset.to_csv(subset_path, index=False)
        written["baseline_subset"] = subset_path

        spp_path = outdir / "baseline_samples_per_project.csv"
        samples_per_project(subset).to_csv(spp_path, index=False)
        written["baseline_samples_per_project"] = spp_path

        sbr_path = outdir / "baseline_subjects_by_response.csv"
        subjects_by_response(subset).to_csv(sbr_path, index=False)
        written["baseline_subjects_by_response"] = sbr_path

        sbs_path = outdir / "baseline_subjects_by_sex.csv"
        subjects_by_sex(subset).to_csv(sbs_path, index=False)
        written["baseline_subjects_by_sex"] = sbs_path

        answer = avg_bcells_melanoma_male_responders_at_baseline(conn)
        answer_path = outdir / "final_answer.txt"
        answer_path.write_text(f"{answer}\n")
        written["final_answer"] = answer_path
    finally:
        conn.close()

    return {name: str(path) for name, path in written.items()}


def main():
    written = run_pipeline()
    print("Pipeline complete. Wrote:")
    for name, path in written.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_pipeline.py -v`
Expected: PASS (2 tests), pristine.

- [ ] **Step 5: Run the full suite**

Run: `python3 -m pytest -v`
Expected: PASS (all Parts 1–5 tests).

- [ ] **Step 6: Run the real pipeline and verify artifacts**

Run: `python3 pipeline.py`
Expected: prints "Pipeline complete." and the 8 artifact paths. Then confirm:

```bash
ls outputs/
cat outputs/final_answer.txt              # -> 10206.15
head -1 outputs/responder_comparison.csv  # header row with population,...,significant_primary,concordant
```

Expected: all 8 files present; `final_answer.txt` is `10206.15`; `responder_comparison.csv` has 5 data rows. If anything is missing or wrong, stop and investigate.

- [ ] **Step 7: Commit**

```bash
git add pipeline.py tests/test_pipeline.py outputs/
git commit -m "feat(part5): pipeline.py generating all Part 2-4 output artifacts"
```

---

### Task 2: Makefile + streamlit dependency

**Files:**
- Create: `Makefile`
- Modify: `requirements.txt` (add `streamlit`)

**Interfaces:**
- Consumes: `pipeline.py` (Task 1); `dashboard/app.py` (created in the next spec — the `dashboard` target references it).
- Produces: a `Makefile` with `.PHONY` targets `setup`, `pipeline`, `dashboard`.

- [ ] **Step 1: Add `streamlit` to `requirements.txt`**

Append this line (keep all existing lines):

```text
streamlit==1.38.0
```

- [ ] **Step 2: Create the `Makefile`**

Note: recipe lines must be indented with a TAB, not spaces.

```makefile
.PHONY: setup pipeline dashboard

setup:
	python3 -m pip install -r requirements.txt

pipeline:
	python3 pipeline.py

dashboard:
	streamlit run dashboard/app.py
```

- [ ] **Step 3: Verify `make setup` installs dependencies**

Run: `cd /Users/zoyagoel/teiko && make setup`
Expected: pip installs everything in `requirements.txt` without error (streamlit + the existing deps). It is fine if most are already satisfied.

- [ ] **Step 4: Verify `make pipeline` runs end-to-end**

Run: `make pipeline`
Expected: rebuilds `cell_counts.db` and prints the 8 artifact paths; `outputs/final_answer.txt` is `10206.15`. Confirm:

```bash
cat outputs/final_answer.txt
```

- [ ] **Step 5: Verify the `dashboard` target is wired (dry run)**

Run: `make -n dashboard`
Expected: prints `streamlit run dashboard/app.py` without executing it. (The dashboard app itself is built in the next spec; do not run `make dashboard` for real yet.)

- [ ] **Step 6: Commit**

```bash
git add Makefile requirements.txt
git commit -m "feat(part5): Makefile with setup/pipeline/dashboard targets"
```

---

## Self-Review

**Spec coverage (against submission requirements + overarching design §8):**
- `make setup` installs deps → Task 2 Makefile. ✓
- `make pipeline` initializes DB + loads (Part 1) + generates all Part 2–4 tables and plots → `run_pipeline` builds DB then writes 8 artifacts (Task 1); `pipeline` target (Task 2). ✓
- `make dashboard` starts the server → Task 2 target (app in next spec). ✓
- Idempotent pipeline → `build_database` unlinks + rebuilds; `outdir.mkdir(exist_ok=True)`; `to_csv` overwrites. ✓
- Outputs are committed deliverables → Task 1 Step 7 adds `outputs/`. ✓
- Reuses analysis functions (no re-implementation) → `run_pipeline` only calls existing functions. ✓

**Placeholder scan:** No TBD/TODO; all code, tests, and Makefile complete. ✓

**Type consistency:** `run_pipeline(csv_path, db_path, outdir) -> dict`; artifact names in the returned dict match the test's `expected` set and the artifact table. `build_database`, `CSV_PATH`, `DB_PATH` imported from `load_data` match Part 1's exports. ✓

**Note:** `make dashboard` references `dashboard/app.py`, which does not exist until the next spec — the target is defined now (verified via `make -n`) and becomes runnable once the dashboard is built. This is an intentional, documented forward reference, not a gap.
