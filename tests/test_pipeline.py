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
