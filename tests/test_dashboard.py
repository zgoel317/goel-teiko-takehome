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


def test_app_runs_and_renders_part3_finding():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
    infos = " ".join(el.value for el in at.info)
    assert "responder" in infos.lower() or "cell population" in infos.lower()


def test_connection_usable_across_threads(tmp_path):
    import threading
    from db.connection import get_connection
    from load_data import build_database

    csv = tmp_path / "cell-count.csv"
    csv.write_text(
        "project,subject,condition,age,sex,treatment,response,sample,sample_type,"
        "time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte\n"
        "prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100\n"
    )
    db = tmp_path / "cell_counts.db"
    build_database(csv_path=csv, db_path=db)
    conn = get_connection(db_path=db, check_same_thread=False)
    errors = []

    def query():
        try:
            conn.execute("SELECT COUNT(*) FROM samples").fetchone()
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    t = threading.Thread(target=query)
    t.start()
    t.join()
    assert not errors


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
