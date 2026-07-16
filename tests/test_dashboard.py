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


def test_app_runs_and_renders_response_finding():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file("dashboard/app.py").run(timeout=120)
    assert not at.exception
    # The Response tab computes and renders a plain-language finding via st.info.
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
