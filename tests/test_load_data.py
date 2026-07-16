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
