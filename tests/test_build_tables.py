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


def test_read_csv_coerces_blank_response_to_none(tmp_path):
    from load_data import read_csv
    csv = tmp_path / "cell-count.csv"
    csv.write_text(
        "project,subject,condition,age,sex,treatment,response,sample,"
        "sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,"
        "nk_cell,monocyte\n"
        "prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100\n"
        "prj2,s1,healthy,40,F,none,,smp1,WB,0,5,5,5,5,5\n"
    )
    df = read_csv(csv)
    responses = dict(zip(df["sample"], df["response"]))
    assert responses["smp0"] == "yes"
    assert responses["smp1"] is None


def test_create_schema_applies_ddl(tmp_path):
    from db.connection import get_connection
    from load_data import create_schema
    conn = get_connection(tmp_path / "t.db")
    create_schema(conn)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"projects", "subjects", "samples", "populations",
            "cell_counts"} <= tables
    conn.close()
