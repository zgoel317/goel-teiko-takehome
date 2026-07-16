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
