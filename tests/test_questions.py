from db.connection import get_connection
from load_data import build_database
from analysis.questions import avg_bcells_melanoma_male_responders_at_baseline

# 3 qualifying rows (melanoma, M, yes, t=0) spanning different sample_type /
# treatment, with b_cell 100, 100, 101 -> mean 100.3333 -> 100.33 (rounding).
# Plus rows excluded by each dimension.
FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,100,10,10,10,10
prj1,s1,melanoma,55,M,phauximab,yes,smp1,WB,0,100,10,10,10,10
prj1,s2,melanoma,60,M,miraclib,yes,smp2,PBMC,0,101,10,10,10,10
prj1,s3,melanoma,50,M,miraclib,yes,smp3,PBMC,7,9000,10,10,10,10
prj1,s4,melanoma,50,F,miraclib,yes,smp4,PBMC,0,9000,10,10,10,10
prj1,s5,melanoma,50,M,miraclib,no,smp5,PBMC,0,9000,10,10,10,10
prj1,s6,carcinoma,50,M,miraclib,yes,smp6,PBMC,0,9000,10,10,10,10
"""


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_average_b_cells_rounded_two_decimals(tmp_path):
    avg = avg_bcells_melanoma_male_responders_at_baseline(_conn(tmp_path))
    assert avg == 100.33  # mean(100, 100, 101) = 100.3333... -> 100.33


def test_returns_float(tmp_path):
    avg = avg_bcells_melanoma_male_responders_at_baseline(_conn(tmp_path))
    assert isinstance(avg, float)
