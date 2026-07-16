from db.connection import get_connection
from load_data import build_database
from analysis.subsets import (
    baseline_subset,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)

# 3 rows that qualify (melanoma+miraclib+PBMC+t0) across 2 projects/both sexes/
# both responses, plus 4 rows that must be excluded (t=7, WB, phauximab, carcinoma).
FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,10,10,10,10
prj1,s1,melanoma,60,F,miraclib,no,smp1,PBMC,0,10,10,10,10,10
prj2,s2,melanoma,55,M,miraclib,yes,smp2,PBMC,0,10,10,10,10,10
prj1,s3,melanoma,50,M,miraclib,yes,smp3,PBMC,7,10,10,10,10,10
prj1,s4,melanoma,50,M,miraclib,yes,smp4,WB,0,10,10,10,10,10
prj1,s5,melanoma,50,M,phauximab,yes,smp5,PBMC,0,10,10,10,10,10
prj1,s6,carcinoma,50,M,miraclib,yes,smp6,PBMC,0,10,10,10,10,10
"""


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_baseline_subset_filters_correctly(tmp_path):
    sub = baseline_subset(_conn(tmp_path))
    assert set(sub["sample"]) == {"smp0", "smp1", "smp2"}
    assert list(sub.columns) == ["sample", "subject_id", "project", "response", "sex"]
    # excluded rows absent
    for excluded in ["smp3", "smp4", "smp5", "smp6"]:
        assert excluded not in set(sub["sample"])


def test_samples_per_project(tmp_path):
    sub = baseline_subset(_conn(tmp_path))
    spp = samples_per_project(sub).set_index("project")["n_samples"].to_dict()
    assert spp == {"prj1": 2, "prj2": 1}


def test_subjects_by_response(tmp_path):
    sub = baseline_subset(_conn(tmp_path))
    sbr = subjects_by_response(sub).set_index("response")["n_subjects"].to_dict()
    assert sbr == {"yes": 2, "no": 1}


def test_subjects_by_sex(tmp_path):
    sub = baseline_subset(_conn(tmp_path))
    sbs = subjects_by_sex(sub).set_index("sex")["n_subjects"].to_dict()
    assert sbs == {"M": 2, "F": 1}
