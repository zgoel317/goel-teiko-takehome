import pandas as pd

from db.connection import get_connection
from load_data import build_database
from analysis.comparison import cohort_frequencies, subject_means, POPULATIONS

# 2 melanoma+miraclib+PBMC subjects (1 responder, 1 non-responder), each 2 timepoints,
# plus rows that must be EXCLUDED: wrong condition, wrong treatment, wrong sample_type,
# and a null-response (healthy) subject.
FIXTURE_CSV = """\
project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
prj1,s0,melanoma,50,M,miraclib,yes,smp0,PBMC,0,10,20,30,40,100
prj1,s0,melanoma,50,M,miraclib,yes,smp1,PBMC,7,20,20,30,40,90
prj1,s1,melanoma,60,F,miraclib,no,smp2,PBMC,0,5,20,30,40,105
prj1,s1,melanoma,60,F,miraclib,no,smp3,PBMC,7,15,20,30,40,95
prj1,s2,melanoma,55,M,miraclib,yes,smp4,WB,0,10,20,30,40,100
prj1,s3,carcinoma,55,M,miraclib,yes,smp5,PBMC,0,10,20,30,40,100
prj1,s4,melanoma,55,M,phauximab,no,smp6,PBMC,0,10,20,30,40,100
prj2,s5,healthy,40,F,none,,smp7,PBMC,0,5,5,5,5,5
"""


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(FIXTURE_CSV)
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_cohort_filters_to_melanoma_miraclib_pbmc_with_response(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    # Only s0 and s1 qualify: 2 subjects x 2 timepoints x 5 populations = 20 rows
    assert set(cohort["subject_id"]) == {"s0", "s1"}
    assert len(cohort) == 20
    assert set(cohort["response"]) == {"yes", "no"}
    # Excluded units must not appear
    assert "s2" not in set(cohort["subject_id"])  # WB
    assert "s3" not in set(cohort["subject_id"])  # carcinoma
    assert "s4" not in set(cohort["subject_id"])  # phauximab
    assert "s5" not in set(cohort["subject_id"])  # healthy / null response


def test_cohort_has_expected_columns(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    for col in ["sample", "total_count", "population", "count", "percentage",
                "subject_id", "time", "response", "sex"]:
        assert col in cohort.columns


def test_timepoint_filter(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path), timepoint=0)
    assert set(cohort["time"]) == {0}
    assert len(cohort) == 10  # 2 subjects x 1 timepoint x 5 populations


def test_subject_means_collapses_to_one_row_per_subject_population(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    sm = subject_means(cohort)
    assert list(sm.columns) == ["subject_id", "response", "population", "percentage"]
    assert len(sm) == 10  # 2 subjects x 5 populations
    # s0 b_cell: sample smp0 total=200 ->5%, smp1 total=200 ->10%; mean=7.5
    s0_b = sm[(sm.subject_id == "s0") & (sm.population == "b_cell")].iloc[0]
    assert abs(s0_b["percentage"] - 7.5) < 1e-9


def test_subject_means_preserves_canonical_population_order(tmp_path):
    cohort = cohort_frequencies(_conn(tmp_path))
    sm = subject_means(cohort)
    s0 = sm[sm["subject_id"] == "s0"]
    assert list(s0["population"]) == POPULATIONS
