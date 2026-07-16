import numpy as np
import pandas as pd

from db.connection import get_connection
from load_data import build_database
from analysis.comparison import compare_responders, POPULATIONS


def _fixture_csv(seed=0):
    """~24 melanoma+miraclib+PBMC subjects (12 responders/12 non), 3 timepoints,
    with a planted cd4_t_cell elevation in responders; other populations flat."""
    rng = np.random.default_rng(seed)
    header = ("project,subject,condition,age,sex,treatment,response,sample,"
              "sample_type,time_from_treatment_start,b_cell,cd8_t_cell,"
              "cd4_t_cell,nk_cell,monocyte")
    lines = [header]
    sample_no = 0
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(12):
            sid = f"{resp}{s:02d}"
            for t in (0, 7, 14):
                cd4 = 300 + 40 * is_r + rng.integers(-10, 10)
                b = 100 + rng.integers(-10, 10)
                cd8 = 250 + rng.integers(-10, 10)
                nk = 150 + rng.integers(-10, 10)
                mono = 200 + rng.integers(-10, 10)
                lines.append(
                    f"prj1,{sid},melanoma,50,M,miraclib,{resp},smp{sample_no:04d},"
                    f"PBMC,{t},{b},{cd8},{cd4},{nk},{mono}")
                sample_no += 1
    return "\n".join(lines) + "\n"


def _conn(tmp_path):
    csv = tmp_path / "cell-count.csv"
    csv.write_text(_fixture_csv())
    db = tmp_path / "test.db"
    build_database(csv_path=csv, db_path=db)
    return get_connection(db)


def test_result_has_one_row_per_population_in_order(tmp_path):
    res = compare_responders(_conn(tmp_path))
    assert list(res["population"]) == POPULATIONS
    assert len(res) == 5


def test_result_has_all_expected_columns(tmp_path):
    res = compare_responders(_conn(tmp_path))
    for col in ["population", "n_responders", "n_nonresponders", "mean_resp",
                "mean_nonresp", "mean_diff", "boot_ci_low", "boot_ci_high",
                "coef", "ci_low", "ci_high", "mixed_p", "q_mixed", "group_var",
                "mw_U", "mw_p", "q_mw", "rank_biserial", "median_resp",
                "median_nonresp", "significant_primary", "concordant"]:
        assert col in res.columns


def test_fdr_and_flags_are_well_formed(tmp_path):
    res = compare_responders(_conn(tmp_path))
    # q-values are valid probabilities and never smaller than their raw p
    assert (res["q_mixed"] >= res["mixed_p"] - 1e-9).all()
    assert (res["q_mixed"] <= 1.0).all() and (res["q_mixed"] >= 0.0).all()
    # flags are boolean and concordant implies primary
    assert res["significant_primary"].dtype == bool
    assert res["concordant"].dtype == bool
    assert (~res["concordant"] | res["significant_primary"]).all()
    assert res["n_responders"].iloc[0] == 12
    assert res["n_nonresponders"].iloc[0] == 12


def test_planted_cd4_effect_is_detected(tmp_path):
    res = compare_responders(_conn(tmp_path)).set_index("population")
    # cd4 was planted higher in responders
    assert res.loc["cd4_t_cell", "mean_diff"] > 0
    assert res.loc["cd4_t_cell", "significant_primary"]
    # guards the group_means pivot label mapping (yes -> mean_resp)
    assert res.loc["cd4_t_cell", "mean_resp"] > res.loc["cd4_t_cell", "mean_nonresp"]
