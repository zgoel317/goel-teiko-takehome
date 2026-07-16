import numpy as np
import pandas as pd

from analysis.comparison import (
    benjamini_hochberg,
    mannwhitney_by_population,
    mixedlm_by_population,
    bootstrap_diff_ci,
)


def _synthetic_cohort(effect=5.0, seed=0, n_per=40):
    """Two populations: 'A' has a planted responder effect, 'B' has none.
    30/40 subjects per group, 3 timepoints each, subject-level random intercept."""
    rng = np.random.default_rng(seed)
    rows = []
    for resp, is_r in [("yes", 1), ("no", 0)]:
        for s in range(n_per):
            sid = f"{resp}_{s}"
            intercept = rng.normal(0, 1)  # subject random effect
            for t in (0, 7, 14):
                rows.append(dict(subject_id=sid, response=resp, time=t,
                                 population="A",
                                 percentage=20 + effect * is_r + intercept + rng.normal(0, 1)))
                rows.append(dict(subject_id=sid, response=resp, time=t,
                                 population="B",
                                 percentage=20 + intercept + rng.normal(0, 1)))
    return pd.DataFrame(rows)


def test_benjamini_hochberg_known_values():
    q = benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
    # BH: q_i = p_i * n / rank, then cumulative-min from the top.
    # All equal to 0.05 here since p_i * 5 / i == 0.05 for every i.
    assert np.allclose(q, [0.05, 0.05, 0.05, 0.05, 0.05])
    # q-values are clipped to [0,1] and monotone non-decreasing in p
    q2 = benjamini_hochberg([0.9, 0.001])
    assert np.all(q2 <= 1.0)
    assert q2[1] <= q2[0]


def test_mannwhitney_detects_planted_effect():
    res = mannwhitney_by_population(_synthetic_cohort(effect=5.0), ["A", "B"])
    a = res[res.population == "A"].iloc[0]
    b = res[res.population == "B"].iloc[0]
    assert a["mw_p"] < 0.05          # A differs
    assert a["rank_biserial"] > 0    # responders higher
    assert b["mw_p"] > 0.05          # B does not
    assert a["n_responders"] == 40 and a["n_nonresponders"] == 40


def test_mixedlm_recovers_effect_and_reports_group_var():
    res = mixedlm_by_population(_synthetic_cohort(effect=5.0), ["A", "B"])
    a = res[res.population == "A"].iloc[0]
    assert abs(a["coef"] - 5.0) < 1.5     # coefficient near planted effect
    assert a["mixed_p"] < 0.05
    assert a["ci_low"] < a["coef"] < a["ci_high"]
    assert "group_var" in res.columns and a["group_var"] >= 0


def test_bootstrap_ci_brackets_diff_and_is_reproducible():
    coh = _synthetic_cohort(effect=5.0)
    r1 = bootstrap_diff_ci(coh, ["A", "B"], n_boot=500, seed=0)
    r2 = bootstrap_diff_ci(coh, ["A", "B"], n_boot=500, seed=0)
    a = r1[r1.population == "A"].iloc[0]
    assert a["boot_ci_low"] < a["mean_diff"] < a["boot_ci_high"]
    assert a["boot_ci_low"] > 0            # whole CI above 0 for a strong effect
    pd.testing.assert_frame_equal(r1, r2)  # seeded => identical
