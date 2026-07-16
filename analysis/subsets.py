"""Part 4: baseline data-subset analysis.

Identifies the baseline (time_from_treatment_start = 0) melanoma / miraclib /
PBMC cohort and breaks it down by project, response, and sex.
"""
import pandas as pd

_BASELINE_QUERY = """
SELECT sa.sample_id               AS sample,
       su.subject_id              AS subject_id,
       su.project_id              AS project,
       su.response               AS response,
       su.sex                    AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
WHERE su.condition = :condition
  AND su.treatment = :treatment
  AND sa.sample_type = :sample_type
  AND sa.time_from_treatment_start = 0
ORDER BY sa.sample_id
"""


def baseline_subset(conn, condition="melanoma", treatment="miraclib",
                    sample_type="PBMC"):
    """All baseline (t=0) samples for the given cohort filter.

    Columns: sample, subject_id, project, response, sex.
    """
    params = {"condition": condition, "treatment": treatment,
              "sample_type": sample_type}
    return pd.read_sql_query(_BASELINE_QUERY, conn, params=params)


def samples_per_project(subset):
    """Number of samples per project in the subset."""
    return (
        subset.groupby("project").size()
        .reset_index(name="n_samples")
        .sort_values("project")
        .reset_index(drop=True)
    )


def subjects_by_response(subset):
    """Number of distinct subjects per response value in the subset."""
    return (
        subset.drop_duplicates("subject_id")
        .groupby("response").size()
        .reset_index(name="n_subjects")
        .sort_values("response")
        .reset_index(drop=True)
    )


def subjects_by_sex(subset):
    """Number of distinct subjects per sex in the subset."""
    return (
        subset.drop_duplicates("subject_id")
        .groupby("sex").size()
        .reset_index(name="n_subjects")
        .sort_values("sex")
        .reset_index(drop=True)
    )
