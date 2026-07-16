"""Part 4: the specific posed question.

Average number of B cells for melanoma male responders at baseline (t=0),
across all sample types and all treatments.
"""
import pandas as pd

_BCELLS_QUERY = """
SELECT cc.count AS b_cell
FROM subjects su
JOIN samples sa      ON sa.subject_id = su.subject_id
JOIN cell_counts cc  ON cc.sample_id = sa.sample_id
JOIN populations p   ON p.population_id = cc.population_id
WHERE su.condition = 'melanoma'
  AND su.sex = 'M'
  AND su.response = 'yes'
  AND sa.time_from_treatment_start = 0
  AND p.name = 'b_cell'
"""


def avg_bcells_melanoma_male_responders_at_baseline(conn):
    """Mean B-cell count for melanoma male responders at t=0 (2 decimals).

    All sample types and treatments are included, per the question.
    """
    df = pd.read_sql_query(_BCELLS_QUERY, conn)
    return round(float(df["b_cell"].mean()), 2)
