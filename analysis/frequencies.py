"""Part 2: per-sample relative frequency of each immune cell population."""
import pandas as pd

# Relative frequency computed in SQL: the per-sample total is a window sum over
# the sample's rows, so each population's percentage is count / total * 100.
FREQUENCY_QUERY = """
SELECT
    s.sample_id                                    AS sample,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
    p.name                                         AS population,
    cc.count                                       AS count,
    100.0 * cc.count
        / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts cc
JOIN samples s     ON s.sample_id = cc.sample_id
JOIN populations p ON p.population_id = cc.population_id
ORDER BY s.sample_id, p.population_id
"""


def sample_frequencies(conn):
    """Return the per-sample relative-frequency table.

    Columns: sample, total_count, population, count, percentage.
    One row per (sample, population); percentages within a sample sum to 100.
    """
    return pd.read_sql_query(FREQUENCY_QUERY, conn)
