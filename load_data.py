"""Part 1: build the SQLite schema and load cell-count.csv.

Run directly:
    python load_data.py
Creates cell_counts.db in the repository root.
"""
from pathlib import Path

import pandas as pd

from db.connection import get_connection

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "db" / "schema.sql"

POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def create_schema(conn, schema_path=SCHEMA_PATH):
    """Apply the DDL in schema_path to an open connection."""
    conn.executescript(Path(schema_path).read_text())


def read_csv(csv_path=CSV_PATH):
    """Read the raw CSV, coercing blank response cells to None (SQL NULL)."""
    df = pd.read_csv(csv_path)
    df["response"] = df["response"].where(df["response"].notna(), None)
    return df


def build_tables(df):
    """Reshape the raw dataframe into the five normalized tables."""
    projects = (
        df[["project"]]
        .drop_duplicates()
        .rename(columns={"project": "project_id"})
        .reset_index(drop=True)
    )

    subjects = (
        df[["subject", "project", "condition", "age", "sex",
            "treatment", "response"]]
        .drop_duplicates(subset="subject")
        .rename(columns={"subject": "subject_id", "project": "project_id"})
        .reset_index(drop=True)
    )

    samples = (
        df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
        .drop_duplicates(subset="sample")
        .rename(columns={"sample": "sample_id", "subject": "subject_id"})
        .reset_index(drop=True)
    )

    populations = pd.DataFrame({
        "population_id": range(1, len(POPULATION_COLUMNS) + 1),
        "name": POPULATION_COLUMNS,
    })

    name_to_id = dict(zip(populations["name"], populations["population_id"]))
    long = df.melt(
        id_vars=["sample"],
        value_vars=POPULATION_COLUMNS,
        var_name="name",
        value_name="count",
    )
    long["population_id"] = long["name"].map(name_to_id)
    cell_counts = (
        long[["sample", "population_id", "count"]]
        .rename(columns={"sample": "sample_id"})
        .reset_index(drop=True)
    )

    return projects, subjects, samples, populations, cell_counts


def build_database(csv_path=CSV_PATH, db_path=DB_PATH, schema_path=SCHEMA_PATH):
    """Build the SQLite database from scratch. Idempotent."""
    db_path = Path(db_path)
    if db_path.exists():
        db_path.unlink()  # rebuild from scratch on every run

    conn = get_connection(db_path)
    try:
        create_schema(conn, schema_path)
        df = read_csv(csv_path)
        projects, subjects, samples, populations, cell_counts = build_tables(df)
        # Insert in FK-dependency order.
        projects.to_sql("projects", conn, if_exists="append", index=False)
        subjects.to_sql("subjects", conn, if_exists="append", index=False)
        samples.to_sql("samples", conn, if_exists="append", index=False)
        populations.to_sql("populations", conn, if_exists="append", index=False)
        cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)
        conn.commit()
    finally:
        conn.close()
    return db_path


def main():
    db_path = build_database()
    conn = get_connection(db_path)
    tables = ["projects", "subjects", "samples", "populations", "cell_counts"]
    counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
              for t in tables}
    conn.close()
    print(f"Built {db_path}")
    for table, n in counts.items():
        print(f"  {table}: {n} rows")


if __name__ == "__main__":
    main()
