import sqlite3
from pathlib import Path

import pytest

from db.connection import get_connection

SCHEMA_PATH = Path(__file__).resolve().parent.parent / "db" / "schema.sql"


def _fresh_db(tmp_path):
    conn = get_connection(tmp_path / "test.db")
    conn.executescript(SCHEMA_PATH.read_text())
    return conn


def test_schema_creates_all_tables(tmp_path):
    conn = _fresh_db(tmp_path)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"projects", "subjects", "samples", "populations",
            "cell_counts"} <= tables
    conn.close()


def test_foreign_keys_enforced(tmp_path):
    conn = _fresh_db(tmp_path)
    # A sample referencing a non-existent subject must fail.
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO samples VALUES ('smp1', 'missing', 'PBMC', 0)")
        conn.commit()
    conn.close()
