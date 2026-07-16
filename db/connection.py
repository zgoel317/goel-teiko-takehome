"""SQLite connection helper."""
import sqlite3
from pathlib import Path

# Repo root is the parent of the db/ package directory.
DEFAULT_DB_PATH = Path(__file__).resolve().parent.parent / "cell_counts.db"


def get_connection(db_path=DEFAULT_DB_PATH, check_same_thread=True):
    """Return a SQLite connection with foreign keys enabled and Row access.

    check_same_thread=False is used by the Streamlit dashboard, which caches one
    read-only connection but reruns on different threads.
    """
    conn = sqlite3.connect(str(db_path), check_same_thread=check_same_thread)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
    return conn
