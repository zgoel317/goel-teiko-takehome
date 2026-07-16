"""Cached data access for the dashboard: a self-building DB connection and
cached wrappers over the analysis functions."""
from pathlib import Path

import pandas as pd
import streamlit as st

from load_data import build_database, CSV_PATH, DB_PATH
from db.connection import get_connection as _open_connection
from analysis.frequencies import sample_frequencies
from analysis.comparison import (
    compare_responders,
    cohort_frequencies,
    responder_boxplots,
)
from analysis.subsets import baseline_subset
from analysis.questions import avg_bcells_melanoma_male_responders_at_baseline

_META_QUERY = """
SELECT sa.sample_id                 AS sample,
       sa.subject_id                AS subject_id,
       su.project_id                AS project,
       su.condition                 AS condition,
       su.treatment                 AS treatment,
       sa.sample_type               AS sample_type,
       sa.time_from_treatment_start AS time,
       su.response                  AS response,
       su.sex                       AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
"""


def ensure_database(db_path=DB_PATH, csv_path=CSV_PATH):
    """Build the SQLite DB from the CSV only if it does not already exist."""
    db_path = Path(db_path)
    if not db_path.exists():
        build_database(csv_path=csv_path, db_path=db_path)
    return db_path


@st.cache_resource
def get_connection():
    ensure_database()
    return _open_connection(check_same_thread=False)


@st.cache_data
def load_sample_metadata():
    return pd.read_sql_query(_META_QUERY, get_connection())


@st.cache_data
def load_frequencies():
    return sample_frequencies(get_connection())


@st.cache_data
def load_frequencies_annotated():
    return load_frequencies().merge(load_sample_metadata(), on="sample", how="left")


@st.cache_data
def load_dataset_summary():
    m = load_sample_metadata()
    return {
        "samples": int(m["sample"].nunique()),
        "subjects": int(m["subject_id"].nunique()),
        "projects": int(m["project"].nunique()),
    }


@st.cache_data
def load_comparison(condition, treatment, sample_type):
    return compare_responders(get_connection(), condition=condition,
                              treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_cohort(condition, treatment, sample_type):
    return cohort_frequencies(get_connection(), condition=condition,
                              treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_boxplot_figure(condition, treatment, sample_type):
    return responder_boxplots(load_cohort(condition, treatment, sample_type))


@st.cache_data
def load_baseline(condition, treatment, sample_type):
    return baseline_subset(get_connection(), condition=condition,
                           treatment=treatment, sample_type=sample_type)


@st.cache_data
def load_final_answer():
    return avg_bcells_melanoma_male_responders_at_baseline(get_connection())
