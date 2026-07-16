"""Subset Explorer tab (Part 4): baseline subset breakdowns + posed answer."""
import plotly.express as px
import streamlit as st

from dashboard.data import (
    load_baseline,
    load_sample_metadata,
    load_final_answer,
)
from dashboard.components import kpi_row
from analysis.subsets import (
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)


def render():
    st.subheader("Baseline subset explorer (time = 0)")
    meta = load_sample_metadata()
    c1, c2 = st.columns(2)
    with c1:
        condition = st.selectbox("Condition",
                                 sorted(meta["condition"].dropna().unique()),
                                 index=_default_index(meta["condition"], "melanoma"),
                                 key="subset_condition")
    with c2:
        treatment = st.selectbox("Treatment",
                                 sorted(meta["treatment"].dropna().unique()),
                                 index=_default_index(meta["treatment"], "miraclib"),
                                 key="subset_treatment")

    subset = load_baseline(condition, treatment, "PBMC")
    if subset.empty:
        st.info("No baseline PBMC samples match this condition/treatment.")
    else:
        kpi_row([("Samples", len(subset)),
                 ("Subjects", subset["subject_id"].nunique())])
        col1, col2, col3 = st.columns(3)
        with col1:
            _bar(samples_per_project(subset), "project", "n_samples",
                 "Samples per project")
        with col2:
            _bar(subjects_by_response(subset), "response", "n_subjects",
                 "Subjects by response")
        with col3:
            _bar(subjects_by_sex(subset), "sex", "n_subjects", "Subjects by sex")

    st.divider()
    answer = load_final_answer()
    st.metric(
        "Avg B cells — melanoma male responders at baseline "
        "(all sample & treatment types)",
        f"{answer:.2f}",
    )


def _default_index(series, value):
    options = sorted(series.dropna().unique())
    return options.index(value) if value in options else 0


def _bar(df, x, y, title):
    fig = px.bar(df, x=x, y=y, title=title)
    st.plotly_chart(fig, use_container_width=True)
