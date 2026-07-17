"""Part 4 · Subset Analysis: baseline breakdowns, the posed answer, and a
degenerate-subset coverage callout."""
import plotly.express as px
import streamlit as st

from dashboard.data import load_baseline, load_final_answer, load_subset_coverage
from dashboard.components import kpi_row
from analysis.subsets import (
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)


def _bar(df, x, y, title):
    st.plotly_chart(px.bar(df, x=x, y=y, title=title), use_container_width=True)


def render():
    st.subheader("Part 4 · Baseline subset analysis")
    st.write("Melanoma PBMC samples at baseline (time = 0) from miraclib-treated "
             "patients, broken down by project, response, and sex.")

    subset = load_baseline("melanoma", "miraclib", "PBMC")
    kpi_row([("Baseline samples", len(subset)),
             ("Subjects", subset["subject_id"].nunique())])

    c1, c2, c3 = st.columns(3)
    with c1:
        _bar(samples_per_project(subset), "project", "n_samples",
             "Samples per project")
    with c2:
        _bar(subjects_by_response(subset), "response", "n_subjects",
             "Subjects by response")
    with c3:
        _bar(subjects_by_sex(subset), "sex", "n_subjects", "Subjects by sex")

    st.divider()
    st.metric(
        "Avg B cells — melanoma male responders at baseline "
        "(all sample & treatment types)",
        f"{load_final_answer():.2f}",
    )

    st.divider()
    st.markdown("#### Data coverage — which subsets are degenerate")
    st.write("Baseline (time = 0) PBMC subset sizes by condition and treatment. "
             "A subset is **degenerate** when it has no responder/non-responder "
             "label — healthy patients and untreated (`none`) arms carry no "
             "response, so they cannot be compared.")
    st.dataframe(load_subset_coverage(), use_container_width=True, hide_index=True)
