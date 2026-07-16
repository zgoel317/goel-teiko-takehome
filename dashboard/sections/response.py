"""Response Analysis tab (Part 3): layered responder vs non-responder stats."""
import streamlit as st

from dashboard.data import (
    load_sample_metadata,
    load_cohort,
    load_comparison,
    load_boxplot_figure,
)
from dashboard.components import finding_summary, style_comparison

_MIN_PER_GROUP = 3


def _default_index(series, value):
    options = sorted(series.dropna().unique())
    return options.index(value) if value in options else 0


def render():
    st.subheader("Responders vs non-responders")
    meta = load_sample_metadata()
    c1, c2, c3 = st.columns(3)
    with c1:
        condition = st.selectbox("Condition",
                                 sorted(meta["condition"].dropna().unique()),
                                 index=_default_index(meta["condition"], "melanoma"),
                                 key="response_condition")
    with c2:
        treatment = st.selectbox("Treatment",
                                 sorted(meta["treatment"].dropna().unique()),
                                 index=_default_index(meta["treatment"], "miraclib"),
                                 key="response_treatment")
    with c3:
        sample_type = st.selectbox("Sample type",
                                   sorted(meta["sample_type"].dropna().unique()),
                                   index=_default_index(meta["sample_type"], "PBMC"),
                                   key="response_sample_type")

    cohort = load_cohort(condition, treatment, sample_type)
    if cohort.empty:
        st.info("No responders/non-responders match this cohort.")
        return
    per_group = cohort.drop_duplicates("subject_id")["response"].value_counts()
    if set(per_group.index) != {"yes", "no"} or per_group.min() < _MIN_PER_GROUP:
        st.warning(
            f"Not enough subjects to compare: need at least {_MIN_PER_GROUP} "
            "responders and non-responders in this cohort.")
        return

    try:
        result = load_comparison(condition, treatment, sample_type)
        fig = load_boxplot_figure(condition, treatment, sample_type)
    except Exception:
        st.warning("This cohort could not be analysed (degenerate or too-small "
                   "groups for the statistical models). Try a broader filter.")
        return

    st.markdown("#### Finding")
    st.info(finding_summary(result))
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("Full statistics table"):
        st.caption("Highlighted rows are significant under the primary "
                   "mixed-effects model (FDR q < 0.05).")
        st.dataframe(style_comparison(result), use_container_width=True)

    with st.expander("Methodology & assumptions"):
        st.markdown(
            "- **Primary — linear mixed-effects model** "
            "(`percentage ~ response + C(time)`, random intercept per subject): "
            "uses all timepoints while respecting repeated measures.\n"
            "- **Cross-check — Mann-Whitney U** on per-subject means: "
            "distribution-free, independent units.\n"
            "- **Multiple comparisons — Benjamini-Hochberg FDR** across the "
            "five populations.\n"
            "- **Effect size — bootstrap 95% CI** of the responder − "
            "non-responder difference.\n"
            "- A population is a **candidate** if significant under the mixed "
            "model, and **confirmed** only if also significant under the "
            "non-parametric test."
        )
